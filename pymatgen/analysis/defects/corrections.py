# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Implementation of defect correction methods.
"""

import logging
import numpy as np
import scipy
import math
from scipy import stats
from pymatgen.core import Structure
from pymatgen.analysis.defects.core import DefectCorrection, Interstitial
from pymatgen.analysis.defects.utils import ang_to_bohr, hart_to_ev, eV_to_k, \
    generate_reciprocal_vectors_squared, QModel, converge, tune_for_gamma, \
    generate_R_and_G_vecs, kumagai_to_V

__author__ = "Danny Broberg, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Mar 15, 2018"

logger = logging.getLogger(__name__)


class FreysoldtCorrection(DefectCorrection):
    """
    A class for FreysoldtCorrection class. Largely adapated from PyCDT code

    If this correction is used, please reference Freysoldt's original paper.
    doi: 10.1103/PhysRevLett.102.016402
    """

    def __init__(self, dielectric_const, q_model=None, energy_cutoff=520, madetol=0.0001, axis=None):
        """
        Initializes the FreysoldtCorrection class
        Args:
            dielectric_const (float or 3x3 matrix): Dielectric constant for the structure
            q_model (QModel): instantiated QModel object or None.
                Uses default parameters to instantiate QModel if None supplied
            energy_cutoff (int): Maximum energy in eV in reciprocal space to perform
                integration for potential correction.
            madeltol(float): Convergence criteria for the Madelung energy for potential correction
            axis (int): Axis to calculate correction.
                If axis is None, then averages over all three axes is performed.
        """
        self.q_model = QModel() if not q_model else q_model
        self.energy_cutoff = energy_cutoff
        self.madetol = madetol
        self.dielectric_const = dielectric_const

        if isinstance(dielectric_const, int) or \
                isinstance(dielectric_const, float):
            self.dielectric = float(dielectric_const)
        else:
            self.dielectric = float(np.mean(np.diag(dielectric_const)))

        self.axis = axis

        self.metadata = {"pot_plot_data": {}, "pot_corr_uncertainty_md": {}}

    def get_correction(self, entry):
        """
        Gets the Freysoldt correction for a defect entry
        Args:
            entry (DefectEntry): defect entry to compute Freysoldt correction on.

                Requires following keys to exist in DefectEntry.parameters dict:

                    axis_grid (3 x NGX where NGX is the length of the NGX grid
                    in the x,y and z axis directions. Same length as planar
                    average lists):
                        A list of 3 numpy arrays which contain the cartesian axis
                        values (in angstroms) that correspond to each planar avg
                        potential supplied.

                    bulk_planar_averages (3 x NGX where NGX is the length of
                    the NGX grid in the x,y and z axis directions.):
                        A list of 3 numpy arrays which contain the planar averaged
                        electrostatic potential for the bulk supercell.

                    defect_planar_averages (3 x NGX where NGX is the length of
                    the NGX grid in the x,y and z axis directions.):
                        A list of 3 numpy arrays which contain the planar averaged
                        electrostatic potential for the defective supercell.

                    initial_defect_structure (Structure) structure corresponding to
                        initial defect supercell structure (uses Lattice for charge correction)

                    defect_frac_sc_coords (3 x 1 array) Fractional co-ordinates of
                        defect location in supercell structure
        Returns:
            FreysoldtCorrection values as a dictionary
        """

        if self.axis is None:
            list_axis_grid = np.array(entry.parameters["axis_grid"])
            list_bulk_plnr_avg_esp = np.array(entry.parameters["bulk_planar_averages"])
            list_defect_plnr_avg_esp = np.array(entry.parameters["defect_planar_averages"])
            list_axes = range(len(list_axis_grid))
        else:
            list_axes = np.array(self.axis)
            list_axis_grid, list_bulk_plnr_avg_esp, list_defect_plnr_avg_esp = [], [], []
            for ax in list_axes:
                list_axis_grid.append(np.array(entry.parameters["axis_grid"][ax]))
                list_bulk_plnr_avg_esp.append(np.array(entry.parameters["bulk_planar_averages"][ax]))
                list_defect_plnr_avg_esp.append(np.array(entry.parameters["defect_planar_averages"][ax]))

        lattice = entry.parameters["initial_defect_structure"].lattice.copy()
        defect_frac_coords = entry.parameters["defect_frac_sc_coords"]

        q = entry.defect.charge

        es_corr = self.perform_es_corr(lattice, entry.charge)

        pot_corr_tracker = []

        for x, pureavg, defavg, axis in zip(list_axis_grid, list_bulk_plnr_avg_esp, list_defect_plnr_avg_esp,
                                            list_axes):
            tmp_pot_corr = self.perform_pot_corr(
                x, pureavg, defavg, lattice, entry.charge, defect_frac_coords,
                axis, widthsample=1.0)
            pot_corr_tracker.append(tmp_pot_corr)

        pot_corr = np.mean(pot_corr_tracker)

        entry.parameters["freysoldt_meta"] = dict(self.metadata)
        entry.parameters["potalign"] = pot_corr / (-q) if q else 0.

        return {"freysoldt_electrostatic": es_corr, "freysoldt_potential_alignment": pot_corr}

    def perform_es_corr(self, lattice, q, step=1e-4):
        """
        Peform Electrostatic Freysoldt Correction
        Args:
            lattice: Pymatgen lattice object
            q (int): Charge of defect
            step (float): step size for numerical integration
        Return:
            Electrostatic Point Charge contribution to Freysoldt Correction (float)
        """
        logger.info("Running Freysoldt 2011 PC calculation (should be " "equivalent to sxdefectalign)")
        logger.debug("defect lattice constants are (in angstroms)" + str(lattice.abc))

        [a1, a2, a3] = ang_to_bohr * np.array(lattice.get_cartesian_coords(1))
        logging.debug("In atomic units, lat consts are (in bohr):" + str([a1, a2, a3]))
        vol = np.dot(a1, np.cross(a2, a3))  # vol in bohr^3

        def e_iso(encut):
            gcut = eV_to_k(encut)  # gcut is in units of 1/A
            return scipy.integrate.quad(lambda g: self.q_model.rho_rec(g * g) ** 2, step, gcut)[0] * (q ** 2) / np.pi

        def e_per(encut):
            eper = 0
            for g2 in generate_reciprocal_vectors_squared(a1, a2, a3, encut):
                eper += (self.q_model.rho_rec(g2) ** 2) / g2
            eper *= (q ** 2) * 2 * round(np.pi, 6) / vol
            eper += (q ** 2) * 4 * round(np.pi, 6) * self.q_model.rho_rec_limit0 / vol
            return eper

        eiso = converge(e_iso, 5, self.madetol, self.energy_cutoff)
        logger.debug("Eisolated : %f", round(eiso, 5))

        eper = converge(e_per, 5, self.madetol, self.energy_cutoff)

        logger.info("Eperiodic : %f hartree", round(eper, 5))
        logger.info("difference (periodic-iso) is %f hartree", round(eper - eiso, 6))
        logger.info("difference in (eV) is %f", round((eper - eiso) * hart_to_ev, 4))

        es_corr = round((eiso - eper) / self.dielectric * hart_to_ev, 6)
        logger.info("Defect Correction without alignment %f (eV): ", es_corr)
        return es_corr

    def perform_pot_corr(self,
                         axis_grid,
                         pureavg,
                         defavg,
                         lattice,
                         q,
                         defect_frac_position,
                         axis,
                         widthsample=1.0):
        """
        For performing planar averaging potential alignment
        Args:
             axis_grid (1 x NGX where NGX is the length of the NGX grid
                    in the axis direction. Same length as pureavg list):
                        A numpy array which contain the cartesian axis
                        values (in angstroms) that correspond to each planar avg
                        potential supplied.
             pureavg (1 x NGX where NGX is the length of the NGX grid in
                    the axis direction.):
                        A numpy array for the planar averaged
                        electrostatic potential of the bulk supercell.
             defavg (1 x NGX where NGX is the length of the NGX grid in
                    the axis direction.):
                        A numpy array for the planar averaged
                        electrostatic potential of the defect supercell.
             lattice: Pymatgen Lattice object of the defect supercell
             q (float or int): charge of the defect
             defect_frac_position: Fracitional Coordinates of the defect in the supercell
             axis (int): axis for performing the freysoldt correction on
             widthsample (float): width (in Angstroms) of the region in between defects
                where the potential alignment correction is averaged. Default is 1 Angstrom.
        Returns:
            Potential Alignment contribution to Freysoldt Correction (float)
        """
        logging.debug("run Freysoldt potential alignment method for axis " + str(axis))
        nx = len(axis_grid)

        # shift these planar averages to have defect at origin
        axfracval = defect_frac_position[axis]
        axbulkval = axfracval * lattice.abc[axis]
        if axbulkval < 0:
            axbulkval += lattice.abc[axis]
        elif axbulkval > lattice.abc[axis]:
            axbulkval -= lattice.abc[axis]

        if axbulkval:
            for i in range(nx):
                if axbulkval < axis_grid[i]:
                    break
            rollind = len(axis_grid) - i
            pureavg = np.roll(pureavg, rollind)
            defavg = np.roll(defavg, rollind)

        # if not self._silence:
        logger.debug("calculating lr part along planar avg axis")
        reci_latt = lattice.reciprocal_lattice
        dg = reci_latt.abc[axis]
        dg /= ang_to_bohr  # convert to bohr to do calculation in atomic units

        # Build background charge potential with defect at origin
        v_G = np.empty(len(axis_grid), np.dtype("c16"))
        v_G[0] = 4 * np.pi * -q / self.dielectric * self.q_model.rho_rec_limit0
        g = np.roll(np.arange(-nx / 2, nx / 2, 1, dtype=int), int(nx / 2)) * dg
        g2 = np.multiply(g, g)[1:]
        v_G[1:] = 4 * np.pi / (self.dielectric * g2) * -q * self.q_model.rho_rec(g2)
        v_G[nx // 2] = 0 if not (nx % 2) else v_G[nx // 2]

        # Get the real space potential by peforming a  fft and grabbing the imaginary portion
        v_R = np.fft.fft(v_G)

        if abs(np.imag(v_R).max()) > self.madetol:
            raise Exception("imaginary part found to be %s", repr(np.imag(v_R).max()))
        v_R /= (lattice.volume * ang_to_bohr ** 3)
        v_R = np.real(v_R) * hart_to_ev

        # get correction
        short = (np.array(defavg) - np.array(pureavg) - np.array(v_R))
        checkdis = int((widthsample / 2) / (axis_grid[1] - axis_grid[0]))
        mid = int(len(short) / 2)

        tmppot = [short[i] for i in range(mid - checkdis, mid + checkdis + 1)]
        logger.debug("shifted defect position on axis (%s) to origin", repr(axbulkval))
        logger.debug("means sampling region is (%f,%f)", axis_grid[mid - checkdis], axis_grid[mid + checkdis])

        C = -np.mean(tmppot)
        logger.debug("C = %f", C)
        final_shift = [short[j] + C for j in range(len(v_R))]
        v_R = [elmnt - C for elmnt in v_R]

        logger.info("C value is averaged to be %f eV ", C)
        logger.info("Potentital alignment energy correction (-q*delta V):  %f (eV)", -q * C)
        self.pot_corr = -q * C

        # log plotting data:
        self.metadata["pot_plot_data"][axis] = {
            "Vr": v_R,
            "x": axis_grid,
            "dft_diff": np.array(defavg) - np.array(pureavg),
            "final_shift": final_shift,
            "check": [mid - checkdis, mid + checkdis + 1]
        }

        # log uncertainty:
        self.metadata["pot_corr_uncertainty_md"][axis] = {"stats": stats.describe(tmppot)._asdict(), "potcorr": -q * C}

        return self.pot_corr

    def plot(self, axis, title=None, saved=False):
        """
        Plots the planar average electrostatic potential against the Long range and
        short range models from Freysoldt. Must run perform_pot_corr or get_correction
        (to load metadata) before this can be used.
        Args:
             axis (int): axis to plot
             title (str): Title to be given to plot. Default is no title.
             saved (bool): Whether to save file or not. If False then returns plot
                object. If True then saves plot as   str(title) + "FreyplnravgPlot.pdf"

        """
        import matplotlib.pyplot as plt
        if not self.metadata["pot_plot_data"]:
            raise ValueError("Cannot plot potential alignment before running correction!")

        x = self.metadata['pot_plot_data'][axis]['x']
        v_R = self.metadata['pot_plot_data'][axis]['Vr']
        dft_diff = self.metadata['pot_plot_data'][axis]['dft_diff']
        final_shift = self.metadata['pot_plot_data'][axis]['final_shift']
        check = self.metadata['pot_plot_data'][axis]['check']

        plt.figure()
        plt.clf()
        plt.plot(x, v_R, c="green", zorder=1, label="long range from model")
        plt.plot(x, dft_diff, c="red", label="DFT locpot diff")
        plt.plot(x, final_shift, c="blue", label="short range (aligned)")

        tmpx = [x[i] for i in range(check[0], check[1])]
        plt.fill_between(tmpx, -100, 100, facecolor="red", alpha=0.15, label="sampling region")

        plt.xlim(round(x[0]), round(x[-1]))
        ymin = min(min(v_R), min(dft_diff), min(final_shift))
        ymax = max(max(v_R), max(dft_diff), max(final_shift))
        plt.ylim(-0.2 + ymin, 0.2 + ymax)
        plt.xlabel(r"distance along axis ($\AA$)", fontsize=15)
        plt.ylabel("Potential (V)", fontsize=15)
        plt.legend(loc=0)
        plt.axhline(y=0, linewidth=0.2, color="black")
        plt.title(str(title) + " defect potential", fontsize=18)
        plt.xlim(0, max(x))
        if saved:
            plt.savefig(str(title) + "FreyplnravgPlot.pdf")
            return
        else:
            return plt


class KumagaiCorrection(DefectCorrection):
    """
    A class for KumagaiCorrection class. Largely adapated from PyCDT code

    If this correction is used, please reference Kumagai and Oba's original paper
    (doi: 10.1103/PhysRevB.89.195205) as well as Freysoldt's original
    paper (doi: 10.1103/PhysRevLett.102.016402)

    NOTE that equations 8 and 9 from Kumagai et al. reference are divided by (4 pi) to get SI units
    """

    def __init__(self,
                 dielectric_tensor,
                 sampling_radius=None,
                 gamma=None):
        """
        Initializes the Kumagai Correction
        Args:
            dielectric_tensor (float or 3x3 matrix): Dielectric constant for the structure

            optional data that can be tuned:
                sampling_radius (float): radius (in Angstrom) which sites must be outside
                    of to be included in the correction. Publication by Kumagai advises to
                    use Wigner-Seitz radius of defect supercell, so this is default value.
                gamma (float): convergence parameter for gamma function.
                    Code will automatically determine this if set to None.
        """
        self.metadata = {"gamma": gamma, "sampling_radius": sampling_radius, "potalign": None}

        if isinstance(dielectric_tensor, int) or \
                isinstance(dielectric_tensor, float):
            self.dielectric = np.identity(3) * dielectric_tensor
        else:
            self.dielectric = np.array(dielectric_tensor)

    def get_correction(self, entry):
        """
        Gets the Kumagai correction for a defect entry
        Args:
            entry (DefectEntry): defect entry to compute Kumagai correction on.

                Requires following parameters in the DefectEntry to exist:

                    bulk_atomic_site_averages (list):  list of bulk structure"s atomic site averaged ESPs * charge,
                        in same order as indices of bulk structure
                        note this is list given by VASP's OUTCAR (so it is multiplied by a test charge of -1)

                    defect_atomic_site_averages (list):  list of defect structure"s atomic site averaged ESPs * charge,
                        in same order as indices of defect structure
                        note this is list given by VASP's OUTCAR (so it is multiplied by a test charge of -1)

                    site_matching_indices (list):  list of corresponding site index values for
                        bulk and defect site structures EXCLUDING the defect site itself
                        (ex. [[bulk structure site index, defect structure"s corresponding site index], ... ]

                    initial_defect_structure (Structure): Pymatgen Structure object representing un-relaxed defect
                        structure

                    defect_frac_sc_coords (array): Defect Position in fractional coordinates of the supercell
                        given in bulk_structure
        Returns:
            KumagaiCorrection values as a dictionary

        """
        bulk_atomic_site_averages = entry.parameters["bulk_atomic_site_averages"]
        defect_atomic_site_averages = entry.parameters["defect_atomic_site_averages"]
        site_matching_indices = entry.parameters["site_matching_indices"]
        defect_sc_structure = entry.parameters["initial_defect_structure"]
        defect_frac_sc_coords = entry.parameters["defect_frac_sc_coords"]

        lattice = defect_sc_structure.lattice
        volume = lattice.volume
        q = entry.defect.charge

        if not self.metadata["gamma"]:
            self.metadata["gamma"] = tune_for_gamma(lattice, self.dielectric)

        prec_set = [25, 28]
        g_vecs, recip_summation, r_vecs, real_summation = generate_R_and_G_vecs(self.metadata["gamma"],
                                                                                prec_set, lattice, self.dielectric)

        pot_shift = self.get_potential_shift(self.metadata["gamma"], volume)
        si = self.get_self_interaction(self.metadata["gamma"])
        es_corr = [(real_summation[ind] + recip_summation[ind] + pot_shift + si) for ind in range(2)]

        # increase precision if correction is not converged yet
        # TODO: allow for larger prec_set to be tried if this fails
        if abs(es_corr[0] - es_corr[1]) > 0.0001:
            logger.debug("Es_corr summation not converged! ({} vs. {})\nTrying a larger prec_set...".format(es_corr[0],
                                                                                                            es_corr[1]))
            prec_set = [30, 35]
            g_vecs, recip_summation, r_vecs, real_summation = generate_R_and_G_vecs(self.metadata["gamma"],
                                                                                    prec_set, lattice, self.dielectric)
            es_corr = [(real_summation[ind] + recip_summation[ind] + pot_shift + si) for ind in range(2)]
            if abs(es_corr[0] - es_corr[1]) < 0.0001:
                raise ValueError("Correction still not converged after trying prec_sets up to 35... serious error.")

        es_corr = es_corr[0] * -(q ** 2.) * kumagai_to_V / 2.  # [eV]

        # if no sampling radius specified for pot align, then assuming Wigner-Seitz radius:
        if not self.metadata["sampling_radius"]:
            wz = lattice.get_wigner_seitz_cell()
            dist = []
            for facet in wz:
                midpt = np.mean(np.array(facet), axis=0)
                dist.append(np.linalg.norm(midpt))
            self.metadata["sampling_radius"] = min(dist)

        # assemble site_list based on matching indices
        # [[defect_site object, Vqb for site], .. repeat for all non defective sites]
        site_list = []
        for bs_ind, ds_ind in site_matching_indices:
            Vqb = -(defect_atomic_site_averages[int(ds_ind)] - bulk_atomic_site_averages[int(bs_ind)])
            site_list.append([defect_sc_structure[int(ds_ind)], Vqb])

        pot_corr = self.perform_pot_corr(defect_sc_structure, defect_frac_sc_coords, site_list,
                                         self.metadata["sampling_radius"], q, r_vecs[0],
                                         g_vecs[0], self.metadata["gamma"])

        entry.parameters["kumagai_meta"] = dict(self.metadata)
        entry.parameters["potalign"] = pot_corr / (-q) if q else 0.

        return {"kumagai_electrostatic": es_corr, "kumagai_potential_alignment": pot_corr}

    def perform_es_corr(self, gamma, prec, lattice, charge):
        """
        Peform Electrostatic Kumagai Correction
        Args:
            gamma (float): Ewald parameter
            prec (int): Precision parameter for reciprical/real lattice vector generation
            lattice: Pymatgen Lattice object corresponding to defect supercell
            charge (int): Defect charge
        Return:
            Electrostatic Point Charge contribution to Kumagai Correction (float)
        """
        volume = lattice.volume

        g_vecs, recip_summation, r_vecs, real_summation = generate_R_and_G_vecs(gamma, [prec], lattice, self.dielectric)
        recip_summation = recip_summation[0]
        real_summation = real_summation[0]

        es_corr = (recip_summation + real_summation + self.get_potential_shift(gamma, volume) +
                   self.get_self_interaction(gamma))

        es_corr *= -(charge ** 2.) * kumagai_to_V / 2.  # [eV]

        return es_corr

    def perform_pot_corr(self, defect_structure, defect_frac_coords, site_list,
                         sampling_radius, q, r_vecs, g_vecs, gamma):
        """
        For performing potential alignment in manner described by Kumagai et al.
        Args:
            defect_structure: Pymatgen Structure object corrsponding to the defect supercell

            defect_frac_coords (array): Defect Position in fractional coordinates of the supercell
                given in bulk_structure

            site_list: list of corresponding site index values for
                bulk and defect site structures EXCLUDING the defect site itself
                (ex. [[bulk structure site index, defect structure"s corresponding site index], ... ]

            sampling_radius (float): radius (in Angstrom) which sites must be outside
                of to be included in the correction. Publication by Kumagai advises to
                use Wigner-Seitz radius of defect supercell, so this is default value.

            q (int): Defect charge

            r_vecs: List of real lattice vectors to use in summation

            g_vecs: List of reciprocal lattice vectors to use in summation

            gamma (float): Ewald parameter

        Return:
            Potential alignment contribution to Kumagai Correction (float)
        """
        volume = defect_structure.lattice.volume
        potential_shift = self.get_potential_shift(gamma, volume)

        pot_dict = {}  # keys will be site index in the defect structure
        for_correction = []  # region to sample for correction

        # for each atom, do the following:
        # (a) get relative_vector from defect_site to site in defect_supercell structure
        # (b) recalculate the recip and real summation values based on this r_vec
        # (c) get information needed for pot align
        for site, Vqb in site_list:
            dist, jimage = site.distance_and_image_from_frac_coords(defect_frac_coords)
            vec_defect_to_site = defect_structure.lattice.get_cartesian_coords(site.frac_coords -
                                                                               jimage - defect_frac_coords)
            dist_to_defect = np.linalg.norm(vec_defect_to_site)
            if abs(dist_to_defect - dist) > 0.001:
                raise ValueError("Error in computing vector to defect")

            relative_real_vectors = [r_vec - vec_defect_to_site for r_vec in r_vecs[:]]

            real_sum = self.get_real_summation(gamma, relative_real_vectors)
            recip_sum = self.get_recip_summation(gamma, g_vecs, volume, r=vec_defect_to_site[:])

            Vpc = (real_sum + recip_sum + potential_shift) * kumagai_to_V * q

            defect_struct_index = defect_structure.index(site)
            pot_dict[defect_struct_index] = {
                "Vpc": Vpc,
                "Vqb": Vqb,
                "dist_to_defect": dist_to_defect
            }

            logger.debug("For atom {}\n\tbulk/defect DFT potential difference = "
                         "{}".format(defect_struct_index, Vqb))
            logger.debug("\tanisotropic model charge: {}".format(Vpc))
            logger.debug("\t\treciprocal part: {}".format(recip_sum * kumagai_to_V * q))
            logger.debug("\t\treal part: {}".format(real_sum * kumagai_to_V * q))
            logger.debug("\t\tself interaction part: {}".format(potential_shift * kumagai_to_V * q))
            logger.debug("\trelative_vector to defect: {}".format(vec_defect_to_site))

            if dist_to_defect > sampling_radius:
                logger.debug("\tdistance to defect is {} which is outside minimum sampling "
                             "radius {}".format(dist_to_defect,
                                                sampling_radius))
                for_correction.append(Vqb - Vpc)
            else:
                logger.debug("\tdistance to defect is {} which is inside minimum sampling "
                             "radius {} (so will not include for correction)"
                             "".format(dist_to_defect, sampling_radius))

        if len(for_correction):
            pot_alignment = np.mean(for_correction)
        else:
            logger.info("No atoms sampled for_correction radius!"
                        " Assigning potential alignment value of 0.")
            pot_alignment = 0.

        self.metadata["potalign"] = pot_alignment
        pot_corr = -q * pot_alignment

        # log uncertainty stats:
        self.metadata["pot_corr_uncertainty_md"] = {
            "stats": stats.describe(for_correction)._asdict(),
            "number_sampled": len(for_correction)
        }
        self.metadata["pot_plot_data"] = pot_dict

        logger.info("Kumagai potential alignment (site averaging): %f", pot_alignment)
        logger.info("Kumagai potential alignment correction energy: %f eV", pot_corr)

        return pot_corr

    def get_real_summation(self, gamma, real_vectors):
        """
        Get real summation term from list of real-space vectors
        """
        real_part = 0
        invepsilon = np.linalg.inv(self.dielectric)
        rd_epsilon = np.sqrt(np.linalg.det(self.dielectric))

        for r_vec in real_vectors:
            if np.linalg.norm(r_vec) > 1e-8:
                loc_res = np.sqrt(np.dot(r_vec, np.dot(invepsilon, r_vec)))
                nmr = scipy.special.erfc(gamma * loc_res)
                real_part += nmr / loc_res

        real_part /= (4 * np.pi * rd_epsilon)

        return real_part

    def get_recip_summation(self, gamma, recip_vectors, volume, r=[0., 0., 0.]):
        """
        Get Reciprocal summation term from list of reciprocal-space vectors
        """
        recip_part = 0

        for g_vec in recip_vectors:
            # dont need to avoid G=0, because it will not be
            # in recip list (if generate_R_and_G_vecs is used)
            Gdotdiel = np.dot(g_vec, np.dot(self.dielectric, g_vec))
            summand = np.exp(-Gdotdiel / (4 * (gamma ** 2))) * np.cos(np.dot(g_vec, r)) / Gdotdiel
            recip_part += summand

        recip_part /= volume

        return recip_part

    def get_self_interaction(self, gamma):
        """
        Args:
            gamma ():

        Returns:
            Self-interaction energy of defect.
        """
        determ = np.linalg.det(self.dielectric)
        return - gamma / (2. * np.pi * np.sqrt(np.pi * determ))

    def get_potential_shift(self, gamma, volume):
        """
        Args:
            gamma (float): Gamma
            volume (float): Volume.

        Returns:
            Potential shift for defect.
        """
        return - 0.25 / (volume * gamma ** 2.)

    def plot(self, title=None, saved=False):
        """
        Plots the AtomicSite electrostatic potential against the Long range and short range models
        from Kumagai and Oba (doi: 10.1103/PhysRevB.89.195205)
        """
        import matplotlib.pyplot as plt
        if "pot_plot_data" not in self.metadata.keys():
            raise ValueError("Cannot plot potential alignment before running correction!")

        sampling_radius = self.metadata["sampling_radius"]
        site_dict = self.metadata["pot_plot_data"]
        potalign = self.metadata["potalign"]

        plt.figure()
        plt.clf()

        distances, sample_region = [], []
        Vqb_list, Vpc_list, diff_list = [], [], []
        for site_ind, site_dict in site_dict.items():
            dist = site_dict["dist_to_defect"]
            distances.append(dist)

            Vqb = site_dict["Vqb"]
            Vpc = site_dict["Vpc"]

            Vqb_list.append(Vqb)
            Vpc_list.append(Vpc)
            diff_list.append(Vqb - Vpc)

            if dist > sampling_radius:
                sample_region.append(Vqb - Vpc)

        plt.plot(distances, Vqb_list,
                 color='r', marker='^', linestyle='None',
                 label='$V_{q/b}$')

        plt.plot(distances, Vpc_list,
                 color='g', marker='o', linestyle='None',
                 label='$V_{pc}$')

        plt.plot(distances, diff_list, color='b', marker='x', linestyle='None',
                 label='$V_{q/b}$ - $V_{pc}$')

        x = np.arange(sampling_radius, max(distances) * 1.05, 0.01)
        y_max = max(max(Vqb_list), max(Vpc_list), max(diff_list)) + .1
        y_min = min(min(Vqb_list), min(Vpc_list), min(diff_list)) - .1
        plt.fill_between(x, y_min, y_max, facecolor='red',
                         alpha=0.15, label='sampling region')
        plt.axhline(y=potalign, linewidth=0.5, color='red',
                    label='pot. align. / -q')

        plt.legend(loc=0)
        plt.axhline(y=0, linewidth=0.2, color='black')

        plt.ylim([y_min, y_max])
        plt.xlim([0, max(distances) * 1.1])

        plt.xlabel(r'Distance from defect ($\AA$)', fontsize=20)
        plt.ylabel('Potential (V)', fontsize=20)
        plt.title(str(title) + " atomic site potential plot", fontsize=20)

        if saved:
            plt.savefig(str(title) + "KumagaiESPavgPlot.pdf")
        else:
            return plt


class BandFillingCorrection(DefectCorrection):
    """
    A class for BandFillingCorrection class. Largely adapted from PyCDT code
    """

    def __init__(self, resolution=0.01):
        """
        Initializes the Bandfilling correction

        Args:
            resolution (float): energy resolution to maintain for gap states
        """
        self.resolution = resolution
        self.metadata = {
            "num_hole_vbm": None,
            "num_elec_cbm": None,
            "potalign": None
        }

    def get_correction(self, entry):
        """
        Gets the BandFilling correction for a defect entry
        Args:
            entry (DefectEntry): defect entry to compute BandFilling correction on.
                Requires following parameters in the DefectEntry to exist:
                    eigenvalues
                        dictionary of defect eigenvalues, as stored in a Vasprun object

                    kpoint_weights (list of floats)
                        kpoint weights corresponding to the dictionary of eigenvalues,
                        as stored in a Vasprun object

                    potalign (float)
                        potential alignment for the defect calculation
                        Only applies to non-zero charge,
                        When using potential alignment correction (freysoldt or kumagai),
                        need to divide by -q

                    cbm (float)
                        CBM of bulk calculation (or band structure calculation of bulk);
                        calculated on same level of theory as the defect
                        (ex. GGA defects -> requires GGA cbm)

                    vbm (float)
                        VBM of bulk calculation (or band structure calculation of bulk);
                        calculated on same level of theory as the defect
                        (ex. GGA defects -> requires GGA vbm)
        Returns:
            Bandfilling Correction value as a dictionary

        """
        eigenvalues = entry.parameters["eigenvalues"]
        kpoint_weights = entry.parameters["kpoint_weights"]
        potalign = entry.parameters["potalign"]
        vbm = entry.parameters["vbm"]
        cbm = entry.parameters["cbm"]

        bf_corr = self.perform_bandfill_corr(eigenvalues, kpoint_weights, potalign, vbm, cbm)

        entry.parameters["bandfilling_meta"] = dict(self.metadata)

        return {"bandfilling_correction": bf_corr}

    def perform_bandfill_corr(self, eigenvalues, kpoint_weights, potalign, vbm, cbm):
        """
        This calculates the band filling correction based on excess of electrons/holes in CB/VB...

        Note that the total free holes and electrons may also be used for a "shallow donor/acceptor"
               correction with specified band shifts:
                +num_elec_cbm * Delta E_CBM (or -num_hole_vbm * Delta E_VBM)
        """
        bf_corr = 0.

        self.metadata["potalign"] = potalign
        self.metadata["num_hole_vbm"] = 0.
        self.metadata["num_elec_cbm"] = 0.

        core_occupation_value = list(eigenvalues.values())[0][0][0][1]  # get occupation of a core eigenvalue
        if len(eigenvalues.keys()) == 1:
            # needed because occupation of non-spin calcs is sometimes still 1... should be 2
            spinfctr = 2. if core_occupation_value == 1. else 1.
        elif len(eigenvalues.keys()) == 2:
            spinfctr = 1.
        else:
            raise ValueError("Eigenvalue keys greater than 2")

        # for tracking mid gap states...
        shifted_cbm = potalign + cbm  # shift cbm with potential alignment
        shifted_vbm = potalign + vbm  # shift vbm with potential alignment

        for spinset in eigenvalues.values():
            for kptset, weight in zip(spinset, kpoint_weights):
                for eig, occu in kptset:  # eig is eigenvalue and occu is occupation
                    if (occu and (eig > shifted_cbm - self.resolution)):  # donor MB correction
                        bf_corr += weight * spinfctr * occu * (eig - shifted_cbm)  # "move the electrons down"
                        self.metadata["num_elec_cbm"] += weight * spinfctr * occu
                    elif (occu != core_occupation_value) and (
                            eig <= shifted_vbm + self.resolution):  # acceptor MB correction
                        bf_corr += weight * spinfctr * (core_occupation_value - occu) * (
                                shifted_vbm - eig)  # "move the holes up"
                        self.metadata["num_hole_vbm"] += weight * spinfctr * (core_occupation_value - occu)

        bf_corr *= -1  # need to take negative of this shift for energetic correction

        return bf_corr


class BandEdgeShiftingCorrection(DefectCorrection):
    """
    A class for BandEdgeShiftingCorrection class. Largely adapted from PyCDT code
    """

    def __init__(self):
        """
        Initializes the BandEdgeShiftingCorrection class
        """
        self.metadata = {
            "vbmshift": 0.,
            "cbmshift": 0.,
        }

    def get_correction(self, entry):
        """
        Gets the BandEdge correction for a defect entry
        Args:
            entry (DefectEntry): defect entry to compute BandFilling correction on.
                Requires some parameters in the DefectEntry to properly function:
                    hybrid_cbm (float)
                        CBM of HYBRID bulk calculation one wishes to shift to

                    hybrid_vbm (float)
                        VBM of HYBRID bulk calculation one wishes to shift to

                    cbm (float)
                        CBM of bulk calculation (or band structure calculation of bulk);
                        calculated on same level of theory as the defect
                        (ex. GGA defects -> requires GGA cbm)

                    vbm (float)
                        VBM of bulk calculation (or band structure calculation of bulk);
                        calculated on same level of theory as the defect
                        (ex. GGA defects -> requires GGA vbm)
        Returns:
            BandfillingCorrection value as a dictionary
        """
        hybrid_cbm = entry.parameters["hybrid_cbm"]
        hybrid_vbm = entry.parameters["hybrid_vbm"]
        vbm = entry.parameters["vbm"]
        cbm = entry.parameters["cbm"]

        self.metadata["vbmshift"] = hybrid_vbm - vbm  # note vbmshift has UPWARD as positive convention
        self.metadata["cbmshift"] = hybrid_cbm - cbm  # note cbmshift has UPWARD as positive convention

        charge = entry.charge
        bandedgeshifting_correction = charge * self.metadata["vbmshift"]
        entry.parameters["bandshift_meta"] = dict(self.metadata)

        return {
            "bandedgeshifting_correction": bandedgeshifting_correction
        }


class LevelShiftCorrection(DefectCorrection):
    """
    A class for DefectLevelShiftingCorrection class. Contains several different
    Level shifting approaches that can be tried:
        'iA' = shift OCCUPIED band edge states 100% with the band edges (only shift occupied states, not holes)
             --> update on 3/8. Not going to shift shallow states because of bad resolution...
        'iB' = shift band edge states 100% with the band edges (do hole shifting approach)
        'ii' = Use pawpyseed % to shift KS levels in gap + use approach iA (TODO = iB could also be tried)
        'iii' = If KS level is localized then dont shift it, otherwise do (ii)
        'iv' = If localized level is found in host and there are occupied states
                    near band edges (like iA) then shift to that level. Otherwise shift 100% with band edges
        'v' = Apply Pawpyseed shift to any midgap states based on Procar percentage on nearest neighbors

    Requires some parameters in the DefectEntry to properly function:
        vbm

        cbm

        ["bandshift_meta"]["cbmshift"]
            Shift to apply to GGA CBM to get to CBM of HYBRID bulk calculation

        ["bandshift_meta"]["vbmshift"]
            Shift to apply to GGA VBM to get to VBM of HYBRID bulk calculation

        num_hole_vbm
            number of free holes that were found in valence band for the defect calculation
            calculated in the metadata of the BandFilling Correction

        num_elec_cbm
            number of free electrons that were found in the conduction band for the defect calculation
            calculated in the metadata of the BandFilling Correction

        potalign

        eigenvalues

        kpoint_weights

        ["defect_ks_delocal_data"]

    """

    def __init__(self, tolerance = 0.01, tol_shallow_defect=.15):
        self.tolerance = tolerance
        self.tol_shallow_defect = tol_shallow_defect
        self.metadata = {}

    def get_correction(self, entry, shift_approach='iA'):
        """
        Gets the ShallowLevel correction for a defect entry
        """
        shift_approach = shift_approach.lower()
        vbm = entry.parameters["vbm"]
        cbm = entry.parameters["cbm"]
        vbmshift = entry.parameters["bandshift_meta"]["vbmshift"]
        cbmshift = entry.parameters["bandshift_meta"]["cbmshift"]

        num_hole_vbm = entry.parameters["num_hole_vbm"]
        num_elec_cbm = entry.parameters["num_elec_cbm"]

        potalign = entry.parameters["potalign"]

        eigenvalues = entry.parameters["eigenvalues"]
        kpoint_weights = entry.parameters["kpoint_weights"]

        core_occupation_value = list(eigenvalues.values())[0][0][0][1]  # get occupation of a core eigenvalue
        if len(eigenvalues.keys()) == 1:
            # needed because occupation of non-spin calcs is sometimes still 1... should be 2
            spinfctr = 2. if core_occupation_value == 1. else 1.
        elif len(eigenvalues.keys()) == 2:
            spinfctr = 1.
        else:
            raise ValueError("Eigenvalue keys greater than 2")

        # for tracking mid gap states...
        shifted_cbm = potalign + cbm  # shift cbm with potential alignment
        shifted_vbm = potalign + vbm  # shift vbm with potential alignment

        # grab KS eigenvalue energy and occupation information
        num_occupied_acceptor = 0.
        occupied_KS_deep_states = [] #grab KS band indices which are occupied and in gap and not "shallow acceptor" type
        all_possible_KS_defect_states = [] #grab all possible defect states (shallow states included)
        for spinset in eigenvalues.values():
            for kptset, weight in zip(spinset, kpoint_weights):
                for bandind, bandset in enumerate(kptset):
                    eig, occu = bandset # eig is eigenvalue and occu is occupation
                    if occu and (eig > shifted_vbm + self.tolerance) and (eig < shifted_cbm - self.tolerance):
                        all_possible_KS_defect_states.append( bandind)
                        if (eig <= shifted_vbm + self.tolerance + self.tol_shallow_defect):
                            num_occupied_acceptor += weight * spinfctr * occu
                        elif (eig < shifted_cbm - self.tolerance - self.tol_shallow_defect):
                            occupied_KS_deep_states.append( bandind)

        occupied_KS_deep_states = list(set(occupied_KS_deep_states))
        all_possible_KS_defect_states = list(set(all_possible_KS_defect_states))

        # grab KS wf localization information
        lbl_dict = entry.parameters['defect_ks_delocal_data']['localized_band_indices']
        localized_band_indices = list(set([int(band_index) for spin_list in lbl_dict.values()
                                           for band_index in spin_list]))

        # If occupied KS levels exist, grab eigenvalues and occupation for their band indices
        # also grab eigenvalue information for any possible localized band indices
        KS_defect_energy_and_occu_dict = {}
        local_band_energy_dict = {}
        all_band_eigen_and_occu = {}
        if len(all_possible_KS_defect_states):
            for spinset in eigenvalues.values():
                for kptset, weight in zip(spinset, kpoint_weights):
                    for bandind, bandset in enumerate(kptset):
                        eig, occu = bandset  # eig is eigenvalue and occu is occupation
                        if bandind in all_band_eigen_and_occu.keys():
                            all_band_eigen_and_occu[bandind].append( [eig, weight, weight * spinfctr * occu])
                        else:
                            all_band_eigen_and_occu[bandind] = [[eig, weight, weight * spinfctr * occu]]

                        if bandind in all_possible_KS_defect_states:
                            # store list of [eigenvalue, kpt weight, weighted occupation]
                            if bandind in KS_defect_energy_and_occu_dict.keys():
                                KS_defect_energy_and_occu_dict[bandind].append( [eig, weight, weight * spinfctr * occu])
                            else:
                                KS_defect_energy_and_occu_dict[bandind] = [[eig, weight, weight * spinfctr * occu]]

                        if bandind in localized_band_indices:
                            # store list of [eigenvalue, kpt weight]
                            if bandind in local_band_energy_dict.keys():
                                local_band_energy_dict[bandind].append( [eig, weight])
                            else:
                                local_band_energy_dict[bandind] = [[eig, weight]]

        # grab Pawpyseed information -> is this needed?? Or the info already exists??
        pawpy_band_proj_md = {} #key is band index, value is dict with keys
        #                           {'perc_vbm_shift', 'perc_cbm_shift', 'wgted_eigen', 'tot_occu', 'shift_corr'}
        # for bandind, eig_occu_list in KS_defect_energy_and_occu_dict.items():
        # for bandind in entry.parameters['defect_ks_delocal_data']['pawpyseed'].keys(): #TODO -> uncomment these two lines if running pawpyseed
        #     # TODO -> (3/1) would be cool to have pawpyseed only applied if it stays in gap after the application? Means I can consider levels outside of gap?
        #     eig_occu_list = all_band_eigen_and_occu[bandind]
        for bandind, tmp_eig_occu_list in all_band_eigen_and_occu.items():
            if 'pawpyseed' in entry.parameters['defect_ks_delocal_data'].keys():
                if str(bandind) not in entry.parameters['defect_ks_delocal_data']['pawpyseed'].keys():
                    continue
                else:
                    eig_occu_list = all_band_eigen_and_occu[bandind]
                    pawpy = entry.parameters['defect_ks_delocal_data']['pawpyseed'][str(bandind)]
                    pawpy_band_proj = [np.array( pawpy[0]).mean(), np.array( pawpy[1]).mean()] #just taking average of two spins
            else:
                eig_occu_list = tmp_eig_occu_list[:]
                pawpy_band_proj = None

            eig_vec = np.array( eig_occu_list)[:, 0]
            wgt_vec = np.array( eig_occu_list)[:, 1]
            wgt_vec /= wgt_vec.sum()
            wgt_avg_eigen = np.dot( eig_vec, wgt_vec)
            occu  =  np.array( eig_occu_list)[:, 2].sum()
            # print("\t{}: wgt_avg_eig = {} , occu = {}".format( bandind, wgt_avg_eigen, occu))

            if pawpy_band_proj is None:
                # a FAKE pawpy parser
                print("WARNING: using fake pawpybandproj...")
                if wgt_avg_eigen <= shifted_vbm + entry.parameters['gap'] / 2.:
                    pawpy_band_proj = [1., 0.]
                else:
                    pawpy_band_proj = [0., 1.]

            wgted_shift = (pawpy_band_proj[0] * vbmshift + pawpy_band_proj[1] * cbmshift)
            this_corr = wgted_shift * occu

            pawpy_band_proj_md.update({bandind: {'perc_vbm_shift': pawpy_band_proj[0],
                                                'perc_cbm_shift': pawpy_band_proj[1],
                                                'wgted_eigen': wgt_avg_eigen,
                                                'wgted_shift': wgted_shift, 'tot_occu': occu, 'shift_corr': this_corr}
                                       })


        total_shift_corr = {}
        corr_type, AB_type = shift_approach[:-1], shift_approach[-1]
        if AB_type not in ["a", "b"]:
            corr_type = shift_approach
        elif AB_type == 'a':
            # if corr_type in ['iv', 'v']:
            #     raise ValueError("SHOULD NOT include free carrier shifts at edge?")
            # occu_acceptor_shift = num_occupied_acceptor * vbmshift
            elec_cbm_shift_correction = num_elec_cbm * cbmshift
            total_shift_corr.update( {
                                       # "occu_acceptor_shift": occu_acceptor_shift,
                                       "elec_cbm_shift_correction": elec_cbm_shift_correction
                                       })
            self.metadata.update( {'num_occupied_acceptor': num_occupied_acceptor})
        elif AB_type == 'b':
            # if corr_type in ['iv', 'v']:
            #     raise ValueError("SHOULD NOT include free carrier shifts at edge?")
            hole_vbm_shift_correction = -1. * num_hole_vbm * vbmshift # negative sign because these are holes
            elec_cbm_shift_correction = num_elec_cbm * cbmshift
            total_shift_corr.update( { "hole_vbm_shift_correction": hole_vbm_shift_correction,
                                       "elec_cbm_shift_correction": elec_cbm_shift_correction
                                       } )

        if corr_type == 'ii':
            # print("Running type ii level shift correction (always do Pawpyseed):")
            tot_corr = 0.
            for bandind in occupied_KS_deep_states: #loop over deep (non shallow state) indices
                vdict = pawpy_band_proj_md[bandind]
                # print("\t{}: eigen (rel to vbm) = {}, occu = {}, corr = {}".format( bandind, vdict["wgted_eigen"]-shifted_vbm,
                #                                                                     vdict["tot_occu"],  vdict["shift_corr"]))
                tot_corr += vdict["shift_corr"]

            # print("\nFinal Pawpyseed correction = {}".format( tot_corr))

            self.metadata.update( {'pawpy_band_proj_md': pawpy_band_proj_md.copy()})
            total_shift_corr.update( { "pawpyseed_always_KS_shift": tot_corr})

        elif corr_type == 'iii':
            # print("Running type iii level shift correction (only do Pawpyseed when not detected as localized):")
            tot_corr = 0.
            for bandind in occupied_KS_deep_states:
                vdict = pawpy_band_proj_md[bandind]
                # print("\t{}: eigen (rel to vbm) = {}, occu = {}, corr = {}".format( bandind, vdict["wgted_eigen"]-shifted_vbm,
                #                                                                     vdict["tot_occu"],  vdict["shift_corr"]))
                # if bandind in localized_band_indices:
                #     print("\tLOCALIZED  LEVEL! will not move.")
                # else:
                #     tot_corr += vdict["shift_corr"]
                if bandind not in localized_band_indices:
                    tot_corr += vdict["shift_corr"]

            # print("\nFinal Pawpyseed correction = {}".format( tot_corr))

            self.metadata.update( {'pawpy_band_proj_md': pawpy_band_proj_md.copy(),
                                   'localized_levels': localized_band_indices[:]})
            total_shift_corr.update( { "pawpyseed_nonlocal_KS_shift": tot_corr})

        elif corr_type == 'iv':
            # print("Running type iv level shift correction (if non-localized levels are not occupied, then try to find localized level "
            #       "for it; if this isnt doable then do Pawpyseed):")
            # first create wgt energy eigen value of all localized levels...
            wgt_avg_eigen_for_each_local_band = {}
            for bandind, eig_occu_list in local_band_energy_dict.items():
                eig_vec = np.array( eig_occu_list)[:, 0]
                wgt_vec = np.array( eig_occu_list)[:, 1]
                wgt_vec /= wgt_vec.sum()
                wgt_avg_eigen = np.dot( eig_vec, wgt_vec)
                wgt_avg_eigen_for_each_local_band[bandind] = wgt_avg_eigen

            #next consider all things that could possibly be shifted:
            # (i) occupied KS levels which are NOT already localized,
            # (ii) occupied acceptor/donors states which are NOT already localized
            track_localized_mapping = wgt_avg_eigen_for_each_local_band.copy() #so I can delete things

            smart_shift_md = {'all_possible_KS_defect_states': all_possible_KS_defect_states[:],
                              'track_for_corr': {}}
            tot_corr = 0.
            needs_shifting = []
            for bandind in all_possible_KS_defect_states:
                eigval = round(pawpy_band_proj_md[bandind]["wgted_eigen"]-shifted_vbm,2)
                if bandind in track_localized_mapping.keys():
                    smart_shift_md['track_for_corr'].update( {bandind: {'type': 'localized', 'corr': 0.}})
                    del track_localized_mapping[bandind] # dont want to recount this level later
                    # print("\t{} ({} above vbm) is localized! No correction needed.".format( bandind, eigval))
                else:
                    # print("\t{} ({} above vbm) needs shifting!".format( bandind, eigval))
                    needs_shifting.append( bandind)

            find_smallest_shifts = []
            for bandind in needs_shifting:
                bands_eigen = pawpy_band_proj_md[bandind]['wgted_eigen']
                find_smallest_shifts.extend( [[abs(val - bands_eigen), tb, bandind] for tb, val in track_localized_mapping.items()])

            find_smallest_shifts.sort()
            for diff, tb, bandind in find_smallest_shifts:
                if (bandind in smart_shift_md['track_for_corr'].keys()):
                    continue
                elif diff > max(abs(vbmshift), abs(cbmshift)):  # if shift is bigger than the host band shifts then use pawpyshift
                    this_corr = pawpy_band_proj_md[bandind]["shift_corr"]
                    smart_shift_md['track_for_corr'].update( {bandind: {'type': 'pawpyseed', 'corr': this_corr}})
                    tot_corr += this_corr
                    # print("\t{} could not be matched to localized level. Using pawpyseed.".format( bandind))
                else:
                    shifting_to = pawpy_band_proj_md[tb]['wgted_eigen'] - pawpy_band_proj_md[bandind]['wgted_eigen']
                    occu = pawpy_band_proj_md[bandind]['tot_occu']
                    this_corr = shifting_to * occu
                    smart_shift_md['track_for_corr'].update( {bandind: {'type': 'shift_to_{}'.format( tb),
                                                                        'corr': this_corr}})
                    # print("\t{} matched to localized level {} away. Using this for shift.".format( bandind,
                    #                                                                                round(pawpy_band_proj_md[tb]['wgted_eigen'],2)))

            #if ran out of localized levels, run through and confirm that all possible ks states matched
            for bandind in all_possible_KS_defect_states:
                if bandind not in smart_shift_md['track_for_corr'].keys():
                    this_corr = pawpy_band_proj_md[bandind]["shift_corr"]
                    smart_shift_md['track_for_corr'].update( {bandind: {'type': 'pawpyseed', 'corr': this_corr}})
                    tot_corr += this_corr
                    print("\t{} could not be matched to localized level. Using pawpyseed.".format( bandind))

            if len( smart_shift_md['track_for_corr'].keys()) != len( all_possible_KS_defect_states):
                raise ValueError("Something went wrong with level matching routine?")

            # print("\nFinal Smart KS shift correction = {}".format( tot_corr))

            self.metadata.update( {'pawpy_band_proj_md': pawpy_band_proj_md.copy(),
                                   'smart_shift_md': smart_shift_md.copy(),
                                   'wgt_avg_eigen_for_each_local_band': wgt_avg_eigen_for_each_local_band.copy()})

            total_shift_corr.update( { "smart_KS_shift": tot_corr})

        elif corr_type == 'v':
            """
            i) Figure out a radius which is half way between first and second nearest neighbors 
            to defect (to within some tolerance)
            """
            if isinstance(entry.defect, Interstitial): #need to get defect site. If interstitial, it doesnt exist in bulk
                bs = entry.parameters['initial_defect_structure']
            else:
                bs = entry.parameters['bulk_sc_structure']
            if isinstance(bs, dict):
                bs = Structure.from_dict(bs)
            index = entry.parameters['defect_index_sc_coords']
            def_site = bs.sites[int(index)]
            vals = [[v[1], v[2]] for v in bs.get_neighbors( def_site, 5., include_index=True)]
            vals.sort()
            #get 1nn and 2nn distances (with 0.1 tolerance)
            range_nn = [vals[0][0]]
            for dist,_ in vals:
                if abs(dist - range_nn[0]) > 0.1 and len(range_nn) == 1:
                    range_nn.append( dist)
            halfway_1nn_2nn_rad = np.mean( range_nn) #radius is half way between two
            defect_site_indices_1nn = [] #list of defect structure site indices to consider in procar analysis
            bi_to_di_site_map = {int(bi): int(di) for bi,di in entry.parameters['site_matching_indices']}
            for dist, bindex in vals:
                if dist < halfway_1nn_2nn_rad:
                    if isinstance(entry.defect, Interstitial):
                        defect_site_indices_1nn.append( int(bindex))
                    else:
                        defect_site_indices_1nn.append( bi_to_di_site_map[int(bindex)])

            """
            ii) From Procar, figure out what fraction of the single particle level exists 
            on each atomic site within that radius
            """
            local_multiply_dict = {} #key is band index, value is multiplier = (1 - local percentage on nns)
            spd = entry.parameters['defect_ks_delocal_data']['stored_procar_dat']
            for bandind in all_possible_KS_defect_states:
                try:
                    tot_procar_occu = 0.
                    nn_procar_occu = 0.
                    for spinind in spd.keys():
                        for kptdict in spd[spinind]:
                            for site_ind, vallist in enumerate(kptdict[str(bandind)]['rad_dat']):
                                dist, occu = vallist
                                tot_procar_occu += occu
                                if dist < halfway_1nn_2nn_rad:
                                    nn_procar_occu += occu
                    local_multiply_dict[bandind] = (1- nn_procar_occu/tot_procar_occu)
                except:
                    #if try fails then likely because bandindex is not in procar set... which means we wont use it
                    local_multiply_dict[bandind] = 0.
                # print("{} has perc = {}".format(bandind, nn_procar_occu/tot_procar_occu))

            """
            iii) multiply the pawpyseed shift for that single particle level by (1 - the fraction from (ii) )
            """
            tot_corr = 0.
            for bandind in all_possible_KS_defect_states:
                tot_corr += local_multiply_dict[bandind] * pawpy_band_proj_md[bandind]['shift_corr']

            self.metadata.update( {'pawpy_band_proj_md': pawpy_band_proj_md.copy(),
                                   'local_multiply_dict': local_multiply_dict.copy()})

            total_shift_corr.update( { "nn_adjusted_KS_shift": tot_corr})



        elif corr_type != 'i':
            raise ValueError('Shift_approach = {} not recognized...'.format( shift_approach))

        entry.parameters["levelshift_meta"] = dict(self.metadata)

        return total_shift_corr
