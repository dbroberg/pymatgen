# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

from monty.json import MSONable
from pymatgen.analysis.defects.corrections import FreysoldtCorrection, KumagaiCorrection, \
    BandFillingCorrection, generate_g_sum, find_optimal_gamma, BandEdgeShiftingCorrection
from pymatgen.analysis.defects.core import Vacancy
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
"""
This module implements DefectCompatibility analysis for consideration of
defects
"""

__author__ = "Danny Broberg, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Mar 15, 2018"


class DefectCompatibility(MSONable):
    """
    The DefectCompatibility class combines a list of DefectEntries for a
    given system and applies corrections / suggests failed jobs that
    should not be considered
    Args:

        defect_entries: List of defect_entries to consider.
        user_defect_qualifiers: a dictionary for specifying the dictionary of qualifiers
                                for corrections and delocalization analysis.
                                Defaults are in dictionary above

    required settings for defect_entry.parameters:
        freysoldt: ["axis_grid", "bulk_planar_averages", "defect_planar_averages", "dielectric"]
        kumagai: ["dim", "bulk_atomic_site_averages", "defect_atomic_site_averages", "site_matching_indices",
                    "dielectric]
        bandfilling: ["eigenvalues", "kpoint_weights", "potalign", "vbm", "cbm"]
        bandshifting: ["hybrid_cbm", "hybrid_vbm", "num_hole_vbm", "num_elec_cbm", "vbm", "cbm"]
        defect relaxation/structure analysis: ["final_defect_structure"]

    """

    def __init__(self,
                 plnr_avg_var_tol=0.1,
                 plnr_avg_minmax_tol=0.1,
                 atomic_site_var_tol=0.1,
                 atomic_site_minmax_tol=0.1,
                 tot_relax_tol=0.1,
                 perc_relax_tol=0.1,
                 defect_tot_relax_tol=0.1,
                 preferred_cc='freysoldt',
                 free_chg_cutoff=2.,
                 use_bandfilling=True,
                 use_bandedgeshift=True):
        # TODO: fine tune qualifiers a bit more...
        self.plnr_avg_var_tol = plnr_avg_var_tol
        self.plnr_avg_minmax_tol = plnr_avg_minmax_tol
        self.atomic_site_var_tol = atomic_site_var_tol
        self.atomic_site_minmax_tol = atomic_site_minmax_tol
        self.tot_relax_tol = tot_relax_tol
        self.perc_relax_tol = perc_relax_tol
        self.defect_tot_relax_tol = defect_tot_relax_tol

        self.preferred_cc = preferred_cc
        self.free_chg_cutoff = free_chg_cutoff
        self.use_bandfilling = use_bandfilling
        self.use_bandedgeshift = use_bandedgeshift

    def process_entry(self, defect_entry):
        """
        Process a given Defect entry with chosen qualifiers.
        Performed Defect Corrections as needed and consider delocalization analysis,
            then update corrections to defect entry

        Corrections are applied based on:
            i) if free charges are more than free_chg_cutoff then will not apply charge correction,
                because it no longer is applicable
            ii) use charge correction set by preferred_cc
            iii) only use BandFilling correction if use_bandfilling is set to True
            iv) only use BandEdgeShift correction if use_bandedgeshift is set to True
        """
        self.perform_corrections(defect_entry)

        self.delocalization_analysis(defect_entry)

        corrections = {}
        if (self.free_chg_cutoff < defect_entry.parameters["num_hole_vbm"]) or (
                self.free_chg_cutoff < defect_entry.parameters["num_elec_cbm"]):
            print('Will not use charge correction because too many free charges')
            # TODO: should the potential alignment correction still be used in this scenario, though?
        elif 'freysoldt' in self.preferred_cc.lower():
            frey_meta = defect_entry.parameters['freysoldt_meta']
            frey_corr = frey_meta["freysoldt_electrostatic"] + frey_meta["freysoldt_potential_alignment_correction"]
            corrections.update({'charge_correction': frey_corr})
        else:
            kumagai_meta = defect_entry.parameters['kumagai_meta']
            kumagai_corr = kumagai_meta["kumagai_electrostatic"] + \
                kumagai_meta["kumagai_potential_alignment_correction"]
            corrections.update({'charge_correction': kumagai_corr})

        if self.use_bandfilling:
            bfc_corr = defect_entry.parameters["bandfilling_meta"]["bandfilling_correction"]
            corrections.update({'bandfilling_correction': bfc_corr})

        if self.use_bandedgeshift:
            bandfill_meta = defect_entry.parameters["bandshift_meta"]
            bes_corr = bandfill_meta["vbm_shift_correction"] + bandfill_meta["hole_vbm_shift_correction"] + \
                bandfill_meta["elec_cbm_shift_correction"]
            corrections.update({'bandedgeshifting_correction': bes_corr})

            # also want to update relevant data for phase diagram
            defect_entry.parameters.update({
                'phasediagram_meta': {
                    'vbm': defect_entry.parameters['hybrid_vbm'],
                    'gap': defect_entry.parameters['hybrid_cbm'] - defect_entry.parameters['hybrid_vbm']
                }
            })
        else:  # if not using bandedge shift -> still want to have vbm and gap ready for phase diagram
            defect_entry.parameters.update({
                'phasediagram_meta': {
                    'vbm': defect_entry.parameters['vbm'],
                    'gap': defect_entry.parameters['cbm'] - defect_entry.parameters['vbm']
                }
            })

        defect_entry.correction.update(corrections)

        return defect_entry


    def perform_all_corrections(self, defect_entry):

        # consider running freysoldt correction
        required_frey_params = ["axis_grid", "bulk_planar_averages", "defect_planar_averages", "dielectric"]
        run_freysoldt = True if len( set(defect_entry.parameters.keys()).intersection(required_frey_params)) \
                                == len(required_frey_params) else False
        if not run_freysoldt:
            print('Insufficient DefectEntry parameters exist for Freysoldt Correction.')
        elif 'freysoldt_meta' not in defect_entry.parameters.keys():
            defect_entry = self.run_freysoldt( defect_entry)


        # consider running kumagai correction
        required_kumagai_params = ["dim", "bulk_atomic_site_averages", "defect_atomic_site_averages",
                                   "site_matching_indices", "dielectric"]
        run_kumagai = True if len( set(defect_entry.parameters.keys()).intersection(required_kumagai_params)) \
                                == len(required_kumagai_params) else False
        if not run_kumagai:
            print('Insufficient DefectEntry parameters exist for Kumagai Correction.')
        elif 'kumagai_meta' not in defect_entry.parameters.keys():
            defect_entry = self.run_kumagai( defect_entry)

        # consider running band filling correction
        required_bandfilling_params = ["eigenvalues", "kpoint_weights", "potalign", "vbm", "cbm"]
        run_bandfilling = True if len( set(defect_entry.parameters.keys()).intersection(required_bandfilling_params)) \
                                == len(required_bandfilling_params) else False
        if not run_bandfilling:
            print('Insufficient DefectEntry parameters exist for BandFilling Correction.')
        elif 'bandfilling_meta' not in defect_entry.parameters.keys():
            # TODO: add ability to modify the potalign value to prefer kumagai or freysoldt?
            defect_entry = self.run_bandfilling( defect_entry)

        # consider running band edge shifting correction
        required_bandedge_shifting_params = ["hybrid_cbm", "hybrid_vbm", "num_hole_vbm", "num_elec_cbm", "vbm", "cbm"]
        run_bandedge_shifting = True if len( set(defect_entry.parameters.keys()).intersection(required_bandedge_shifting_params)) \
                                == len(required_bandedge_shifting_params) else False
        if not run_bandedge_shifting:
            print('Insufficient DefectEntry parameters exist for BandShifting Correction.')
        elif 'bandedgeshift_meta' not in defect_entry.parameters.keys():
            defect_entry = self.run_band_edge_shifting( defect_entry)

        return defect_entry

    def perform_freysoldt(self, defect_entry):
        FC = FreysoldtCorrection(defect_entry.parameters['dielectric'])
        freycorr = FC.get_correction(defect_entry)

        freysoldt_meta = FC.metadata.copy()
        freysoldt_meta["freysoldt_potalign"] = defect_entry.parameters["potalign"]
        freysoldt_meta["freysoldt_electrostatic"] = freycorr["freysoldt_electrostatic"]
        freysoldt_meta["freysoldt_potential_alignment_correction"] = freycorr["freysoldt_potential_alignment"]
        defect_entry.parameters.update({'freysoldt_meta': freysoldt_meta})
        return defect_entry

    def perform_kumagai(self, defect_entry):
        # can save alot of time if gamma or g_sum in defect_entry.parameters, so check if they exist
        gamma = defect_entry.parameters['gamma'] if 'gamma' in defect_entry.parameters.keys() else None
        g_sum = defect_entry.parameters['g_sum'] if 'g_sum' in defect_entry.parameters.keys() else None

        if not gamma:
            defect_struct_sc = defect_entry.defect_sc_structure.copy()
            gamma = find_optimal_gamma(defect_struct_sc.lattice, defect_entry.parameters["dielectric"])

        if not g_sum:
            defect_struct_sc = defect_entry.defect_sc_structure.copy()
            g_sum = generate_g_sum(defect_struct_sc.lattice, defect_entry.parameters["dielectric"],
                                   defect_entry.parameters['dim'], gamma)

        KC = KumagaiCorrection(defect_entry.parameters['dielectric'], gamma=gamma, g_sum=g_sum)
        kumagaicorr = KC.get_correction(defect_entry)

        kumagai_meta = {k: v for k, v in KC.metadata.items() if k != 'g_sum'}
        kumagai_meta["kumagai_potalign"] = defect_entry.parameters["potalign"]
        kumagai_meta["kumagai_electrostatic"] = kumagaicorr["kumagai_electrostatic"]
        kumagai_meta["kumagai_potential_alignment_correction"] = kumagaicorr["kumagai_potential_alignment"]
        defect_entry.parameters.update({'kumagai_meta': kumagai_meta})
        return defect_entry

    def run_bandfilling(self, defect_entry):
        BFC = BandFillingCorrection()
        bfc_dict = BFC.get_correction(defect_entry)

        bandfilling_meta = defect_entry.parameters.copy()
        bandfilling_meta["bandfilling_correction"] = bfc_dict['bandfilling']
        defect_entry.parameters.update({'bandfilling_meta': bandfilling_meta,
                                        # also update free holes and electrons for band edge shifting correction...
                                        'num_hole_vbm': bandfilling_meta["num_hole_vbm"],
                                        'num_elec_cbm': bandfilling_meta["num_elec_cbm"]})
        return defect_entry

    def run_band_edge_shifting(self, defect_entry):
        BEC = BandEdgeShiftingCorrection()
        bec_dict = BEC.get_correction(defect_entry)

        bandshift_meta = BEC.metadata.copy()
        bandshift_meta.update(bec_dict)

        defect_entry.parameters.update({"bandedgeshift_meta": bandshift_meta})

        return defect_entry


    def delocalization_analysis(self, defect_entry):
        """
        Do delocalization analysis. To do this, one considers:
            i) sampling region of planar averaged electrostatic potential (freysoldt approach)
            ii) sampling region of atomic site averaged potentials (kumagai approach)
            iii) structural relaxation amount outside of radius considered in kumagai approach (default is wigner seitz radius)
            iv) if defect is not a vacancy type -> track to see how much the defect has moved

        calculations that fail delocalization get "is_compatibile" set to False in parameters
        also parameters recieves a "delocalization_meta" with following dict:
            plnr_avg = {'is_compatible': True/False, 'metadata': metadata used for determining this}
            atomic_site = {'is_compatible': True/False, 'metadata': metadata used for determining this}
            structure_relax = {'is_compatible': True/False, 'metadata': metadata used for determining this}
            defectsite_relax = {'is_compatible': True/False, 'metadata': metadata used for determing this}
        """

        if 'freysoldt_meta' in defect_entry.parameters.keys():
            defect_entry = self.is_freysoldt_delocalized(defect_entry)
        else:
            print('Insufficient information provided for performing Freysoldt '
                  'correction delocalization analysis.\n'
                  'Cannot perform planar averaged electrostatic potential '
                  'compatibility analysis.')


        if 'kumagai_meta' in defect_entry.parameters.keys():
            defect_entry = self.is_kumagai_delocalized(defect_entry)
        else:
            print('Insufficient information provided for performing Kumagai '
                  'correction delocalization analysis.\n'
                  'Cannot perform atomic site averaged electrostatic '
                  'potential compatibility analysis.')


        if ('final_defect_structure' in defect_entry.parameters.keys()) and \
                ('final_structure' in defect_entry.parameters.keys()):
            defect_entry = self.is_final_relaxed_structure_delocalized(defect_entry)
        else:
            print('final_defect_structure does not exist in defect_entry.parameters. '
                  'Cannot perform full structure site relaxation compatibility analysis.')

        #TODO: overview check to see if is_compatible is set for delocalization...
        defect_entry.parameters.update({'is_compatible': is_compatible})

        return defect_entry

    def is_freysoldt_delocalized(self, defect_entry):
        plnr_avg_analyze_meta = {}
        plnr_avg_allows_compatible = True
        for ax in range(3):
            freystats = defect_entry.parameters['freysoldt_meta']['pot_corr_uncertainty_md'][ax]['stats']

            frey_variance_compatible = True if freystats['variance'] <= self.plnr_avg_var_tol else False
            frey_window = abs(freystats['minmax'][1] - freystats['minmax'][0])
            frey_minmax_compatible = True if frey_window <= self.plnr_avg_minmax_tol else False

            plnr_avg_analyze_meta[ax].update({'frey_variance_compatible': frey_variance_compatible,
                                              'frey_variance': freystats['variance'],
                                              'plnr_avg_var_tol': self.plnr_avg_var_tol,
                                              'frey_minmax_compatible': frey_minmax_compatible,
                                              'frey_minmax_window': frey_window,
                                              'plnr_avg_minmax_tol': self.plnr_avg_minmax_tol})

            if (not frey_variance_compatible) or (not frey_minmax_compatible):
                plnr_avg_allows_compatible = False

        defect_entry.parameters.update({'delocalization_meta': {'plnr_avg':
                                                                    {'is_compatible': plnr_avg_allows_compatible,
                                                                     'metadata': plnr_avg_analyze_meta}}})
        return defect_entry

    def is_kumagai_delocalized(self, defect_entry):
        atomic_site_analyze_meta = {}
        kumagaistats = defect_entry.parameters['kumagai_meta']['pot_corr_uncertainty_md']['stats']

        kumagai_variance_compatible = True if kumagaistats['variance'] <= self.atomic_site_var_tol else False
        kumagai_window = abs(kumagaistats['minmax'][1] - kumagaistats['minmax'][0])
        kumagai_minmax_compatible = True if kumagai_window <= self.atomic_site_minmax_tol else False

        atomic_site_analyze_meta.update({'kumagai_variance_compatible': kumagai_variance_compatible,
                                         'kumagai_variance': kumagaistats['variance'],
                                         'atomic_site_var_tol': self.atomic_site_var_tol,
                                         'kumagai_minmax_compatible': kumagai_minmax_compatible,
                                         'kumagai_minmax_window': kumagai_window,
                                         'plnr_avg_minmax_tol': self.atomic_site_minmax_tol})

        atomic_site_allows_compatible = True if (
            kumagai_variance_compatible and kumagai_minmax_compatible) else False

        defect_entry.parameters.update({'delocalization_meta': {'atomic_site':
                                                                    {'is_compatible': atomic_site_allows_compatible,
                                                                     'metadata': atomic_site_analyze_meta}}})
        return defect_entry

    def is_final_relaxed_structure_delocalized(self, defect_entry):
        structure_relax_analyze_meta = {}
        sc_scale = defect_entry.parameters[
            'scaling_matrix'] if 'scaling_matrix' in defect_entry.parameters.keys() else 1
        initial_defect_structure = defect_entry.defect.generate_defect_structure(sc_scale)
        final_defect_structure = defect_entry.parameters["final_structure"]
        radius_to_sample = defect_entry.parameters["kumagai_meta"]['sampling_radius']

        initsites, finalsites = [], []
        if type(defect_entry.defect) != Vacancy:
            sga = SpacegroupAnalyzer(initial_defect_structure)
            periodic_struc = sga.get_symmetrized_structure()
            poss_deflist = sorted(
                periodic_struc.get_sites_in_sphere(defect_entry.defect.site.coords, 2, include_index=True), key=lambda x: x[1])
            defindex = poss_deflist[0][2]
        else:
            defindex = None

        for site_ind in range(len(initial_defect_structure)):
            initsites.append(initial_defect_structure[site_ind].frac_coords)
            finalsites.append(final_defect_structure[site_ind].frac_coords)

        distmatrix = initial_defect_structure.lattice.get_all_distance(finalsites, initsites)

        distdata = []
        totpert = 0.
        for ind in range(len(finalsites)):
            if ind == defindex:
                continue
            else:
                totpert += distmatrix[ind, ind]
                # append [distance to defect, distance traveled, index in structure]
                distdata.append([distmatrix[ind, defindex], distmatrix[ind, ind], ind])

        distdata.sort()
        tot_relax_outside_wsrad = 0.
        perc_relax_outside_wsrad = 0.
        for newind in range(len(distdata)):
            distdata[newind].append(100 * distdata[newind][1] / totpert)  # append percentage for relaxation in
            if distdata[newind][0] > radius_to_sample:
                tot_relax_outside_wsrad += distdata[newind][1]
                perc_relax_outside_wsrad += distdata[newind][3]

        structure_tot_relax_compatible = True if tot_relax_outside_wsrad <= self.tot_relax_tol else False
        structure_perc_relax_compatible = True if perc_relax_outside_wsrad <= self.perc_relax_tol else False
        structure_relax_analyze_meta.update({'structure_tot_relax_compatible': structure_tot_relax_compatible,
                                             'tot_relax_outside_wsrad': tot_relax_outside_wsrad,
                                             'tot_relax_tol': self.tot_relax_tol,
                                             'structure_perc_relax_compatible': structure_perc_relax_compatible,
                                             'perc_relax_outside_wsrad': perc_relax_outside_wsrad,
                                             'perc_relax_tol': self.perc_relax_tol,
                                             'full_structure_relax_data': distdata,
                                             'defect_index': defindex})

        structure_relax_allows_compatible = True if (
            structure_tot_relax_compatible and structure_perc_relax_compatible) else False

        defect_entry.parameters.update({'delocalization_meta':{'structure_relax':
                                                                   {'is_compatible': structure_relax_allows_compatible,
                                                                    'metadata': structure_relax_analyze_meta}}})

        #NEXT: do single defect delocalization analysis (requires similar data, so might as well run in tandem
        # with structural delocalizaiton
        defectsite_relax_analyze_meta = {}
        if type(defect_entry.defect) == Vacancy:
            defectsite_relax_allows_compatible = True
            defectsite_relax_analyze_meta.update({'relax_amount': None,
                                                  'defect_tot_relax_tol': self.defect_tot_relax_tol})
        else:
            defect_relax_amount = distmatrix[defindex, defindex]
            defectsite_relax_allows_compatible = True if defect_relax_amount <= self.defect_tot_relax_tol else False
            defectsite_relax_analyze_meta.update({'relax_amount': defect_relax_amount,
                                                  'defect_tot_relax_tol': self.defect_tot_relax_tol})

        defect_entry.parameters.update({'delocalization_meta': {'defectsite_relax':
                                                                    {'is_compatible': defectsite_relax_allows_compatible,
                                                                     'metadata': defectsite_relax_analyze_meta}}})
        return defect_entry