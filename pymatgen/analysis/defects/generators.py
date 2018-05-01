# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
"""
This module defines classes to generate point defect structures
"""

__author__ = "Danny Broberg, Shyam Dwaraknath, Bharat Medasani, Nils E. R. Zimmermann, Geoffroy Hautier"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Mar 15, 2018"

import six
import logging
from abc import ABCMeta, abstractmethod

from monty.json import MSONable

from pymatgen.core import PeriodicSite
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.defects.core import Vacancy, Interstitial, Substitution
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.defects.point_defects import StructureMotifInterstitial, TopographyAnalyzer

logger = logging.getLogger(__name__)


class DefectGenerator(six.with_metaclass(ABCMeta, MSONable)):
    """
    Abstract class for point defects
    Implements generator pattern
    """

    def __iter__(self):
        """
        Return self as this should be an iterator
        """
        return self

    @abstractmethod
    def __next__(self):
        """
        Abstract method to return defects
        """
        return


class VacancyGenerator(DefectGenerator):
    """
    Simple generator for vacancies based on periodically
    equivalent sites
    """

    def __init__(self, structure, include_bv_charge=False):
        """
        Initializes a Vacancy Generator
        Args:
            structure(Structure): pymatgen structure object
        """
        self.structure = structure
        self.include_bv_charge = include_bv_charge

        # Find equivalent site list
        sga = SpacegroupAnalyzer(self.structure)
        self.symm_structure = sga.get_symmetrized_structure()
        self.equiv_site_seq = list(self.symm_structure.equivalent_sites)

        self.struct_valences = None
        if self.include_bv_charge:
            bv = BVAnalyzer()
            self.struct_valences = bv.get_valences(self.structure)

    def __next__(self):
        """
        Returns the next vacancy in the sequence or
        raises StopIteration
        """
        if len(self.equiv_site_seq) > 0:
            vac_site = self.equiv_site_seq.pop(0)
            charge = 0.0
            if self.struct_valences:
                site_index = self.structure.get_sites_in_sphere(vac_site[0].coords, 0.1, include_index=True)[0][2]
                charge = -1 * self.struct_valences[site_index]

            return Vacancy(self.structure, vac_site[0], charge=charge)
        else:
            raise StopIteration


class InterstitialGenerator(DefectGenerator):
    """
    Generator for interstitials at positions
    where the interstitialcy is coordinated by nearest neighbors
    in a way that resembles basic structure motifs
    (e.g., tetrahedra, octahedra).  The algorithm is called InFiT
    (Interstitialcy Finding Tool), it was introducted by
    Nils E. R. Zimmermann, Matthew K. Horton, Anubhav Jain,
    and Maciej Haranczyk (Front. Mater., 4, 34, 2017),
    and it is used by the Python Charged Defect Toolkit
    (PyCDT: D. Broberg et al., Comput. Phys. Commun., in press, 2018).
    """

    def __init__(self, structure, element):
        """
        Initializes an Interstitial generator using structure motifs
        Args:
            structure (Structure): pymatgen structure object
            element (str or Element or Specie): element for the interstitial
        """
        self.structure = structure
        self.element = element
        interstitial_finder = StructureMotifInterstitial(self.structure, self.element)
        self.defect_sites = list(interstitial_finder.enumerate_defectsites())

        #for multiplicity, neccessary to get prim_structure
        spa = SpacegroupAnalyzer(self.structure, symprec=1e-2)
        prim_struct = spa.get_primitive_standard_structure()
        conv_prim_rat = int(self.structure.num_sites/prim_struct.num_sites)
        self.multiplicities = [int(interstitial_finder.get_defectsite_multiplicity(def_ind) / conv_prim_rat)
                               for def_ind in range(len(self.defect_sites))]

    def __next__(self):
        """
        Returns the next interstitial or
        raises StopIteration
        """
        if len(self.defect_sites) > 0:
            int_site = self.defect_sites.pop(0)
            mult = self.multiplicities.pop(0)

            return Interstitial(self.structure, int_site, multiplicity=mult)
        else:
            raise StopIteration


class VoronoiInterstitialGenerator(DefectGenerator):
    """
    Generator for interstitials based on a simple Voronoi analysis
    """

    def __init__(self, structure, element):
        """
        Initializes an Interstitial generator using Voronoi sites
        Args:
            structure (Structure): pymatgen structure object
            element (str or Element or Specie): element for the interstitial
        """
        self.structure = structure
        self.element = element

        framework = list(self.structure.symbol_set)
        get_voronoi = TopographyAnalyzer(self.structure, framework, [], check_volume=False)
        get_voronoi.cluster_nodes()
        get_voronoi.remove_collisions()

        #trim equivalent nodes with symmetry analysis
        struct_to_trim = self.structure.copy()
        for poss_inter in get_voronoi.vnodes:
            struct_to_trim.append( self.element, poss_inter.frac_coords, coords_are_cartesian=False)

        symmetry_finder = SpacegroupAnalyzer(struct_to_trim, symprec=1e-1)
        equiv_sites_list = symmetry_finder.get_symmetrized_structure().equivalent_sites

        self.equiv_site_seq = []
        for poss_site_list in equiv_sites_list:
            if poss_site_list[0] not in self.structure:
                self.equiv_site_seq.append(poss_site_list)


    def __next__(self):
        """
        Returns the next interstitial or
        raises StopIteration
        """
        if len(self.equiv_site_seq) > 0:
            inter_site_list = self.equiv_site_seq.pop(0)

            return Interstitial(self.structure, inter_site_list[0], multiplicity=len(inter_site_list))
        else:
            raise StopIteration