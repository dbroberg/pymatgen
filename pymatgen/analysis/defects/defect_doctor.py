"""
A set of classes and functions that are useful for defect workflow management
"""

from pymatgen.io.vasp.sets import MPRelaxSet
# from pymatgen.analysis.defects.defect_compatibility import DefectCompatibility
# from pymatgen.analysis.defects.defect_builder import TaskDefectBuilder
# from pymatgen.analysis.defects.thermodynamics import DefectPhaseDiagram
# from pymatgen.analysis.structure_matcher import StructureMatcher

from atomate.vasp.fireworks.core import TransmuterFW

from fireworks import Workflow



def get_fw_from_defect( defect, supercell_size,
                        defect_input_set=None,
                        job_type='normal', db_file='>>db_file<<', vasp_cmd='>>vasp_cmd<<'):

    chgdef_trans = ["DefectTransformation"]
    chgdef_trans_params = [{"scaling_matrix": supercell_size,
                            "defect": defect.copy()}]

    test_bulk_struct = defect.bulk_structure.copy()
    test_bulk_struct.make_supercell( supercell_size)
    num_atoms = len(test_bulk_struct)

    def_tag = "{}:{}_{}_{}_{}atoms".format(test_bulk_struct.composition.reduced_formula,
                                           job_type, defect.name, defect.charge, num_atoms)


    if (job_type == 'normal') and (defect_input_set is None):
        defect_sc = defect.generate_defect_structure(supercell=supercell_size)

        reciprocal_density = 50 if job_type == 'hse' else 100
        kpoints_settings = {"reciprocal_density": reciprocal_density}
        stdrd_defect_incar_settings = {"EDIFF": 0.0001, "EDIFFG": 0.001, "IBRION": 2, "ISMEAR": 0, "SIGMA": 0.05,
                                       "ISPIN": 2, "ISYM": 2, "LVHAR": True, "LVTOT": True, "NSW": 100,
                                       "NELM": 60, "ISIF": 2, "LAECHG": False, "LWAVE": True}
        defect_input_set = MPRelaxSet(defect_sc,
                                      user_incar_settings=stdrd_defect_incar_settings.copy(),
                                      user_kpoints_settings=kpoints_settings,
                                      use_structure_charge=True)

    elif defect_input_set is None:
        raise ValueError("job_type = {} and no defect input set is specified...need to specify input set".format( job_type))


    fw = TransmuterFW(name=def_tag, structure=defect.bulk_structure.copy(),
                      transformations=chgdef_trans,
                      transformation_params=chgdef_trans_params,
                      vasp_input_set=defect_input_set,
                      vasp_cmd=vasp_cmd,
                      copy_vasp_outputs=False,
                      db_file=db_file,
                      job_type=job_type,
                      bandstructure_mode="auto")

    return fw




class DefectResubber(object):
    """
    Takes a defect builder result and creates follow up
    workflow outputs as needed (from atomate)
    for completing the defect thermodynamics desired.

    Does this based on supercell size scaling, delocalization metrics...
    """

    def __init__(self, lpad):
        self.lpad = lpad

    def additional_charges(self, defect_phase_diagram, supercell_size, name_wf='chg_defect_wf', submit = False):
        """
        Submit additional charges for a given defect_phase_diagram

        :param defect_phase_diagram:
        :param name_wf:
        :return:
        """
        print("Considering additional charges to run for defect_phase_diagram")
        fws = []

        print('Finished charges:')
        for k, v in defect_phase_diagram.finished_charges.items():
            print('\t{}: {}'.format( k, v))
        print('\nSTABLE charges:')
        for k, v in defect_phase_diagram.stable_charges.items():
            print('\t{}: {}\t(t.l.: {} ) '.format( k, v, defect_phase_diagram.transition_levels[k]))

        print('\nNow consider charges for follow up..')
        rec_dict = defect_phase_diagram.suggest_charges()
        for defname, charge_list in rec_dict.items():
            defect_template = defect_phase_diagram.stable_entries[defname][0].defect.copy()
            for charge in charge_list:
                defect = defect_template.copy()
                defect.set_charge( charge)
                fws.append( get_fw_from_defect( defect, supercell_size) )
                print('\trerunning ', defect.name, defect.charge)

        print("Created a total of {} fireworks".format( len(fws)))
        # load Workflow to lpad
        if submit:
            wf = Workflow( fws, name=name_wf)
            self.lpad.add_wf( wf)
            print('Submitted!')

        return

    # def larger_supercells(self, entries, max_atoms=800):
    #     """
    #     Use compatibility class in a phaes diagram without filtering to determine whether
    #     larger supercells should be run
    #     :param entries:
    #     :return:
    #     """
    #     vbm, band_gap,
    #     DefectPhaseDiagram( entries, vbm, band_gap, filter_compatible=False)
    #     #TODO -> based on compatibility, allow for larger supercells to be run...
    #     return



if __name__ == "__main__":
    pass
    #get database
    # from atomate.vasp.database import VaspCalcDb
    # db = VaspCalcDb.from_db_file("db.json", admin=True)
    # tasks = db.get_collections('tasks')
    # defects?? #TODO make a defect collection for storing defect entries...
    #
    # drs = DefectResubber(tasks, defects, lpad)
    # drs.resubmit( bulk_struct, name_wf)



