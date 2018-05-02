import numpy as np
import copy

import src.cctbx_dep.unit_cell as uc
from src.rietveld_phases import RietveldPhases

test_phase = RietveldPhases("./data/cifs/1000032.cif")

def test_assemble_lattice_params():
   assert 'uc_mask' in test_phase.phase_settings
   assert len(test_phase.phase_settings["uc_mask"]) == 6

def test_unit_cell_parameter_gen():
   tp_uc_gen = uc.unit_cell_parameter_gen(test_phase.phase_settings)
   uc_params = test_phase.phase_settings["unit_cell"].parameters()
   assert next(tp_uc_gen)[0] == 'uc_a'
   assert np.isclose(next(tp_uc_gen)[1],uc_params[2])
   uc_mask = copy.deepcopy(test_phase.phase_settings["uc_mask"])
   test_phase.phase_settings["uc_mask"] = np.ones(6,dtype=bool)
   for i, x in enumerate(
      uc.unit_cell_parameter_gen(test_phase.phase_settings)):
      assert np.isclose(uc_params[i],x[1])
