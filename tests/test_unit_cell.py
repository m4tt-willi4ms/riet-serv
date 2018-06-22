import numpy as np
import copy
import pytest

import src.cctbx_dep.unit_cell as uc
from src.rietveld_phases import RietveldPhases

@pytest.fixture(scope="module")
def test_phase():
   return RietveldPhases("./data/cifs/1000032.cif")

def test_assemble_lattice_parameters(test_phase):
   test_phase_settings = test_phase.phase_settings
   assert 'uc_mask' in test_phase_settings
   assert len(test_phase_settings["uc_mask"]) == 6

   assert test_phase_settings["crystal_system_trigonal"] == 'H'
   assert 'chemical_name' in test_phase_settings
   assert test_phase_settings["chemical_name"] == 'Corundum'

def test_unit_cell_parameter_gen(test_phase):
   tp_uc_gen = uc.unit_cell_parameter_gen(test_phase.phase_settings)
   uc_params = test_phase.phase_settings["unit_cell"].parameters()
   assert next(tp_uc_gen)[0] == 'uc_a'
   assert np.isclose(next(tp_uc_gen)[1], uc_params[2])
   uc_mask = copy.deepcopy(test_phase.phase_settings["uc_mask"])
   test_phase.phase_settings["uc_mask"] = np.ones(6,dtype=bool)
   for i, x in enumerate(
      uc.unit_cell_parameter_gen(test_phase.phase_settings)):
      assert np.isclose(uc_params[i], x[1])

def test_update_unit_cell(test_phase):
   settings = test_phase.phase_settings
   lattice_parameters = 12*np.random.random(size=6)
   uc.update_unit_cell(settings, lattice_parameters)
   new_uc = settings["unit_cell"].parameters()
   uc_mask = settings["uc_mask"]
   from itertools import compress
   for x, y  in zip(compress(list(new_uc), uc_mask), list(lattice_parameters)):
      assert x, y

   settings['crystal_system'] = 'Triclinic'
   lattice_parameters = 20*np.random.random(size=6)+80
   uc.update_unit_cell(settings, lattice_parameters)
   new_uc = settings["unit_cell"].parameters()
   for x, y in zip(list(new_uc), list(lattice_parameters)):
      assert x == y

@pytest.fixture(scope="module")
def rutile_phase():
   return RietveldPhases("./data/cifs/9015662-rutile.cif")

def test_inverse_filter(rutile_phase):
   uc_params = rutile_phase.phase_settings["unit_cell"].parameters()
   assert rutile_phase.phase_settings["crystal_system"] == 'Tetragonal'
   assert rutile_phase.phase_settings["uc_mask"] == \
         [True, False, True, False, False, False]
   assert uc_params == (4.5937, 4.5937, 2.9587, 90, 90, 90)
