import numpy as np

import src.cctbx_dep.phase_from_cif as pfc
import src.cctbx_dep.unit_cell as unit_cell
from src.rietveld_phases import RietveldPhases as Rp

test_file_path = "./data/cifs/1000032.cif"
uc = (4.7605, 4.7605, 12.9956, 90, 90, 120)

def test_load_cif():
   #load_cif() is called in RietveldPhases' __init__
   test_phase = Rp(test_file_path)
   test_phase_settings = test_phase.phase_settings
   assert 'structure' in test_phase_settings
   labels = ['Al1','O1']
   element_symbols = ['Al','O']
   occupancy = 1.0
   for x,l,s in zip(test_phase_settings["structure"].scatterers(),
                    labels, element_symbols):
      assert x.label == l
      assert np.isclose(x.occupancy,occupancy)
      assert x.element_symbol() == s
      print x.report_details(test_phase_settings["unit_cell"],'')

   assert 'unit_cell' in test_phase_settings
   print test_phase_settings["unit_cell"].parameters()
   for x,y in zip(test_phase_settings["unit_cell"].parameters(), uc):
      assert np.isclose(x,y)

   assert 'crystal_system' in test_phase_settings
   assert type(test_phase_settings["crystal_system"]) is str
   d_min = test_phase_settings["d_min"]
   pfc.load_cif(test_phase.file_path_cif, d_min = d_min)
   assert test_phase_settings["crystal_system"] == 'Trigonal'
   unit_cell.assemble_lattice_params(test_phase.phase_settings, Rp.custom_dtype)
   assert test_phase_settings["crystal_system_trigonal"] == 'H'

   assert 'chemical_name' in test_phase_settings
   assert test_phase_settings["chemical_name"] == 'Corundum'