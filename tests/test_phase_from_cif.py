import src.cctbx_dep.phase_from_cif as pfc

test_file_path = "./data/cifs/1000032.cif"
uc = (4.7605, 4.7605, 12.9956, 90, 90, 120)

def test_load_cif():
   #load_cif() is called in RietveldPhases' __init__
   phase_structure_dict = pfc.load_cif(test_file_path)
   assert 'structure' in phase_structure_dict
   labels = ['Al1','O1']
   element_symbols = ['Al','O']
   occupancy = 1.0
   for x,l,s in zip(phase_structure_dict["structure"].scatterers(),
                    labels,element_symbols):
      assert x.label == l
      assert np.isclose(x.occupancy,occupancy)
      assert x.element_symbol() == s
      print x.report_details(phase_structure_dict["unit_cell"],'')

   assert 'unit_cell' in phase_structure_dict
   print phase_structure_dict["unit_cell"].parameters()
   for x,y in zip(phase_structure_dict["unit_cell"].parameters(), uc):
      assert np.isclose(x,y)

   assert 'crystal_system' in phase_structure_dict
   assert type(phase_structure_dict["crystal_system"]) is str
   assert phase_structure_dict["crystal_system"] == 'Trigonal'
   # test_phase.assemble_lattice_params()
   # assert test_phase.crystal_system == 'H'

   assert 'chemical_name' in phase_structure_dict
   assert phase_structure_dict["chemical_name"] == 'Corundum'