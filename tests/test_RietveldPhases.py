from __future__ import division
import numpy as np
import time
import matplotlib.pyplot as plt
from random import randrange
import profile, pstats
import pytest
import copy

import sys, os
sys.path.append(os.path.abspath(".."))

import src.RietveldPhases as RietveldPhases
from src.RietveldPhases import RietveldPhases as Rp
from src.RietveldRefinery import RietveldRefinery

from cctbx.eltbx import wavelengths

# class RietveldPhasesTest(unittest.TestCase):

#    def setUp(self):
Rp.set_profile(
   r"./data/profiles/Jade-Al2O3-Sim.xye",number_of_columns=2)

def test_data_read_in():
   assert len(Rp.two_theta) == 4250
   assert np.isclose(Rp.d_min,1.08934)
   assert np.isclose(Rp.d_max,17.58881)
   assert len(Rp.I) == 4250
   assert np.isclose(Rp.two_theta[-1],90)
   assert len(Rp.sigma) == 4250
   assert np.isclose(Rp.sigma[-1],np.sqrt(Rp.I[-1]))

def test_global_parameters_exist():
   Rp_dict = Rp.__dict__
   assert 'two_theta_0' in Rp_dict
   assert 'bkgd' in Rp_dict

test_phase = Rp("./data/cifs/1000032.cif")
uc = (4.7605, 4.7605, 12.9956, 90, 90, 120)

def test_parameters_exist():
   tp_dict = test_phase.__dict__
   print tp_dict
   assert 'U' in tp_dict
   assert 'V' in tp_dict
   assert 'W' in tp_dict
   assert 'Scale' in tp_dict
   assert 'eta' in tp_dict
   assert 'lattice_parameters' in tp_dict

def test_set_bkgd_order():
   Rp.set_bkgd_order(3)
   assert len(Rp.bkgd) == 3
   assert len(test_phase.bkgd) == 3

def test_bkgd_param_gen():
   gen = Rp.bkgd_param_gen(order=3)
   for n in xrange(0,3):
      assert next(gen)['labels'] == 'bkgd_' + str(n)
   assert pytest.raises(StopIteration,next,gen)

def test_background_polynomial():
   Rp.bkgd['values'][0] = 1
   assert np.all(np.isclose(Rp.background_polynomial(),
      np.ones(len(Rp.two_theta),dtype=float)))

   Rp.bkgd['values'] = np.array([1,2,3],dtype=float)
   assert np.all(np.isclose(Rp.background_polynomial(),
      1.0+2.0*Rp.two_theta+3.0*Rp.two_theta_powers[2,:]))

   Rp.bkgd['values'] = np.array([0,0,0],dtype=float)
   assert np.all(np.isclose(Rp.background_polynomial(),
      np.zeros(len(Rp.two_theta),dtype=float)))

def test_U_default():
   assert test_phase.U.dtype == RietveldPhases.custom_dtype
   assert test_phase.U['labels'] == RietveldPhases.default_U['labels']
   assert np.isclose(test_phase.U['values'], RietveldPhases.default_U['values'])
   assert np.isclose(test_phase.U['l_limits'],
      RietveldPhases.default_U['l_limits'])
   assert np.isclose(test_phase.U['u_limits'],
      RietveldPhases.default_U['u_limits'])

def test_set_vertical_offset():
   assert test_phase.vertical_offset == False
   assert 'cos_theta' not in test_phase.__dict__
   test_phase.set_vertical_offset(True)
   assert test_phase.vertical_offset == True
   assert 'cos_theta' in Rp.__dict__
   assert np.isclose(Rp.cos_theta[-1],-360/np.pi/np.sqrt(2))

def test_LP_intensity_scaling():
   assert len(Rp.two_theta) == len(Rp.LP_intensity_scaling())

   Rp.two_theta_0['values'] = 0.1
   two_theta = Rp.two_theta - Rp.two_theta_0['values']
   assert np.all(np.isclose(Rp.LP_intensity_scaling(),
      (1+np.cos(np.pi/180*two_theta)**2) \
         /np.sin(np.pi/360*two_theta) \
         /np.sin(np.pi/180*two_theta)))
   Rp.two_theta_0['values'] = 0.0

def test_assemble_global_x():
   #assemble_global_x() is called in RietveldPhases' __init__
   assert 'global_x' in Rp.__dict__
   assert 'global_x_no_bkgd_mask' in Rp.__dict__
   assert len(Rp.global_x) == len(Rp.bkgd) + 1 # + 1 for two_theta_0

def test_update_global_x():
   mask = np.zeros(len(Rp.global_x),dtype=bool)
   mask[0] = True
   new_x = copy.deepcopy(Rp.global_x)
   new_x['values'][0] = 0.3
   Rp.update_global_x(new_x,mask)
   assert np.isclose(Rp.global_x['values'][0],0.3)

def test_global_param_gen():
   gen = Rp.global_param_gen()
   assert next(gen)['labels'] == 'two_theta_0'
   assert np.all(np.char.startswith(next(gen)['labels'],'bkgd_'))
   assert pytest.raises(StopIteration,next,gen)

def test_phase_param_gen():
   count = 0
   for x in test_phase.phase_param_gen():
      count += 1
   assert count == 6

def test_load_cif():
   #load_cif() is called in RietveldPhases' __init__
   tp_dict = test_phase.__dict__
   assert 'structure' in tp_dict
   labels = ['Al1','O1']
   element_symbols = ['Al','O']
   occupancy = 1.0
   for x,l,s in zip(test_phase.structure.scatterers(),labels,element_symbols):
      assert x.label == l
      assert np.isclose(x.occupancy,occupancy)
      assert x.element_symbol() == s
      print x.report_details(test_phase.unit_cell,'')

   assert 'unit_cell' in tp_dict
   print test_phase.unit_cell.parameters()
   for x,y in zip(test_phase.unit_cell.parameters(), uc):
      assert np.isclose(x,y)

   assert 'crystal_system' in tp_dict
   assert type(test_phase.crystal_system) is str
   test_phase.load_cif(test_phase.fn_cif,d_min = Rp.d_min)
   assert test_phase.crystal_system == 'Trigonal'
   test_phase.assemble_lattice_params()
   assert test_phase.crystal_system == 'HTrigonal'

   assert 'chemical_name' in tp_dict
   assert test_phase.chemical_name == 'Corundum'

def test_assemble_lattice_params():
   assert 'uc_mask' in test_phase.__dict__
   assert len(test_phase.uc_mask) == 6

def test_unit_cell_parameter_gen():
   tp_uc_gen = test_phase.unit_cell_parameter_gen()
   assert next(tp_uc_gen)[0] == 'uc_a'
   assert np.isclose(next(tp_uc_gen)[1],test_phase.unit_cell.parameters()[2])

   uc_mask = copy.deepcopy(test_phase.uc_mask)
   test_phase.uc_mask = np.ones(6,dtype=bool)
   uc_params = test_phase.unit_cell.parameters()
   for i,x in enumerate(test_phase.unit_cell_parameter_gen()):
      assert np.isclose(uc_params[i],x[1])
   test_phase.uc_mask = uc_mask

def test_assemble_phase_x():
   #assemble_phase_x() is called in RietveldPhases' __init__
   tp_dict = test_phase.__dict__
   assert 'phase_x' in tp_dict
   assert 'global_and_phase_x' in tp_dict
   assert 'global_mask_no_bkgd' in tp_dict
   assert 'phase_mask' in tp_dict

   x = len(test_phase.global_and_phase_x)
   assert x == len(test_phase.phase_x) + 1
   assert x == len(test_phase.global_mask_no_bkgd)
   assert x == len(test_phase.phase_mask)

   assert str(type(test_phase.U))[7:17] == 'numpy.void'
   assert len(test_phase.U) == 4

   assert len(test_phase.eta) == test_phase.eta_order
   assert len(test_phase.lattice_parameters) == 2

def test_update_params():
   tmp_x = np.copy(test_phase.phase_x)
   mask = np.char.startswith(tmp_x['labels'],'Sca')
   test_val = 100
   tmp_x['values'][mask] = 100
   test_phase.update_params(tmp_x,mask=mask)
   assert np.isclose(test_phase.Scale['values'],test_val)
   mask = np.char.startswith(tmp_x['labels'],'uc_a')
   print tmp_x['values'][mask]
   tmp_x['values'][mask] = 2
   print tmp_x['values'][mask]
   test_phase.update_params(tmp_x)
   print test_phase.recompute_peak_positions
   print test_phase.phase_x
   print test_phase.lattice_parameters
   print test_phase.unit_cell.parameters()
   #assert False TODO: fix lattice_parameter part of test

def test_eta_polynomial():
   assert 'eta' in test_phase.__dict__
   test_phase.set_eta_order(2)
   test_phase.assemble_phase_x()
   test_eta = np.array([0.5,0.005])
   test_phase.eta['values'] = test_eta
   tmp = test_phase.eta_polynomial()
   assert np.isclose(tmp[-1],test_eta[0]+Rp.two_theta[-1]*test_eta[1])
   test_phase.eta['values'] = np.array([0.5,0])

def test_compute_relative_intensities():
   tp_dict = test_phase.__dict__
   assert 'f_miller_set' in tp_dict
   assert 'crystal_density' in tp_dict
   assert 'd_spacings' in tp_dict
   assert 'relative_intensities' in tp_dict

   assert np.all(np.isclose(test_phase.d_spacings, np.array(
      [ 3.48114434,  2.55177292,  2.38025   ,  2.0860757 ,  1.74057217,
        1.60196736,  1.54715716,  1.51527741,  1.51135808,  1.40499663,
        1.37423798,  1.3364588 ,  1.27588646,  1.23944057,  1.23455034,
        1.19354167,  1.190125  ,  1.14760201,  1.12613194,  1.12452033,
        1.09932894])))
   assert np.all(np.isclose(test_phase.relative_intensities, np.array(
      [ 0.19906132,  0.61069362,  0.29855677,  1.18211397,  0.86296333,
        2.11656703,  0.05664854,  0.08304719,  0.18443051,  1.14926644,
        1.85168419,  0.03841942,  0.0587068 ,  0.69865033,  0.35418387,
        0.03049309,  0.29644002,  0.24097501,  0.18603062,  0.14235153,
        0.40296412])))

def exercise_RietveldPhases():
   cifs = ["1000032.cif","1507774.cif"]

   # Testing LP_intensity_scaling
   assert np.isclose(Rt[0].LP_Intensity_Scaling(20.0),31.7054214503)

   # Testing two_theta outputs
   assert np.all(np.isclose(Rt[0].two_theta_peaks[0:10], np.array(
      [[ 25.56751716],
       [ 35.13879293],
       [ 37.76313851],
       [ 43.33853666],
       [ 46.16161879],
       [ 52.53263976],
       [ 57.47974186],
       [ 59.71809269],
       [ 61.10698338],
       [ 61.28247987]])))

   #Testing Weighted Intensity outputs
   assert np.all(np.isclose(Rt[0].weighted_intensities[0:10], np.array(
      [[ 0.19906132],
       [ 0.61069362],
       [ 0.29855677],
       [ 1.18211397],
       [ 0.01740687],
       [ 0.86296333],
       [ 2.11656703],
       [ 0.05664854],
       [ 0.08304719],
       [ 0.18443051]])))

   #Testing PseudoVoigtProfile
   if display_plots:
      for RV in Rt:
         #Select a random peak:
         if len(RV.two_theta_peaks) != 0:
            rnd_index = randrange(0,len(RV.two_theta_peaks),1)
            RV.showPVProfilePlot("Test",rnd_index, autohide=False)

   # #Testing Refinery
   # RR = RietveldRefinery(Rt,minimizer_input_string,
   #    store_intermediate_state=True) #,show_plots=False)

   # # RR.display(RR.minimize_Scale_Offset)
   # RR.display(RR.minimize_Scale_Offset_W)
   # RR.display(RR.minimize_All)



def run():
   # exercise_RietveldPhases()
   print "OK"

if (__name__ == "__main__"):
   # pr = profile.Profile()
   # pr.enable()
   # run()
   unittest.main(warnings='ignore')
   # profile.run('run(); print')
   # pr.disable()
   # s = StringIO.StringIO()
   # sortby = 'cumulative'
   # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
   # ps.print_stats()
   # print s.getvalue()