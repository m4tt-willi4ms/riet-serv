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

def test_parameters_exist():
   tp_dict = test_phase.__dict__
   print tp_dict
   assert 'U' in tp_dict
   assert 'V' in tp_dict
   assert 'W' in tp_dict
   assert 'Amplitude' in tp_dict
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

Rp.assemble_global_x()

def test_assemble_global_x():
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
   order = 4 
   Rp.set_bkgd_order(order)
   gen = Rp.global_param_gen()
   assert next(gen)['labels'] == 'two_theta_0'
   assert np.all(np.char.startswith(next(gen)['labels'],'bkgd_'))
   assert pytest.raises(StopIteration,next,gen)


def exercise_RietveldPhases():
   # RietveldPhase.fromstring(input_string)
   cifs = ["1000032.cif","1507774.cif"]
   Rt = []

   #Testing Read-in from input_string
   assert np.isclose(Rt[0].x['values'][Rt[0].U_index], 0.0)
   assert Rt[0].x['labels'][Rt[0].U_index] == 'U'
   assert np.isclose(Rt[1].x['values'][Rt[1].U_index], 0.00)
   assert Rt[1].x['labels'][Rt[1].U_index] == 'U'
   assert np.isclose(Rt[0].x['values'][Rt[0].V_index], 0.0)
   assert Rt[0].x['labels'][Rt[0].V_index] == 'V'
   assert np.isclose(Rt[0].x['values'][Rt[0].W_index],  0.01)
   # assert np.isclose(Rt[0].x['values'][Rt[0].Amplitude_index],  1)
   assert np.isclose(RietveldPhases.x['values'] \
      [RietveldPhases.two_theta_0_index] , np.array([0.001]))
   assert len(RietveldPhases.x['values'] \
      [np.char.startswith(RietveldPhases.x['labels'],"Bkgd")]) == 3
   assert np.all(np.isclose(RietveldPhases.x['values'] \
      [np.char.startswith(RietveldPhases.x['labels'],"Bkgd")], \
      np.array([0.,0.,0.])))

   # Testing iotbx output from .cif card
   np.set_printoptions(threshold=None)
   assert str(Rt[0].unit_cell.parameters()) == \
   """(4.7605, 4.7605, 12.995599999999998, 90.0, 90.0, 120.00000000000001)"""
   assert str(Rt[1].unit_cell.parameters()) == \
   """(4.1769, 4.7218, 58.3385, 89.44, 89.634, 75.854)"""

   #Testing d-spacing outputs
   assert np.all(np.isclose(Rt[0].d_spacings[0:30], np.array(
      [  3.48114434e+00,   2.55177292e+00,   2.38025000e+00,
          2.08607570e+00,   1.96485460e+00,   1.74057217e+00,
          1.60196736e+00,   1.54715716e+00,   1.51527741e+00,
          1.51135808e+00,   1.40499663e+00,   1.37423798e+00,
          1.33645880e+00,   1.27588646e+00,   1.23944057e+00,
          1.23455034e+00,   1.19354167e+00,   1.19012500e+00,
          1.16038145e+00,   1.16038145e+00,   1.14760201e+00,
          1.13903464e+00,   1.12613194e+00,   1.12452033e+00,
          1.09932894e+00])))
   assert np.all(np.isclose(Rt[1].d_spacings[0:30], np.array(
      [ 14.58380434,  11.66704347,   9.72253623,   8.33360248,
         7.29190217,   6.48169082,   5.83352174,   5.30320158,
         4.86126811,   4.57845246,   4.56140325,   4.5289481 ,
         4.51721072,   4.48732441,   4.4650281 ,   4.35771092,
         4.27435487,   4.24980397,   4.16680124,   4.1557422 ,
         4.02715497,   4.00946117,   3.99846875,   3.96837552,
         3.90666115,   3.89836103,   3.89236416,   3.88901449,
         3.83110776,   3.75460982])))

   #Testing Relative Intensity Outputs
   assert np.all(np.isclose(Rt[0].relative_intensities[0:30], np.array(
      [ 0.19906132,  0.61069362,  0.29855677,  1.18211397,  0.01740687,
        0.86296333,  2.11656703,  0.05664854,  0.08304719,  0.18443051,
        1.14926644,  1.85168419,  0.03841942,  0.0587068 ,  0.69865033,
        0.35418387,  0.03049309,  0.29644002,  0.02099137,  0.02108226,
        0.24097501,  0.01391577,  0.18603062,  0.14235153,  0.40296412])))
   # print repr(Rt[1].relative_intensities[0:30])
   assert np.all(np.isclose(Rt[1].relative_intensities[0:30], np.array(
      [ 0.00991668,  0.01230013,  0.00716603,  0.00855396,  0.00426061,
        0.00530974,  0.00224645,  0.00295588,  0.00115554,  0.02576297,
        0.00280861,  0.00032264,  0.02274331,  0.00189005,  0.00051361,
        0.00231178,  0.00158587,  0.00066494,  0.00048172,  0.00081688,
        0.00166894,  0.00044647,  0.00100012,  0.00176042,  0.00089117,
        0.00087474,  0.00128861,  0.0011674 ,  0.0561287 ,  0.00185326])))

   #Testing test-data read-in
   if is_Sim_data:
      np.set_printoptions(threshold=None)
      assert str(tst_two_theta) == \
      """[  5.     5.02   5.04 ...,  89.96  89.98  90.  ]"""
      assert str(tst_y) == \
      """[ 0.01  0.01  0.01 ...,  0.02  0.02  0.02]"""

   # Testing Background_Polynomial
   RietveldPhases.x['values'] \
      [np.char.startswith(RietveldPhases.x['labels'],"Bkgd")]  \
      = np.array([1.0,2.0,3.0])
   tst_bkgd_two_theta = np.array([0,1,2.3,100.5])
   assert np.all(np.isclose(RietveldPhases.Background_Polynomial( \
      tst_bkgd_two_theta) , np.array([  1,6,21.47,   3.05027500e+04])))
   #Revert back to flat background
   RietveldPhases.x['values'] \
      [np.char.startswith(RietveldPhases.x['labels'],"Bkgd")]  \
      = np.array([0.0,0.0,0.0])

   # Testing Eta_Polynomial
   Rt0_eta_mask = np.isin(np.array(range(0,len(RietveldPhases.x))), \
      np.array(range(Rt[0].eta_0_index,Rt[0].eta_0_index+Rt[0].eta_rank)))
   RietveldPhases.x['values'][Rt0_eta_mask] = np.array([0.5,0.005])
   tst_eta_two_theta = np.array([0,1,2.3,100.5])
   assert np.all(np.isclose(Rt[0].eta_Polynomial( \
      tst_eta_two_theta) , np.array([ 0.5, 0.505, 0.5115, 1.00249999])))
   #Revert back to unbiased eta
   RietveldPhases.x['values'][Rt0_eta_mask] = np.array([0.5,0.00])

   # Testing LP_Intensity_Scaling
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

   # # RR.display(RR.minimize_Amplitude_Offset)
   # RR.display(RR.minimize_Amplitude_Offset_W)
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