from __future__ import division
from scitbx import lbfgsb
# from scitbx.array_family import flex
# import scitbx.math
# from libtbx.test_utils import approx_equal, eps_eq, Exception_expected
# import sys
import numpy as np
import time
import matplotlib.pyplot as plt
from cctbx.eltbx import wavelengths
from random import randrange
from libtbx import test_utils
import libtbx.load_env

import sys, os
sys.path.append(os.path.abspath(".."))

from RietveldPhases import RietveldPhases
from RietveldRefinery import RietveldRefinery

input_strings = ["""\
U              0.0   -0.1   0.1
V              0.0   -0.1   0.1
W              0.0006   -0.1   0.1
Amplitude         0.001 -inf      inf
eta:           2
""",
"""\
U              0.2   -0.1   0.1
V              0.3   -0.1   0.1
W              0.0008   -0.1   0.1
Amplitude         0.001 -inf      inf
eta:           2
"""]

global_input_string = """\
Bkgd:          3
two_theta_0       0.001      -2.0  2.0
"""

tst_two_theta = []
tst_y = []

is_Sim_data = True #: Should be False unless simulated data 
   #: (e.g. "Jade-AL2O3-Sim.xye") is used
display_plots = True #: Only use to see sample plots
# with open(r"17_05_23_0014_NIST SRM 1976b.xye") as file:
# with open(r"16_01_07_0010_Aspirin_HighRez.xye") as file:
# with open(r"16_03_09_0015_Silver Behenate.xye") as file:
with open(r"Jade-Al2O3-Sim.xye") as file:
   for line in file.readlines():#[4:]:
      # two_thetatmp, ytmp, ztmp = line.split()
      two_thetatmp, ytmp = line.split()
      # if float(two_thetatmp) < 15:
      tst_two_theta.append(float(two_thetatmp))
      tst_y.append(float(ytmp))
tst_two_theta = np.array(tst_two_theta)
tst_y = np.array(tst_y)

def exercise_RietveldPhases():
   # RietveldPhase.fromstring(input_string)
   cifs = ["1000032.cif","1507774.cif"]
   Rt = []
   for cif,input_string in zip(cifs,input_strings):
      Rt.append(RietveldPhases(cif,input_string))
   RietveldPhases.global_params_from_string(global_input_string)

   #Testing Read-in from input_string
   assert np.isclose(Rt[0].x['values'][Rt[0].U_index], 0.0)
   assert Rt[0].x['labels'][Rt[0].U_index] == 'U'
   assert np.isclose(Rt[1].x['values'][Rt[1].U_index], 0.2)
   assert Rt[1].x['labels'][Rt[1].U_index] == 'U'
   assert np.isclose(Rt[0].x['values'][Rt[0].V_index], 0.0)
   assert Rt[0].x['labels'][Rt[0].V_index] == 'V'
   assert np.isclose(Rt[0].x['values'][Rt[0].W_index],  0.0006)
   assert np.isclose(Rt[0].x['values'][Rt[0].Amplitude_index],  0.001)
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
   assert np.all(np.isclose(Rt[0].Compute_Relative_Intensities(), \
      np.array([[  3.48114434e+00,   2.55177292e+00,   2.38025000e+00,
          2.08607570e+00,   1.74057217e+00,   1.60196736e+00,
          1.54715716e+00,   1.51527741e+00,   1.51135808e+00,
          1.40499663e+00,   1.37423798e+00,   1.33645880e+00,
          1.27588646e+00,   1.23944057e+00,   1.23455034e+00,
          1.19354167e+00,   1.19012500e+00,   1.14760201e+00,
          1.12613194e+00,   1.12452033e+00,   1.09932894e+00,
          1.08296667e+00,   1.07858495e+00,   1.04303785e+00,
          1.01795211e+00],
       [  1.29494431e+04,   3.97271668e+04,   1.94218740e+04,
          7.68995083e+04,   5.61379507e+04,   1.37688046e+05,
          3.68513103e+03,   5.40242989e+03,   1.19976722e+04,
          7.47626929e+04,   1.20456747e+05,   2.49928046e+03,
          3.81902637e+03,   4.54489735e+04,   2.30405577e+04,
          1.98365266e+03,   1.92841745e+04,   1.56760345e+04,
          1.21017630e+04,   9.26032788e+03,   2.62138368e+04,
          6.14885686e+03,   3.14784442e+04,   6.34002271e+04,
          8.10831625e+03]])))
   assert str(Rt[1].unit_cell.parameters()) == \
   """(4.1769, 4.7218, 58.3385, 89.44, 89.634, 75.854)"""
   assert np.all(np.isclose(Rt[1].Compute_Relative_Intensities()[:,0:30], \
      np.array([[  5.83352174e+01,   2.91676087e+01,   1.94450725e+01,
          1.45838043e+01,   1.16670435e+01,   9.72253623e+00,
          8.33360248e+00,   7.29190217e+00,   6.48169082e+00,
          5.83352174e+00,   5.30320158e+00,   4.86126811e+00,
          4.57845246e+00,   4.56140325e+00,   4.51721072e+00,
          4.48732441e+00,   4.35771092e+00,   4.27435487e+00,
          4.24980397e+00,   4.15574220e+00,   4.02715497e+00,
          3.99846875e+00,   3.96837552e+00,   3.90666115e+00,
          3.89836103e+00,   3.89236416e+00,   3.88901449e+00,
          3.83110776e+00,   3.75460982e+00,   3.74425360e+00],
       [  2.21750091e+04,   1.55537642e+04,   1.93755054e+04,
          1.23439317e+04,   1.53672613e+04,   8.86039694e+03,
          1.06743619e+04,   5.35494938e+03,   6.61111185e+03,
          2.82384481e+03,   3.66921150e+03,   1.44698106e+03,
          3.19586085e+04,   3.51129634e+03,   2.82920352e+04,
          2.35225575e+03,   2.87566370e+03,   1.97206178e+03,
          8.31181187e+02,   1.02073195e+03,   2.07178668e+03,
          1.24567393e+03,   2.19348434e+03,   1.10152547e+03,
          1.09517013e+03,   1.61191370e+03,   1.45277408e+03,
          6.98694572e+04,   2.30685538e+03,   8.95049478e+03]])))

   #Testing test-data read-in
   if is_Sim_data:
      np.set_printoptions(threshold=None)
      assert str(tst_two_theta) == \
      """[  5.     5.02   5.04 ...,  89.96  89.98  90.  ]"""
      assert str(tst_y) == \
      """[ 1.  1.  1. ...,  2.  2.  2.]"""

   # Testing Background_Polynomial
   RietveldPhases.x['values'] \
      [np.char.startswith(RietveldPhases.x['labels'],"Bkgd")]  \
      = np.array([1.0,2.0,3.0])
   tst_bkgd_two_theta = np.array([0,1,2.3,100.5])
   assert np.all(np.isclose(RietveldPhases.Background_Polynomial( \
      tst_bkgd_two_theta) , np.array([  1,6,21.47,   3.05027500e+04])))

   # Testing Eta_Polynomial
   Rt0_eta_mask = np.isin(np.array(range(0,len(RietveldPhases.x))), \
      np.array(range(Rt[0].eta_0_index,Rt[0].eta_0_index+Rt[0].eta_rank)))
   RietveldPhases.x['values'][Rt0_eta_mask] = np.array([0.5,0.005])
   tst_eta_two_theta = np.array([0,1,2.3,100.5])
   assert np.all(np.isclose(Rt[0].eta_Polynomial( \
      tst_bkgd_two_theta) , np.array([ 0.5, 0.505, 0.5115, 1.00249999])))

   # Testing LP_Intensity_Scaling
   assert np.isclose(Rt[0].LP_Intensity_Scaling(20.0,0.45),14.2674396526)

   #Testing Compile_Weighted_Peak_Intensities
   Rt[0].Compile_Weighted_Peak_Intensities()
   assert np.all(np.isclose(Rt[0].two_theta_peaks[0:10], np.array(
      [[ 25.56751716],
       [ 35.13879293],
       [ 37.76313851],
       [ 43.33853666],
       [ 52.53263976],
       [ 57.47974186],
       [ 59.71809269],
       [ 61.10698338],
       [ 61.28247987],
       [ 66.49316001]])))
   assert np.all(np.isclose(Rt[0].weighted_intensities[0:10], np.array(
      [[ 245948.79616911],
       [ 381572.99664677],
       [ 159248.63470755],
       [ 463968.43986271],
       [ 218966.12149661],
       [ 437750.54517789],
       [  10750.80346893],
       [  14972.33441136],
       [  33039.19755184],
       [ 172366.48175023]])))

   CU_wavelength = wavelengths.characteristic("CU").as_angstrom()
   d_min = CU_wavelength/2/np.sin(np.pi/360*tst_two_theta[-1])
   # print "two_theta_max: " + str(tst_two_theta[-1])
   # print "d-min: "+ str(d_min)
   
   for RV in Rt:
      RV.Compute_Relative_Intensities(d_min=d_min)
      RV.Compile_Weighted_Peak_Intensities()

   #Testing PseudoVoigtProfile
   if display_plots:
      for RV in Rt:
         #Select a random peak:
         if len(RV.two_theta_peaks) != 0:
            rnd_index = randrange(0,len(RV.two_theta_peaks),1)
            tst_two_theta_peak = RV.two_theta_peaks[rnd_index]
            tst_weighted_intensity = RV.weighted_intensities[rnd_index]
            # print rnd_index, tst_two_theta_peak, tst_weighted_intensity
            delta_theta = 5
            # mask = np.ones(len(tst_two_theta),dtype=bool)
            mask = np.abs(tst_two_theta-tst_two_theta_peak) < delta_theta
            RV.showPVProfilePlot("Test",rnd_index,tst_two_theta[mask], \
               tst_y[mask], autohide=False)

   #Testing Refinery
   RR = RietveldRefinery(Rt,tst_two_theta,tst_y)
   if display_plots:
      RR.show_multiplot("Sum of Phases", \
         two_theta_roi=30, \
         delta_theta=10, \
         autohide=False)
   t0 = time.time()
   RR.minimize()
   t1 = time.time()
   print "Time elapse: " + str(t1-t0)
   if display_plots:
      print RietveldPhases.x['labels']
      print RietveldPhases.x['values']
      RR.show_multiplot("Sum of Phases", \
         two_theta_roi=30, \
         delta_theta=10, \
         autohide=False)

def run():
   exercise_RietveldPhases()
   print "OK"

if (__name__ == "__main__"):
   run()