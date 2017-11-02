from __future__ import division
import numpy as np
import time
import matplotlib.pyplot as plt
from random import randrange
import profile, pstats

import sys, os
sys.path.append(os.path.abspath(".."))

from RietveldPhases import RietveldPhases
from RietveldRefinery import RietveldRefinery

from cctbx.eltbx import wavelengths


input_strings = ["""\
U              0.0   -0.1   0.1
V              0.0   -0.1   0.1
W              0.0006   -0.1   0.1
Amplitude         1 0      inf
eta:           2
""",
"""\
U              0.2   -0.1   0.1
V              0.3   -0.1   0.1
W              0.0008   -0.1   0.1
Amplitude         1 0      inf
eta:           2
"""]

global_input_string = """\
Bkgd:          3
two_theta_0       0.001      -2.0  2.0
"""

minimizer_input_string = """\
factr       1e4
maxiter     100
iprint      1
m           10
pgtol       1e-5
epsilon     1e-8
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
tst_y = 0.01*np.array(tst_y)

def exercise_RietveldPhases():
   # RietveldPhase.fromstring(input_string)
   cifs = ["1000032.cif","1507774.cif"]
   Rt = []

   CU_wavelength = wavelengths.characteristic("CU").as_angstrom()
   d_min = CU_wavelength/2/np.sin(np.pi/360*tst_two_theta[-1])
   d_max = CU_wavelength/2/np.sin(np.pi/360*tst_two_theta[0])
   # print "two_theta_max: " + str(tst_two_theta[-1])
   # print "d-min: "+ str(d_min)

   RietveldPhases.global_params_from_string(global_input_string,
      tst_two_theta,tst_y)
   for cif,input_string in zip(cifs,input_strings):
      Rt.append(RietveldPhases(cif,input_string,d_min,d_max, 
         delta_theta = 2.0,Intensity_Cutoff = 0.005))

   #Testing Read-in from input_string
   assert np.isclose(Rt[0].x['values'][Rt[0].U_index], 0.0)
   assert Rt[0].x['labels'][Rt[0].U_index] == 'U'
   assert np.isclose(Rt[1].x['values'][Rt[1].U_index], 0.2)
   assert Rt[1].x['labels'][Rt[1].U_index] == 'U'
   assert np.isclose(Rt[0].x['values'][Rt[0].V_index], 0.0)
   assert Rt[0].x['labels'][Rt[0].V_index] == 'V'
   assert np.isclose(Rt[0].x['values'][Rt[0].W_index],  0.0006)
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

   #Testing Refinery
   RR = RietveldRefinery(Rt,minimizer_input_string,store_intermediate_state=True)

   RR.display(RR.minimize_Amplitude_Offset)
   RR.display(RR.minimize_Amplitude_Offset_W)
   RR.display(RR.minimize_All)


def run():
   exercise_RietveldPhases()
   print "OK"

if (__name__ == "__main__"):
   # pr = profile.Profile()
   # pr.enable()
   run()
   # profile.run('run(); print')
   # pr.disable()
   # s = StringIO.StringIO()
   # sortby = 'cumulative'
   # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
   # ps.print_stats()
   # print s.getvalue()