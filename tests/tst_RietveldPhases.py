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
from random import randint
from libtbx import test_utils
import libtbx.load_env

import sys, os
sys.path.append(os.path.abspath(".."))

from RietveldPhases import RietveldPhases
from RietveldRefinery import RietveldRefinery

input_strings = ["""\
U              0.0   0.0   0.0
V              0.0   0.0   0.0
W              0.0006   0.0   0.0
Amplitude         0.001 0.0      0.0
eta:           2
""",
"""\
U              0.2   0.0   0.0
V              0.3   0.0   0.0
W              0.0008   0.0   0.0
Amplitude         0.05 0.0      0.0
eta:           4
"""]

global_input_string = """\
Bkgd:          3
two_theta_0       0.0      -2.0  2.0
"""

tst_two_theta = []
tst_y = []
# with open(r"17_05_23_0014_NIST SRM 1976b.xye") as file:
# with open(r"16_01_07_0010_Aspirin_HighRez.xye") as file:
# with open(r"16_03_09_0015_Silver Behenate with Baffle and Anti_Scatter \
   # Slit_highrez.xye") as file:
is_Sim_data = True # Should be False unless "Jade-AL2O3-Sim.xye" is used
display_plots = True # Only use to see sample plots
with open(r"Jade-Al2O3-Sim.xye") as file:
   for line in file.readlines():
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
   RietveldPhases.global_params_from_string(global_input_string)
   for cif,input_string in zip(cifs,input_strings):
      Rt.append(RietveldPhases(cif,input_string))

   #Testing Read-in from input_string
   assert np.isclose(Rt[0].x['values'][Rt[0].U_index], 0.0)
   assert Rt[0].x['labels'][Rt[0].U_index] == 'U'
   assert np.isclose(Rt[1].x['values'][Rt[1].U_index], 0.2)
   assert Rt[1].x['labels'][Rt[1].U_index] == 'U'
   assert np.isclose(Rt[0].x['values'][Rt[0].V_index], 0.0)
   assert Rt[0].x['labels'][Rt[0].V_index] == 'V'
   assert np.isclose(Rt[0].x['values'][Rt[0].W_index],  0.0006)
   assert np.isclose(Rt[0].x['values'][Rt[0].Amplitude_index],  0.001)
   assert RietveldPhases.two_theta_0 == 0.0
   # assert len(RietveldPhases.Bkgd) == 3
   # assert np.array_equal(RietveldPhases.Bkgd, np.array([0.,0.,0.]))

   # Testing iotbx output from .cif card
   np.set_printoptions(threshold=None)
   assert str(Rt[0].unit_cell.parameters()) == \
   """(4.7605, 4.7605, 12.995599999999998, 90.0, 90.0, 120.00000000000001)"""
   assert np.all(np.isclose(Rt[0].Compute_Relative_Intensities(), \
      np.array([[  3.48114434e+00,   2.55177292e+00,   2.38025000e+00,
          2.16593333e+00,   2.08607570e+00,   1.96485460e+00,
          1.74057217e+00,   1.60196736e+00,   1.54715716e+00,
          1.51527741e+00,   1.51135808e+00,   1.40499663e+00,
          1.37423798e+00,   1.33645880e+00,   1.27588646e+00,
          1.23944057e+00,   1.23455034e+00,   1.19354167e+00,
          1.19012500e+00,   1.16038145e+00,   1.16038145e+00,
          1.14760201e+00,   1.13903464e+00,   1.12613194e+00,
          1.12452033e+00,   1.09932894e+00,   1.08296667e+00,
          1.07858495e+00,   1.04662972e+00,   1.04303785e+00,
          1.01795211e+00],
       [  1.29494431e+04,   3.97271668e+04,   1.94218740e+04,
          6.41107374e+02,   7.68995083e+04,   1.13236117e+03,
          5.61379507e+04,   1.37688046e+05,   3.68513103e+03,
          5.40242989e+03,   1.19976722e+04,   7.47626929e+04,
          1.20456747e+05,   2.49928046e+03,   3.81902637e+03,
          4.54489735e+04,   2.30405577e+04,   1.98365266e+03,
          1.92841745e+04,   1.36554203e+03,   1.37145458e+03,
          1.56760345e+04,   9.05255830e+02,   1.21017630e+04,
          9.26032788e+03,   2.62138368e+04,   6.14885686e+03,
          3.14784442e+04,   6.91894468e+02,   6.34002271e+04,
          8.10831625e+03]])))
   assert str(Rt[1].unit_cell.parameters()) == \
   """(4.1769, 4.7218, 58.3385, 89.44, 89.634, 75.854)"""
   assert str(Rt[1].Compute_Relative_Intensities()) == \
   """[[  5.83352174e+01   2.91676087e+01   1.94450725e+01 ...,   1.00030333e+00
    1.00019992e+00   1.00015760e+00]
 [  2.21750091e+04   1.55537642e+04   1.93755054e+04 ...,   2.75225955e+01
    2.64709614e+02   1.93820027e+02]]"""

   #Testing test-data read-in
   if is_Sim_data:
      np.set_printoptions(threshold=None)
      assert str(tst_two_theta) == \
      """[  5.     5.02   5.04 ...,  89.96  89.98  90.  ]"""
      assert str(tst_y) == \
      """[ 1.  1.  1. ...,  2.  2.  2.]"""

   # Testing Background_Polynomial
   RietveldPhases.Bkgd = np.array([1.0,2.0,3.0])
   tst_bkgd_two_theta = np.array([0,1,2.3,100.5])
   assert np.all(np.isclose(RietveldPhases.Background_Polynomial( \
      tst_bkgd_two_theta) , np.array([  1,6,21.47,   3.05027500e+04])))

   # Testing LP_Intensity_Scaling
   assert np.isclose(Rt[0].LP_Intensity_Scaling(20.0,0.45),14.2674396526)

   #Testing Compile_Weighted_Peak_Intensities
   Rt[0].Compile_Weighted_Peak_Intensities()
   assert np.all(np.isclose(Rt[0].two_theta_peaks[0:10], np.array(
      [[ 25.56751716],
       [ 35.13879293],
       [ 37.76313851],
       [ 41.66463806],
       [ 43.33853666],
       [ 46.16161879],
       [ 52.53263976],
       [ 57.47974186],
       [ 59.71809269],
       [ 61.10698338]])))
   assert np.all(np.isclose(Rt[0].weighted_intensities[0:10], np.array(
      [[ 245948.79616911],
       [ 381572.99664677],
       [ 159248.63470755],
       [   4225.19160504],
       [ 463968.43986271],
       [   5925.64874038],
       [ 218966.12149661],
       [ 437750.54517789],
       [  10750.80346893],
       [  14972.33441136]])))

   CU_wavelength = wavelengths.characteristic("CU").as_angstrom()
   d_min = CU_wavelength/2/np.sin(np.pi/360*tst_two_theta[-1])
   Rt[0].Compute_Relative_Intensities(d_min=d_min)
   Rt[0].Compile_Weighted_Peak_Intensities()

   #Testing PseudoVoigtProfile
   if display_plots:
      #Select a random peak:
      rnd_index = randint(0,len(Rt[0].two_theta_peaks)-1)
      tst_two_theta_peak = Rt[0].two_theta_peaks[rnd_index]
      tst_weighted_intensity = Rt[0].weighted_intensities[rnd_index]
      # print rnd_index, tst_two_theta_peak, tst_weighted_intensity
      delta_theta = 5
      # mask = np.ones(len(tst_two_theta),dtype=bool)
      mask = np.abs(tst_two_theta-tst_two_theta_peak) < delta_theta
      showPVProfilePlot("Test",Rt[0],rnd_index,tst_two_theta[mask],tst_y[mask] \
         ,np.zeros(1))

   #Testing Refinery
   RR = RietveldRefinery(Rt)
   print "here"
   print RietveldPhases.x_global
   print Rt[0].x


def showPVProfilePlot(plottitle,Rt,index,two_theta,y,Peak_Intensity, delta_theta=0.5):
   plt.ion()
   fig = plt.figure(figsize=(12, 8))
   # plt.subplot(3,1,1)
   plt.scatter(two_theta,y,label='Data',s=1, color='red')
   plt.title(plottitle)

   plt.plot(two_theta,np.sum(Rt.PseudoVoigtProfile(two_theta),axis=0), \
      label=r'Total $I_{\rm calc}$')
   plt.plot(two_theta,Rt.PseudoVoigtProfile(two_theta)[index], \
      label=r'$I_{\rm calc}$')

   plt.legend(bbox_to_anchor=(.8,.7))
   plt.ylabel(r"$I$")

   # fig, ax = plt.subplots()
   plt.show()
   fig.canvas.flush_events()
   time.sleep(1)
   plt.close('all')

# def showplot(plottitle,two_theta,x,y,Rel_Peak_Intensity,delta_theta):
#    plt.figure(figsize=(12, 8))
#    plt.subplot(3,1,1)
#    # for i in xrange(0,len(two_theta)):
#    #    if delta_theta[i]:
#    #       color = 'blue'
#    #    else: color = 'green'
#    #    plt.scatter(two_theta[i],y[i], color=color,s=1)
#    plt.scatter(two_theta,y,label='Data',s=1, color='red')
#    plt.title(plottitle)
#    # plt.axis([20,60,0,1500])

#    plt.plot(two_theta,Profile_Calc(x,two_theta,Rel_Peak_Intensity,delta_theta), label=r'$I_{\rm calc}$')
#    plt.legend(bbox_to_anchor=(.8,.7))
#    plt.ylabel(r"$I$")

#    plt.subplot(3,1,2)
#    # for i in xrange(0,len(two_theta)):
#    #    if delta_theta[i]:
#    #       color = 'blue'
#    #    else: color = 'green'
#    #    plt.scatter(two_theta[i],y[i], color=color,s=1)
#    pltmask = np.logical_and(two_theta >Rel_Peak_Intensity[0,0]-1, two_theta < Rel_Peak_Intensity[0,0]+1)
#    plt.scatter(two_theta[pltmask],y[pltmask],label='Data',s=1, color='red')
#    # plt.title(r"Profile: $Al_2 O_3$")
#    # plt.axis([20,60,0,1500])

#    plt.plot(two_theta[pltmask],Profile_Calc(x,two_theta,Rel_Peak_Intensity,delta_theta)[pltmask], label=r'$I_{\rm calc}$')
#    plt.legend(bbox_to_anchor=(.8,.7))
#    plt.ylabel(r"$I$")

#    plt.subplot(3,1,3)
#    zf = flex.double(len(two_theta))
#    y_calc = Profile_Calc(x,two_theta,Rel_Peak_Intensity,delta_theta)
#    for i in xrange(len(two_theta)):
#       zf[i] = 1/y[i]*(y[i]-y_calc[i])**2
#    plt.scatter(two_theta,zf)
#    plt.ylabel(r"$\frac{1}{I} \, (I-I_{\rm calc})^2$")
#    plt.xlabel(r'$2\,\theta$')

#    plt.show()#block=False)

def run():
   exercise_RietveldPhases()
   print "OK"

if (__name__ == "__main__"):
   run()