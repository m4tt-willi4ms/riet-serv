from __future__ import division
# from scitbx.array_family import flex
# import scitbx.math
# from libtbx.test_utils import approx_equal, eps_eq, Exception_expected
# import sys
import numpy as np
import time
import matplotlib.pyplot as plt
from random import randrange

import sys, os
sys.path.append(os.path.abspath("../.."))

from RietveldPhases import RietveldPhases
from RietveldRefinery import RietveldRefinery

from scitbx import lbfgsb
from cctbx.eltbx import wavelengths
from libtbx import test_utils
import libtbx.load_env

input_strings = ["""\
U              0.0   -0.1   0.1
V              0.0   -0.1   0.1
W              0.0006   -0.1   0.1
Amplitude         0.1 0      inf
eta:           1
""",
"""\
U              0.2   -0.1   0.1
V              0.3   -0.1   0.1
W              0.0008   -0.1   0.1
Amplitude         0.000001 0      inf
eta:           1
""",
"""\
U              0.0   -0.1   0.1
V              0.0   -0.1   0.1
W              0.0006   -0.1   0.1
Amplitude         0.1 0      inf
eta:           1
""",
"""\
U              0.0   -0.1   0.1
V              0.0   -0.1   0.1
W              0.0006   -0.1   0.1
Amplitude         0.1 0      inf
eta:           1
""",
"""\
U              0.0   -0.1   0.1
V              0.0   -0.1   0.1
W              0.0006   -0.1   0.1
Amplitude         0.1 0      inf
eta:           1
""",
"""\
U              0.0   -0.1   0.1
V              0.0   -0.1   0.1
W              0.0006   -0.1   0.1
Amplitude         0.1 0      inf
eta:           1
""",
"""\
U              0.0   -0.1   0.1
V              0.0   -0.1   0.1
W              0.0006   -0.1   0.1
Amplitude         0.1 0      inf
eta:           1
""",
"""\
U              0.0   -0.1   0.1
V              0.0   -0.1   0.1
W              0.0006   -0.1   0.1
Amplitude         0.1 0      inf
eta:           1
"""]

global_input_string = """\
Bkgd:          3
two_theta_0       0.      -0.5  0.5
"""

bkgd_minimizer_input_string = """\
approx_grad True
factr       1e2
iprint      1
maxiter     150
m           7
pgtol       1e-5
epsilon     1e-5
"""

minimizer_input_string = """\
approx_grad True
factr       1e3
iprint      1000
maxiter     100
m           7
pgtol       1e-5
epsilon     1e-6
"""

tst_two_theta = []
tst_y = []

is_Sim_data = True #: Should be False unless simulated data 
   #: (e.g. "Jade-AL2O3-Sim.xye") is used
display_plots = True #: Only use to see sample plots
# with open(r"17_05_23_0014_NIST SRM 1976b.xye") as file:
# with open(r"16_01_07_0010_Aspirin_HighRez.xye") as file:
# with open(r"16_03_09_0015_Silver Behenate.xye") as file:
# os.path.dirname(__file__) + r
with open("cement_15_03_11_0028.xye") as file:
   for line in file.readlines()[1:]:
      two_thetatmp, ytmp, ztmp = line.split()
      # two_thetatmp, ytmp = line.split()
      # if float(two_thetatmp) < 15:
      tst_two_theta.append(float(two_thetatmp))
      tst_y.append(float(ytmp))
tst_two_theta = np.array(tst_two_theta)
# mask = np.ones(len(tst_two_theta),dtype=bool) 
mask = np.logical_and(tst_two_theta >25,tst_two_theta < 45)
tst_two_theta = tst_two_theta[mask]
tst_y = np.array(tst_y)[mask]

def exercise_Rietveld_Refinery_Cement():
   # RietveldPhase.fromstring(input_string) 
   cifs = ["1540705-Alite.cif", \
      "9012789-Belite.cif", \
      "1200009-Ferrite.cif", \
      "1011094-FreeLime.cif", \
      "1000039-AluminateCubic.cif", \
      "9014308-AluminateOrtho.cif", \
      "1000053-Periclase.cif", \
      "9007569-Arcanite.cif"]
   full_cifs = list(map( \
      lambda x: os.path.dirname(__file__) + "\\" + x ,cifs))
   Rt = []


   CU_wavelength = wavelengths.characteristic("CU").as_angstrom()
   d_min = CU_wavelength/2/np.sin(np.pi/360*tst_two_theta[-1])
   d_max = CU_wavelength/2/np.sin(np.pi/360*tst_two_theta[0])
   # print "two_theta_max: " + str(tst_two_theta[-1])
   # print "d-min: "+ str(d_min)
   tst_y_max = np.amax(tst_y)/len(cifs)

   RietveldPhases.global_params_from_string(global_input_string)

   for cif, input_string in zip(cifs,input_strings):
      Rt.append(RietveldPhases(cif,input_string,d_min,d_max, \
         tst_y_max,delta_theta=5.0,Intensity_Cutoff = 0.005))

   #First fit the background
   RR = RietveldRefinery(Rt,tst_two_theta,tst_y,bkgd_minimizer_input_string, \
      use_bkgd_mask=True,bkgd_delta_theta=0.1)
   # RR.display(RR.minimize_Bkgd_0)
   RR.display(RR.minimize_Bkgd)

   #Now use the full dataset
   RR = RietveldRefinery(Rt,tst_two_theta,tst_y,minimizer_input_string)

   RR.show_multiplot("Before: ", \
      two_theta_roi=32.5, \
      delta_theta=3, \
      autohide=False)

   RR.display(RR.minimize_Amplitude_Offset)
   RR.display(RR.minimize_Amplitude_Bkgd_Offset)
   RR.display(RR.minimize_Amplitude_Bkgd_Offset_W)
   # RR.display(RR.minimize_only_Alite)
   # RR.display(RR.minimize_All)


def run():
   exercise_Rietveld_Refinery_Cement()
   print "OK"

if (__name__ == "__main__"):
   run()