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
Amplitude         0.1 0      inf
eta:           1
""",
"""\
U              0.2   -0.1   0.1
V              0.3   -0.1   0.1
W              0.0008   -0.1   0.1
Amplitude         0.000001 0      inf
eta:           4
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
two_theta_0       0.      -2.0  2.0
"""

tst_two_theta = []
tst_y = []

is_Sim_data = True #: Should be False unless simulated data 
   #: (e.g. "Jade-AL2O3-Sim.xye") is used
display_plots = True #: Only use to see sample plots
# with open(r"17_05_23_0014_NIST SRM 1976b.xye") as file:
# with open(r"16_01_07_0010_Aspirin_HighRez.xye") as file:
# with open(r"16_03_09_0015_Silver Behenate.xye") as file:
with open(r"cement_15_03_11_0028.xye") as file:
   for line in file.readlines()[1:]:
      two_thetatmp, ytmp, ztmp = line.split()
      # two_thetatmp, ytmp = line.split()
      # if float(two_thetatmp) < 15:
      tst_two_theta.append(float(two_thetatmp))
      tst_y.append(float(ytmp))
tst_two_theta = np.array(tst_two_theta)
tst_y = np.array(tst_y)

def exercise_Rietveld_Refinery_SinglePhase():
   # RietveldPhase.fromstring(input_string)
   cifs = ["CSD#1001Alite.cif", \
      "CSD#1002Belite.cif", \
      "CSD#1003Ferrite.cif", \
      "CSD#1004FreeLime.cif", \
      "CSD#1005AluminateCubic.cif", \
      "CSD#1006AluminateOrtho.cif", \
      "CSD#1007Periclase.cif", \
      "CSD#1008Arcanite.cif"]
   Rt = []
   for cif, input_string in zip(cifs,input_strings):
      Rt.append(RietveldPhases(cif,input_string))
   RietveldPhases.global_params_from_string(global_input_string)

   CU_wavelength = wavelengths.characteristic("CU").as_angstrom()
   d_min = CU_wavelength/2/np.sin(np.pi/360*tst_two_theta[-1])
   # print "two_theta_max: " + str(tst_two_theta[-1])
   # print "d-min: "+ str(d_min)
   
   for RV in Rt:
      RV.Compute_Relative_Intensities(d_min=d_min)
      RV.Compile_Weighted_Peak_Intensities()

   RR = RietveldRefinery(Rt,tst_two_theta,tst_y)
   t0 = time.time()
   RR.minimize_Amplitude_and_Offset(display_plots = display_plots)
   t1 = time.time()
   print "Time elapsed: " + str(t1-t0)

   t0 = time.time()
   RR.display(RR.minimize_All)
   t1 = time.time()
   print "Time elapsed: " + str(t1-t0)

def exercise_Rietveld_Refinery_Multiphase():
   pass

def run():
   exercise_Rietveld_Refinery_SinglePhase()
   exercise_Rietveld_Refinery_Multiphase()
   print "OK"

if (__name__ == "__main__"):
   run()