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
eta:           2
""",
"""\
U              0.2   -0.1   0.1
V              0.3   -0.1   0.1
W              0.0008   -0.1   0.1
Amplitude         0.000001 0      inf
eta:           4
"""]

global_input_string = """\
Bkgd:          3
two_theta_0       0.      -2.0  2.0
"""

minimizer_input_string = """\
approx_grad True
factr       1e10
iprint      10
m           5
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
tst_y = np.array(tst_y)

def exercise_Rietveld_Refinery_SinglePhase():
   # RietveldPhase.fromstring(input_string)
   cifs = ["1000032.cif","1507774.cif"]
   Rt = []
   Rt.append(RietveldPhases(cifs[0],input_strings[0]))
   RietveldPhases.global_params_from_string(global_input_string)

   CU_wavelength = wavelengths.characteristic("CU").as_angstrom()
   d_min = CU_wavelength/2/np.sin(np.pi/360*tst_two_theta[-1])
   # print "two_theta_max: " + str(tst_two_theta[-1])
   # print "d-min: "+ str(d_min)
   
   for RV in Rt:
      RV.Compute_Relative_Intensities(d_min=d_min)
      RV.Compile_Weighted_Peak_Intensities()

   RR = RietveldRefinery(Rt,tst_two_theta,tst_y,minimizer_input_string)
   # RR.minimize_Amplitude_and_Offset(display_plots = display_plots)
   RR.display(RR.minimize_Amplitude_Offset)
   RR.display(RR.minimize_All)

def exercise_Rietveld_Refinery_Multiphase():
   pass

def run():
   exercise_Rietveld_Refinery_SinglePhase()
   exercise_Rietveld_Refinery_Multiphase()
   print "OK"

if (__name__ == "__main__"):
   run()