from __future__ import division
import numpy as np
import time
import matplotlib.pyplot as plt
from random import randrange

import sys, os
sys.path.append(os.path.abspath(".."))

from src.RietveldPhases import RietveldPhases
from src.RietveldRefinery import RietveldRefinery

from cctbx.eltbx import wavelengths
from libtbx import test_utils
import libtbx.load_env

input_strings = ["""\
U              0.0   -0.1   0.1
V              0.0   -0.1   0.1
W              0.0006   -0.1   0.1
Amplitude         0.1 0      inf
eta:           4
""",
"""\
U              0.2   -0.1   0.1
V              0.3   -0.1   0.1
W              0.0008   -0.1   0.1
Amplitude         0.000001 0      inf
eta:           2
"""]

global_input_string = """\
Bkgd:          3
two_theta_0       0.      -2.0  2.0
"""

minimizer_input_string = """\
approx_grad True
factr       1e9
maxiter     150
iprint      -1
m           10
pgtol       1e-5
epsilon     1e-10
"""

tst_two_theta = []
tst_y = []

is_Sim_data = True #: Should be False unless simulated data 
   #: (e.g. "Jade-AL2O3-Sim.xye") is used
display_plots = True #: Only use to see sample plots
# with open(r"17_05_23_0014_NIST SRM 1976b.xye") as file:
# with open(r"16_01_07_0010_Aspirin_HighRez.xye") as file:
# with open(r"16_03_09_0015_Silver Behenate.xye") as file:
with open(r"..//data//profiles//Jade-Al2O3-Sim.xye") as file:
   for line in file.readlines():#[4:]:
      # two_thetatmp, ytmp, ztmp = line.split()
      two_thetatmp, ytmp = line.split()
      # if float(two_thetatmp) < 15:
      tst_two_theta.append(float(two_thetatmp))
      tst_y.append(float(ytmp))
tst_two_theta = np.array(tst_two_theta)
tst_y = np.array(tst_y)
cifs = ["1000032.cif","1507774.cif"]

CU_wavelength = wavelengths.characteristic("CU").as_angstrom()
d_min = CU_wavelength/2/np.sin(np.pi/360*tst_two_theta[-1])
d_max = CU_wavelength/2/np.sin(np.pi/360*tst_two_theta[0])

def exercise_Rietveld_Refinery_SinglePhase():
   RietveldPhases.global_params_from_string(global_input_string,
      tst_two_theta,tst_y)

   Rt = []
   Rt.append(RietveldPhases(r"..//data//cifs//" + cifs[0],
      input_strings[0],d_min,d_max, delta_theta=2.0,Intensity_Cutoff=0.005))

   RR = RietveldRefinery(Rt,minimizer_input_string,
      store_intermediate_state=True)
   RR.display(RR.minimize_Amplitude_Offset)
   RR.display(RR.minimize_All)
   RietveldPhases.empty_x()

def exercise_Rietveld_Refinery_Multiphase():
   RietveldPhases.global_params_from_string(global_input_string,
      tst_two_theta,tst_y)

   Rt = []
   for cif,input_string in zip(cifs,input_strings):
      Rt.append(RietveldPhases(r"..//data//cifs//" + cif,
         input_string,d_min,d_max, delta_theta=2.0,Intensity_Cutoff=0.005))

   RR = RietveldRefinery(Rt,minimizer_input_string,
      store_intermediate_state=True)
   RR.display(RR.minimize_Amplitude_Offset)
   RR.display(RR.minimize_All)
   RietveldPhases.empty_x()

def run():
   exercise_Rietveld_Refinery_SinglePhase()
   exercise_Rietveld_Refinery_Multiphase()
   print "OK"

if (__name__ == "__main__"):
   run()