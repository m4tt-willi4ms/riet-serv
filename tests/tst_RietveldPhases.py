from __future__ import division
from scitbx import lbfgsb
# from scitbx.array_family import flex
# import scitbx.math
# from libtbx.test_utils import approx_equal, eps_eq, Exception_expected
# import sys
import numpy as np

import sys, os
sys.path.append(os.path.abspath(".."))

from RietveldPhases import RietveldPhases

input_string = """\
two_theta_0       0.0      -2.0  2.0
U              0.0   0.0   0.0
V              0.0   0.0   0.0
W              0.0006   0.0   0.0
Amplitude         0.001 0.0      0.0
eta:           2
Bkgd:          3
"""

def exercise_RietveldPhases():
   # RietveldPhase.fromstring(input_string)
   cifs = ["1000032.cif","1507774.cif"]
   Rt = []
   RietveldPhases.fromstring(input_string)
   for cif in cifs:
      Rt.append(RietveldPhases(cif))
   # Rt = RietveldPhases("test")
   # Rt.fromfile("RefinementParameters.txt")
   assert RietveldPhases.U == 0.0
   assert RietveldPhases.V == 0.0
   assert RietveldPhases.W == 0.0006
   assert RietveldPhases.Amplitude == 0.001
   assert RietveldPhases.two_theta_0 == 0.0
   assert len(RietveldPhases.Bkgd) == 3
   assert np.array_equal(RietveldPhases.Bkgd, np.array([0.,0.,0.]))
   # Rt.load_cif("1000032.cif")
   assert str(Rt[0].unit_cell.parameters()) == \
   """(4.7605, 4.7605, 12.995599999999998, 90.0, 90.0, 120.00000000000001)"""
   assert str(Rt[0].Compute_Relative_Intensities()) == \
   """[[  3.48114434e+00   2.55177292e+00   2.38025000e+00   2.16593333e+00
    2.08607570e+00   1.96485460e+00   1.74057217e+00   1.60196736e+00
    1.54715716e+00   1.51527741e+00   1.51135808e+00   1.40499663e+00
    1.37423798e+00   1.33645880e+00   1.27588646e+00   1.23944057e+00
    1.23455034e+00   1.19354167e+00   1.19012500e+00   1.16038145e+00
    1.16038145e+00   1.14760201e+00   1.13903464e+00   1.12613194e+00
    1.12452033e+00   1.09932894e+00   1.08296667e+00   1.07858495e+00
    1.04662972e+00   1.04303785e+00   1.01795211e+00]
 [  1.29494431e+04   3.97271668e+04   1.94218740e+04   6.41107374e+02
    7.68995083e+04   1.13236117e+03   5.61379507e+04   1.37688046e+05
    3.68513103e+03   5.40242989e+03   1.19976722e+04   7.47626929e+04
    1.20456747e+05   2.49928046e+03   3.81902637e+03   4.54489735e+04
    2.30405577e+04   1.98365266e+03   1.92841745e+04   1.36554203e+03
    1.37145458e+03   1.56760345e+04   9.05255830e+02   1.21017630e+04
    9.26032788e+03   2.62138368e+04   6.14885686e+03   3.14784442e+04
    6.91894468e+02   6.34002271e+04   8.10831625e+03]]"""
   assert str(Rt[1].unit_cell.parameters()) == \
   """(4.1769, 4.7218, 58.3385, 89.44, 89.634, 75.854)"""
   assert str(Rt[1].Compute_Relative_Intensities()) == \
   """[[  5.83352174e+01   2.91676087e+01   1.94450725e+01 ...,   1.00030333e+00
    1.00019992e+00   1.00015760e+00]
 [  2.21750091e+04   1.55537642e+04   1.93755054e+04 ...,   2.75225955e+01
    2.64709614e+02   1.93820027e+02]]"""

def run():
   exercise_RietveldPhases()
   print "OK"

if (__name__ == "__main__"):
   run()
