from __future__ import division
import numpy as np
import time
import matplotlib.pyplot as plt
from random import randrange

import sys, os

from src.RietveldPhases import RietveldPhases
from src.RietveldRefinery import RietveldRefinery, RietveldPlot

from cctbx.eltbx import wavelengths
from libtbx import test_utils
import libtbx.load_env

RietveldPhases.set_profile(r".//data//profiles//Jade-Al2O3-Sim.xye",number_of_columns=2)

cifs = ["1000032.cif","1507774.cif"]

test_phase_list = []
test_phase_list.append(RietveldPhases(r".//data//cifs//" + cifs[0],
          delta_theta=2.0,intensity_cutoff=0.005))
test_plot = RietveldPlot()

test_refinery = RietveldRefinery(test_phase_list,test_plot)

def test_read_in():
   tr_dict = test_refinery.__dict__
   assert 'phase_list' in tr_dict
   assert 'x' in tr_dict
   assert 'mask' in tr_dict

def exercise_Rietveld_Refinery_SinglePhase():
   RietveldPhases.global_params_from_string(global_input_string,
      tst_two_theta,tst_y)



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
      Rt.append(RietveldPhases(r"..//data//cifs//" + cif, d_min,d_max,
         input_string_or_file_name=input_string,
         delta_theta=2.0,Intensity_Cutoff=0.005))

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