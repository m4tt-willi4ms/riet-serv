from __future__ import division
import numpy as np
import time
import matplotlib.pyplot as plt

from src.RietveldPhases import RietveldPhases
from src.RietveldRefinery import RietveldRefinery

def exercise_Rietveld_Refinery_Cement():

   RietveldPhases.set_profile(r"data//profiles//17_11_15_0004_CEMI425R_d6.xye",
      min_two_theta=22)
   # RietveldPhases.set_profile(r"data//profiles//cement_15_03_11_0028.xye",
   #    min_two_theta=25)
   RietveldPhases.set_bkgd_order(3)

   cifs = [
      "1540705-Alite.cif", 
      "1000039-AluminateCubic.cif", 
      "9014308-AluminateOrtho.cif", 
      "9004096-anhydrite.cif",
      "9007569-Arcanite.cif",
      "9005521-bassanite.cif",
      "9012789-Belite.cif", 
      "9009667-calcite.cif",
      "1200009-Ferrite.cif", 
      "1011094-FreeLime.cif", 
      "1000053-Periclase.cif", 
      "9000113-portlandite.cif",
      ]

   print "cifs: \n" 
   for p in cifs:
      print p

   Rt = []
   for cif in cifs:
   #    tt0 = time.time()
      Rt.append(RietveldPhases(r"data//cifs//Cement/"+cif,
         delta_theta=0.5,Intensity_Cutoff = 0.005))
   #    tt1 = time.time()
   #    print "TIME TO READ IN: " +str(tt1-tt0) + " seconds"

   # First fit the background
   RR = RietveldRefinery(Rt,bkgd_refine=True,factr=1e3,
      store_intermediate_state=False, show_plots=True)
   RR.display(RR.minimize_bkgd)

   #Now use the full dataset
   RR = RietveldRefinery(Rt, 
      factr=1e9,store_intermediate_state=False, show_plots=True)
   # RR.display(RR.minimize_Bkgd)
   # RR.display(RR.minimize_Bkgd_Offset)
   # RR.display(RR.minimize_Amplitude)
   # RR.display(RR.minimize_Amplitude)
   # RR.display(RR.minimize_Amplitude_Offset)
   # # RR.display(RR.minimize_Amplitude_Offset_unit_cell)
   # RR.display(RR.minimize_unit_cell)
   # # RR.display(RR.minimize_First_n_Phases)
   RR.display(RR.minimize_Amplitude_Bkgd_Offset)
   # # RR.display(RR.minimize_First_n_Phases,n=3)
   # RR.display(RR.minimize_Amplitude_Offset_W)

   RR = RietveldRefinery(Rt, 
      factr=1e5,store_intermediate_state=False, show_plots=True)
   # RR.display(RR.minimize_eta)
   RR.display(RR.minimize_Amplitude_Bkgd_Offset_W)
   # # RR.display(RR.minimize_only_Alite)
   # RR.display(RR.minimize_W)
   # RR.display(RR.minimize_All)
   # # RR.display(RR.minimize_All)
   # # RR.display(RR.minimize_All)
   # # RR.display(RR.minimize_All)
   # # RR.display(RR.minimize_All)

   #For fine-tuning
   RR2 = RietveldRefinery(RR.phase_list, factr=1e3,maxiter=10,
      store_intermediate_state=False,show_plots=True)
   RR2.display(RR2.minimize_All)
   RR2.display(RR2.minimize_All)
   # RR2.display(RR2.minimize_All)
   # RR2.display(RR2.minimize_All)


def run():
   exercise_Rietveld_Refinery_Cement()
   print "OK"

if (__name__ == "__main__"):
   run()