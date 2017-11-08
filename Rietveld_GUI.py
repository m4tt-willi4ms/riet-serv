from __future__ import division
# from scitbx.array_family import flex
# import scitbx.math
# from libtbx.test_utils import approx_equal, eps_eq, Exception_expected
# import sys
import numpy as np
import time
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import \
   FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib import style

from random import randrange

import sys, os
# sys.path.append(os.path.abspath("../.."))

from src.RietveldPhases import RietveldPhases
from src.RietveldRefinery import RietveldRefinery

from scitbx import lbfgsb
from cctbx.eltbx import wavelengths
from libtbx import test_utils
import libtbx.load_env

import Tkinter as tk

LARGE_FONT = ("Verdana", 12)
# style.use("ggplot")

class RietveldGUI(tk.Tk):
   def __init__(self, *args, **kwargs):
      tk.Tk.__init__(self, *args, **kwargs)

      tk.Tk.iconbitmap(self, default=u"doc//ProtoLogo_R.ico")
      tk.Tk.wm_title(self, "Rietveld Refinement")
      # tk.Tk.wm_geometry(self,'1000x800')

      container = tk.Frame(self)

      container.grid_rowconfigure(0, weight=1)
      container.grid_columnconfigure(0, weight=1)

      self.frames = {}

      # frame = Example(container,self)
      frame = RefinementParameterControl(container,self)
      self.frames[RefinementParameterControl] = frame
      frame.grid(row=0, column=0, sticky="nsew")

      self.show_frame(RefinementParameterControl)

   def show_frame(self, cont):
      frame = self.frames[cont]
      frame.tkraise()

class RefinementParameterControl(tk.Frame):
   def __init__(self, parent, controller, *args, **kwargs):
      tk.Frame.__init__(self, parent)
      self.checkbutton = tk.Checkbutton(self, command=self.checkbutton_clicked,
         *args, **kwargs)
      self.checkbutton.pack(side="left")

   def checkbutton_clicked(self):
      print "here"

class Example(tk.Frame):

   def __init__(self, parent,controller):
      tk.Frame.__init__(self, parent)
      self.RPControl = RefinementParameterControl(self, controller,
         text="Param #1", font = LARGE_FONT)
      self.RPControl.pack(side="top", fill="both", expand=True)

if __name__ == "__main__":
   root = RietveldGUI()
   root.mainloop()

input_strings = ["""\
U              0.00    0     0.1
V              -0.00   -0.1   0
W              0.01   0.0001     1
Amplitude         0.1 0      inf
eta:           3
unit_cell_a    0.01
unit_cell_b    0.01
unit_cell_c    0.01
unit_cell_alpha   0.005
unit_cell_beta    0.005
unit_cell_gamma   0.005
""",
"""\
U              0.00    0     0.1
V              -0.00   -0.1   0
W              0.01   0.0001     1
Amplitude         0.000001 0      inf
eta:           2
""",
"""\
U              0.00    0     0.1
V              -0.00   -0.1   0
W              0.01   0.0001     1
Amplitude         0.1 0      inf
eta:           2
""",
"""\
U              0.00    0     0.1
V              -0.00   -0.1   0
W              0.01   0.0001     1
Amplitude         0.1 0      inf
eta:           2
""",
"""\
U              0.00    0     0.1
V              -0.00   -0.1   0
W              0.01   0.0001     1
Amplitude         0.1 0      inf
eta:           2
""",
"""\
U              0.00    0     0.1
V              -0.00   -0.1   0
W              0.01   0     1
Amplitude         0.1 0      inf
eta:           2
""",
"""\
U              0.00    0     0.1
V              -0.00   -0.1   0
W              0.01   0.0001     1
Amplitude         0.1 0      inf
eta:           2
""",
"""\
U              0.00    0     0.1
V              -0.00   -0.1   0
W              0.01   0.0001     1
Amplitude         0.1 0      inf
eta:           2
"""]

global_input_string = """\
Bkgd:          3
two_theta_0       0.      -0.5  0.5
"""

bkgd_minimizer_input_string = """\
factr       1e10
iprint      -1
maxiter     150
m           10
pgtol       1e-5
epsilon     1e-13
"""

minimizer_input_string = """\
factr       1e2
iprint      -1
maxiter     150
m           10
pgtol       1e-5
epsilon     1e-13
"""

fine_minimizer_input_string = """\
factr       1e2
iprint      1
maxiter     150
m           15
pgtol       1e-5
epsilon     1e-13
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
with open(r"data//profiles//cement_15_03_11_0028.xye") as file:
   for line in file.readlines()[1:]:
      two_thetatmp, ytmp, ztmp = line.split()
      # two_thetatmp, ytmp = line.split()
      # if float(two_thetatmp) < 15:
      tst_two_theta.append(float(two_thetatmp))
      tst_y.append(float(ztmp)**2)
tst_two_theta = np.array(tst_two_theta)
# mask = np.ones(len(tst_two_theta),dtype=bool)
mask = tst_two_theta > 20 
# mask = np.logical_and(tst_two_theta >25,np.logical_or(tst_two_theta<33.75,
#    tst_two_theta>34.3))
# mask = np.logical_or(tst_two_theta<33.75,tst_two_theta>34.3)
tst_two_theta = tst_two_theta[mask]
tst_y = np.array(tst_y)[mask]

def exercise_Rietveld_Refinery_Cement():
   # RietveldPhase.fromstring(input_string) 
   cifs = [
      "1540705-Alite.cif", 
      "9012789-Belite.cif", 
      "1200009-Ferrite.cif", 
      "1000039-AluminateCubic.cif", 
      "9014308-AluminateOrtho.cif", 
      "9007569-Arcanite.cif",
      "1011094-FreeLime.cif", 
      "1000053-Periclase.cif", 
      ]
   Rt = []

   print "cifs: \n" 
   for p in cifs:
      print p
   print "\nInput String: \n"
   for i,p in enumerate(input_strings):
      print "Phase " + str(i+1) + ": \n" + p
   print "Global Input String: \n" + global_input_string
   print "Minimizer Input String: \n" + minimizer_input_string
   print "Fine Minimizer Input String: \n" + fine_minimizer_input_string


   CU_wavelength = wavelengths.characteristic("CU").as_angstrom()
   d_min = CU_wavelength/2/np.sin(np.pi/360*tst_two_theta[-1])
   d_max = CU_wavelength/2/np.sin(np.pi/360*tst_two_theta[0])
   tst_y_max = np.amax(tst_y)/len(cifs)

   RietveldPhases.global_params_from_string(global_input_string,
      tst_two_theta,tst_y)

   for cif, input_string in zip(cifs,input_strings):
   #    tt0 = time.time()
      Rt.append(RietveldPhases(cif,input_string,d_min,d_max, \
         I_max = tst_y_max, delta_theta=1.5,Intensity_Cutoff = 0.005))
   #    tt1 = time.time()
   #    print "TIME TO READ IN: " +str(tt1-tt0) + " seconds"

   # for i,Rp in enumerate(Rt):
   #    tmp = Rp.two_theta_peaks[np.abs(Rp.two_theta_peaks-34) <0.5]
   #    tmp2 = Rp.weighted_intensities[np.abs(Rp.two_theta_peaks-34) <0.5]
   #    print str(i) + ": " + str(tmp)
   #    print str(i) + ": " + str(tmp2)

   # numpeaks = 0
   # for i,Rp in enumerate(Rt):
   #    print Rp.two_theta_peaks.shape[0]
   #    numpeaks += Rp.two_theta_peaks.shape[0]
   # print numpeaks

   # First fit the background
   # RR = RietveldRefinery(Rt,bkgd_minimizer_input_string, \
   #    use_bkgd_mask=False,bkgd_delta_theta=0.05,
   #    store_intermediate_state=True, show_plots=True)
   # RR.display(RR.minimize_Bkgd)

   # #Now use the full dataset
   # RR = RietveldRefinery(Rt,minimizer_input_string,
   #    store_intermediate_state=True, show_plots=True)

   # # RR.display(RR.minimize_Bkgd)
   # # RR.display(RR.minimize_Bkgd_Offset)
   # # RR.display(RR.minimize_Amplitude)
   # # RR.display(RR.minimize_Amplitude)
   # RR.display(RR.minimize_Amplitude_Offset)
   # # RR.display(RR.minimize_Amplitude_Offset_unit_cell)
   # RR.display(RR.minimize_unit_cell)
   # # RR.display(RR.minimize_First_n_Phases)
   # # RR.display(RR.minimize_First_n_Phases,n=3)
   # # RR.display(RR.minimize_Amplitude_Offset_W)
   # RR.display(RR.minimize_Amplitude_Bkgd_Offset_W)
   # # RR.display(RR.minimize_Amplitude_Bkgd_Offset)
   # # RR.display(RR.minimize_only_Alite)
   # # RR.display(RR.minimize_All)
   # # RR.display(RR.minimize_All)
   # # RR.display(RR.minimize_All)
   # # RR.display(RR.minimize_All)
   # # RR.display(RR.minimize_All)

   # #For fine-tuning
   # RR2 = RietveldRefinery(RR.Phase_list,
   #    fine_minimizer_input_string,store_intermediate_state=True,show_plots=True)
   # RR2.display(RR2.minimize_All)
   # # RR2.display(RR2.minimize_All)
   # # RR2.display(RR2.minimize_All)
   # # RR2.display(RR2.minimize_All)


# def run():
#    exercise_Rietveld_Refinery_Cement()
#    print "OK"

# if (__name__ == "__main__"):
#    run()