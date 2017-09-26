from __future__ import division
import os, random, math
import iotbx.cif, cctbx.miller
import scitbx
from cctbx import xray
from cctbx import crystal
from cctbx.array_family import flex
from cctbx.eltbx import wavelengths
import time
import sys, subprocess
import numpy as np
import matplotlib.pyplot as plt
from scitbx import lbfgsb
import jsonpickle
from libtbx import easy_pickle

import sys, os
sys.path.append(os.path.abspath(".."))

from RietveldPhases import RietveldPhases

class RietveldRefinery:
   """
      This class is used to assemble and organize the inputs required to run a
      series of Rietveld refinements, as per some specifications loaded from
      a file/string.
   """

   def __init__(self,Phase_list):
      """
         Given some list of phases, an instance of the refinery is initialized
         to readily run a refinement process.
      """
      self.Phases = Phase_list
      self.x = []
      self.l = []
      self.u = []
      self.labels = []

      #Copy global variables into x first
      self.x.append(RietveldPhases.two_theta_0)
      self.labels.append("two_theta_0")
      for i,b in np.ndenumerate(RietveldPhases.Bkgd):
         self.x.append(b)
         self.labels.append("Bkgd %d" % i[0])
         # print "Bkgd %d: %f" % (i[0],b)

      for i, Rp in enumerate(Phase_list,start=1):
         self.x.append(Rp.U)
         self.l.append(Rp.U_lower)
         self.u.append(Rp.U_upper)
         self.labels.append("U Phase {0}".format(str(i)))

         # print self.labels[i-1], str(self.x[i-1])