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
from scipy.optimize import minimize

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
      self.x = RietveldPhases.x_global
      self.d = len(Phase_list)
      self.x['l_limits'][1] = -4.0
      print tuple(self.x)
      print len(tuple(self.x))
      print len(tuple(tuple(Phase_list[i].x) for i in xrange(0,len(Phase_list))))
      print tuple(Phase_list[i].x for i in xrange(0,len(Phase_list)))
      print np.concatenate(tuple(self.x)+tuple(Phase_list[i].x for i in \
         xrange(0,self.d)))

      self.x = np.concatenate((self.x,Phase_list[0].x))
      self.x['l_limits'][5] = -6.0
      print self.x['values']
      # print len(Phase_list)
      # # print self.x.all()
      # for i in xrange(0,len(Phase_list)):
      #    self.x = np.append(self.x,Phase_list[i].x)
      # print self.x
      # # print np.append(RietveldPhases.x_global,Phase_list[0].x)
      # self.x['l_limits'][0] = -3.0
      # print self.x

      # print tuple(Phase_list)
      # tmp = ()

      # print RietveldPhases.x_global
      # print np.concatenate((RietveldPhases.x_global,Phase_list[0].x))
      # print Phase_list[0].x
      # tmp += (Phase_list[0].x)
      # print tmp  
      # print (RietveldPhases.x_global)+tuple(Phase_list)
      # x = np.concatenate((RietveldPhases.x_global)+tuple(Phase_list))
      # print x
      # #Copy global variables into x first
      # self.x.append(RietveldPhases.two_theta_0)
      # self.labels.append("two_theta_0")
      # for i,b in np.ndenumerate(RietveldPhases.Bkgd):
      #    self.x.append(b)
      #    self.labels.append("Bkgd %d" % i[0])
      #    # print "Bkgd %d: %f" % (i[0],b)

      # for i, Rp in enumerate(Phase_list,start=1):
      #    self.x.append(Rp.U)
      #    self.l.append(Rp.U_lower)
      #    self.u.append(Rp.U_upper)
      #    self.labels.append("U Phase {0}".format(str(i)))

      #    # print self.labels[i-1], str(self.x[i-1])