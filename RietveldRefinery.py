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

from RietveldPhases import RietveldPhases
from scipy.optimize import fmin_l_bfgs_b as minimizer

class RietveldRefinery:
   """
      This class is used to assemble and organize the inputs required to run a
      series of Rietveld refinements, as per some specifications loaded from
      a file/string.
   """

   def __init__(self,Phase_list,two_theta,y,input_string, \
      use_bkgd_mask=False,bkgd_delta_theta=0.5):
      """
         Given some list of phases, an instance of the refinery is initialized
         to readily run a refinement process.
      """
      self.x = RietveldPhases.x
      self.Phase_list = Phase_list
      self.two_theta = two_theta
      self.y = y
      self.use_bkgd_mask = use_bkgd_mask
      self.bkgd_delta_theta = bkgd_delta_theta

      self.mask = np.ones(len(self.x),dtype=bool)

      if input_string[-4:] == ".txt":
         self.params_from_file(input_string)
      else:
         self.params_from_string(input_string)

      if self.use_bkgd_mask:
         bkgd_mask = np.invert(np.any( \
            [self.Phase_list[i].bkgd_mask(self.two_theta, \
               bkgd_delta_theta=self.bkgd_delta_theta) \
               for i in xrange(0,len(self.Phase_list))] \
               ,axis=0))
         self.two_theta = self.two_theta[bkgd_mask]
         self.y = self.y[bkgd_mask]

   def params_from_file(self,filename):
      """
         Reads in a set of optimization parameters from a file
         
         :param str filename: the name of the file from which to read parameters

      """
      with open(filename) as file:
         self.params_from_string(file.read())

   def params_from_string(self,input_string):
      """
         Reads in a set of optimization parameters from a string
      
         :param str input_string: a string containing input parameters
      """
      for line in input_string.splitlines():
         if line.split()[0] == "approx_grad":
            self.approx_grad = bool(line.split()[1])
         if line.split()[0] == "factr":
            self.factr = float(line.split()[1])
         if line.split()[0] == "iprint":
            self.iprint = int(line.split()[1])
         if line.split()[0] == "m":
            self.m = int(line.split()[1])
         if line.split()[0] == "pgtol":
            self.pgtol = float(line.split()[1])
         if line.split()[0] == "epsilon":
            self.epsilon = float(line.split()[1])

   def TotalProfile(self):
      # np.set_printoptions(threshold=None)
      # print RietveldPhases.Background_Polynomial(self.two_theta)
      # for i in xrange(0,len(self.Phase_list)):
         # print i 
         # t0 = time.time()
         # print self.Phase_list[i].Phase_Profile(self.two_theta)
         # t1 = time.time()
         # print t1-t0
      # print np.sum(self.Phase_list[i].Phase_Profile(self.two_theta) \
      #    for i in xrange(0,len(self.Phase_list)))
      if self.use_bkgd_mask:
         return RietveldPhases.Background_Polynomial(self.two_theta)
      else:
         return RietveldPhases.Background_Polynomial(self.two_theta) \
            + np.sum( self.Phase_list[i].Phase_Profile(self.two_theta) \
               for i in xrange(0,len(self.Phase_list)))

   def Weighted_Squared_Errors(self):
      # print (self.y - self.TotalProfile())**2/self.y
      # print RietveldPhases.Background_Polynomial(self.two_theta)
      return (self.y - self.TotalProfile())**2/self.y

   def Weighted_Sum_of_Squares(self,x):
      # print self.mask
      self.x['values'][self.mask] = x
      return np.sum(self.Weighted_Squared_Errors())

   def minimize(self):
      t0 = time.time()
      minimizer(self.Weighted_Sum_of_Squares,self.x['values'][self.mask],
         approx_grad = self.approx_grad, \
         factr = self.factr,
         iprint = self.iprint,
         m = self.m,
         pgtol = self.pgtol,
         epsilon = self.epsilon, \
         bounds = zip(self.x['l_limits'][self.mask], \
            self.x['u_limits'][self.mask]))
      t1 = time.time()
      print "Time elapsed: " + str(t1-t0)

   def minimize_Amplitude(self,display_plots = True):
      self.mask = np.char.startswith(self.x['labels'],"Amp")
      self.minimize()

   def minimize_Amplitude_Offset(self,display_plots = True):
      self.mask = np.logical_or( \
         np.char.startswith(self.x['labels'],"Amp"),
         np.char.startswith(self.x['labels'],"two_"))
      self.minimize()

   def minimize_Bkgd(self,display_plots = True):
      self.mask = np.char.startswith(self.x['labels'],"Bkgd")
      self.minimize()

   def minimize_Bkgd_0(self,display_plots = True):
      self.mask = np.char.startswith(self.x['labels'],"Bkgd_0")
      self.minimize()

   def minimize_Amplitude_Bkgd_Offset(self,display_plots = True):
      self.mask = np.logical_or( \
         np.char.startswith(self.x['labels'],"Amp"), \
         np.char.startswith(self.x['labels'],"Bkgd"), \
         np.char.startswith(self.x['labels'],"two_"))
      self.minimize()

   def minimize_All(self,display_plots = True):
      self.mask = np.ones(len(self.x),dtype=bool)
      self.minimize()

   def display(self,fn):
      self.show_multiplot("Before: " + fn.__name__, \
      two_theta_roi=30, \
      delta_theta=10, \
      autohide=False)

      fn(self)

      self.show_multiplot("After: " + fn.__name__, \
      two_theta_roi=30, \
      delta_theta=10, \
      autohide=False)

   def show_plot(self,plottitle,scale_factor=1,autohide=True):
      if autohide:
         plt.ion()
      fig = plt.figure(figsize=(12, 8))
      plt.scatter(self.two_theta,self.y,label='Data',s=1, color='red')
      plt.title(plottitle)

      plt.plot(self.two_theta,scale_factor*self.TotalProfile(), \
         label=r'$I_{\rm calc}$')

      plt.legend(bbox_to_anchor=(.8,.7))
      plt.ylabel(r"$I$")

      plt.show()
      if autohide:
         fig.canvas.flush_events()
         time.sleep(1)
         plt.close('all')

   def show_multiplot(self,plottitle,two_theta_roi=45, \
         delta_theta=0.5,autohide=True):
      if autohide:
         plt.ion()
      plt.figure(figsize=(12, 8))

      plt.subplot(3,1,1)
      plt.scatter(self.two_theta,self.y,label='Data',s=1, color='red')
      plt.title(plottitle)

      plt.plot(self.two_theta,self.TotalProfile(), label=r'$I_{\rm calc}$')
      plt.legend(bbox_to_anchor=(.8,.7))
      plt.ylabel(r"$I$")

      plt.subplot(3,1,2)
      pltmask = np.abs(self.two_theta - two_theta_roi) \
         < delta_theta
      plt.scatter(self.two_theta[pltmask],self.y[pltmask],label='Data',s=1, \
         color='red')

      plt.plot(self.two_theta[pltmask],self.TotalProfile()[pltmask], \
         label=r'$I_{\rm calc}$')
      plt.legend(bbox_to_anchor=(.8,.7))
      plt.ylabel(r"$I$")

      plt.subplot(3,1,3)
      plt.scatter(self.two_theta,self.Weighted_Squared_Errors())
      plt.ylabel(r"$\frac{1}{I} \, (I-I_{\rm calc})^2$")
      plt.xlabel(r'$2\,\theta$')

      plt.show()
      if autohide:
         fig.canvas.flush_events()
         time.sleep(1)
         plt.close('all')