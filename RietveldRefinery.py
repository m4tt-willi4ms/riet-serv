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

   def __init__(self,Phase_list,two_theta,y):
      """
         Given some list of phases, an instance of the refinery is initialized
         to readily run a refinement process.
      """
      self.x = RietveldPhases.x
      self.Phase_list = Phase_list
      self.two_theta = two_theta
      self.y = y

   def TotalProfile(self,delta_theta=0.5):
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
      return RietveldPhases.Background_Polynomial(self.two_theta) \
         + np.sum(self.Phase_list[i].Phase_Profile(self.two_theta)  \
            for i in xrange(0,len(self.Phase_list)))

   def Weighted_Squared_Errors(self):
      return (self.y - self.TotalProfile())**2/self.y

   def Weighted_Sum_of_Squares(self,x):
      self.x['values'] = x
      return np.sum(self.Weighted_Squared_Errors())

   def minimize(self):
      print self.x['labels']
      print minimizer(self.Weighted_Sum_of_Squares,self.x['values'],
         approx_grad = True, \
         factr = 1e10, \
         bounds = zip(self.x['l_limits'],self.x['u_limits']))
      

   def show_plot(self,plottitle,scale_factor=1,autohide=True):
      if autohide:
         plt.ion()
      fig = plt.figure(figsize=(12, 8))
      # plt.subplot(3,1,1)
      plt.scatter(self.two_theta,self.y,label='Data',s=1, color='red')
      plt.title(plottitle)

      plt.plot(self.two_theta,scale_factor*self.TotalProfile(), \
         label=r'$I_{\rm calc}$')

      plt.legend(bbox_to_anchor=(.8,.7))
      plt.ylabel(r"$I$")

      # fig, ax = plt.subplots()
      plt.show()
      if autohide:
         fig.canvas.flush_events()
         time.sleep(1)
         plt.close('all')

   def show_multiplot(self,plottitle,two_theta_roi, \
         delta_theta=0.5,autohide=True):
      # print self.x['values'][self.x['labels'] == 'Amplitude']
      # self.x['values'][self.x['labels'] == 'Amplitude'] \
         # = np.array([0.001,1e-9])
      # print self.x['values'][self.x['labels'] == 'Amplitude']
      if autohide:
         plt.ion()
      plt.figure(figsize=(12, 8))
      plt.subplot(3,1,1)
      # for i in xrange(0,len(two_theta)):
      #    if delta_theta[i]:
      #       color = 'blue'
      #    else: color = 'green'
      #    plt.scatter(two_theta[i],y[i], color=color,s=1)
      plt.scatter(self.two_theta,self.y,label='Data',s=1, color='red')
      plt.title(plottitle)
      # plt.axis([20,60,0,1500])

      plt.plot(self.two_theta,self.TotalProfile(), label=r'$I_{\rm calc}$')
      plt.legend(bbox_to_anchor=(.8,.7))
      plt.ylabel(r"$I$")

      plt.subplot(3,1,2)
      # for i in xrange(0,len(two_theta)):
      #    if delta_theta[i]:
      #       color = 'blue'
      #    else: color = 'green'
      #    plt.scatter(two_theta[i],y[i], color=color,s=1)
      pltmask = np.abs(self.two_theta - two_theta_roi) \
         < delta_theta
      plt.scatter(self.two_theta[pltmask],self.y[pltmask],label='Data',s=1, \
         color='red')
      # plt.title(r"Profile: $Al_2 O_3$")
      # plt.axis([20,60,0,1500])

      plt.plot(self.two_theta[pltmask],self.TotalProfile()[pltmask], \
         label=r'$I_{\rm calc}$')
      plt.legend(bbox_to_anchor=(.8,.7))
      plt.ylabel(r"$I$")

      plt.subplot(3,1,3)
      plt.scatter(self.two_theta,self.Weighted_Squared_Errors())
      plt.ylabel(r"$\frac{1}{I} \, (I-I_{\rm calc})^2$")
      plt.xlabel(r'$2\,\theta$')

      plt.show()#block=False)
      if autohide:
         fig.canvas.flush_events()
         time.sleep(1)
         plt.close('all')