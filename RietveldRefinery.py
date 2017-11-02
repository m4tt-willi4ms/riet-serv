from __future__ import division
import os, random, math
import time
import sys, subprocess
import numpy as np
import itertools
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.optimize import approx_fprime
from scipy.optimize import fmin_l_bfgs_b as minimizer
import json,codecs

from RietveldPhases import RietveldPhases

class RietveldRefinery:
   """
      This class is used to assemble and organize the inputs required to run a
      series of Rietveld refinements, as per some specifications loaded from
      a file/string.
   """

   def __init__(self,Phase_list,input_string, \
      use_bkgd_mask=False,
      bkgd_delta_theta=0.5,
      store_intermediate_state=False,
      show_plots=True):
      """
         Given some list of phases, an instance of the refinery is initialized
         to readily run a refinement process.
      """
      self.x = RietveldPhases.x
      self.Phase_list = Phase_list
      self.two_theta = RietveldPhases.two_theta
      self.y = RietveldPhases.I
      self.use_bkgd_mask = use_bkgd_mask
      self.bkgd_delta_theta = bkgd_delta_theta
      self.store_intermediate_state = store_intermediate_state
      self.show_plots = show_plots

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

      self.del_mask = np.ones(len(self.x['values']),dtype = bool)

      self.sum_y = np.sum(self.y)

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
         if line.split()[0] == "maxiter":
            self.maxiter = float(line.split()[1])
         if line.split()[0] == "epsilon":
            self.epsilon = float(line.split()[1])

   def TotalProfile(self):
      # np.set_printoptions(threshold=None)
      # print RietveldPhases.Background_Polynomial(self.two_theta)
      # for i in xrange(0,len(self.Phase_list)):
      #    print i 
      #    t0 = time.time()
      #    print self.Phase_list[i].Phase_Profile(self.two_theta)
      #    t1 = time.time()
      #    print "Time elapsed, phase " + str(i) +": " + str(t1-t0)
      # print np.sum(self.Phase_list[i].Phase_Profile(self.two_theta) \
      #    for i in xrange(0,len(self.Phase_list)))
      if self.use_bkgd_mask:
         if self.store_intermediate_state:
            self.TotalProfile_state = \
               RietveldPhases.Background_Polynomial(self.two_theta)
         else:
            return RietveldPhases.Background_Polynomial(self.two_theta)
      else:
         # t0 = time.time()
         # print RietveldPhases.Background_Polynomial(self.two_theta) \
         #    + np.sum( self.Phase_list[i].Phase_Profile(self.two_theta) \
         #       for i in xrange(0,len(self.Phase_list)))
         # t1 = time.time()
         # print "TotalProfile Time: " + str(t1-t0)
         # print map(lambda x: np.amax(self.Phase_list[x].Phase_Profile()), 
         #       range(0,len(self.Phase_list)))
         # print np.sum( x \
         #       for x in self.Profile_generator())
         if self.store_intermediate_state:
            self.TotalProfile_state = \
               RietveldPhases.Background_Polynomial(self.two_theta) \
               + np.sum( x for x in self.Phase_Profiles())
            # + np.sum(map(lambda x: self.Phase_list[x].Phase_Profile(), 
            #    range(0,len(self.Phase_list))))
         else: 
            return RietveldPhases.Background_Polynomial(self.two_theta) \
               + np.sum( x for x in self.Phase_Profiles())
      return self.TotalProfile_state

   def Phase_Profiles(self):
      for i in xrange(0,len(self.Phase_list)):
         yield self.Phase_list[i].Phase_Profile()

   def Weighted_Squared_Errors(self):
      self.TotalProfile_state = self.TotalProfile()
      if self.store_intermediate_state:
         self.Weighted_Squared_Errors_state = \
         (self.y - self.TotalProfile_state)**2/self.y
         try:
            json.dump(self.TotalProfile().tolist(), 
               codecs.open("current_profile.json", 'w', encoding='utf-8'), 
               separators=(',', ':'), 
               # sort_keys=True, 
               indent=4)
            json.dump(self.x.tolist(), 
               codecs.open("xparams.json", 'w', encoding='utf-8'), 
               separators=(',', ':'), 
               # sort_keys=True, 
               indent=4)
         except IOError:
            pass
      else:
         return (self.y - self.TotalProfile_state)**2/self.y
      return self.Weighted_Squared_Errors_state

   def Relative_Differences(self):
      return (self.TotalProfile_state-self.y)/self.y

   def Update_state(self):
      self.TotalProfile_state = self.TotalProfile()
      self.Relative_Differences_state = self.Relative_Differences()
      # self.Relative_Differences_state.shape = \
      #    (self.Relative_Differences_state.shape[0],1)
      self.Weighted_Squared_Errors_state = self.Weighted_Squared_Errors()

   def Weighted_Sum_of_Squares(self,x):
      # print self.mask
      self.x['values'][self.mask] = x
      self.Update_state()
      return np.sum(self.Weighted_Squared_Errors())

   def Weighted_Sum_of_Squares_Grad(self,x):
      self.x['values'][self.mask] = x

      result = np.zeros(len(self.x[self.mask]))

      bkgd_mask = np.logical_and(self.mask, 
         np.char.startswith(self.x['labels'],"Bkgd"))[self.mask]
      if np.any(bkgd_mask):
         result[bkgd_mask] = \
            np.sum(2*np.multiply(RietveldPhases.two_theta_powers,
               self.Relative_Differences_state),axis=1)

      for i in xrange(0,len(self.Phase_list)):
         mask = np.logical_and(self.mask,np.logical_or(
            RietveldPhases.get_global_parameters_mask(include_bkgd=False),
            self.Phase_list[i].get_phase_parameters_mask())
            )
         if np.any(mask):
            result[mask[self.mask]] += np.sum(2*np.multiply(self.Phase_list[i].
               Phase_Profile_Grad(mask, epsilon=self.epsilon),
               self.Relative_Differences_state),axis=1)
      return result

   def minimize(self):
      self.t0 = time.time()
      self.mask = np.logical_and(self.del_mask,self.mask)
      self.display_parameters()
      minimizer(self.Weighted_Sum_of_Squares,self.x['values'][self.mask],
         callback = self.callback, \
         fprime = self.Weighted_Sum_of_Squares_Grad, \
         factr = self.factr,
         iprint = self.iprint,
         m = self.m,
         maxiter = self.maxiter,
         pgtol = self.pgtol,
         epsilon = self.epsilon, \
         bounds = zip(self.x['l_limits'][self.mask], \
            self.x['u_limits'][self.mask]))
      self.t1 = time.time()
      print "\n"

   def Update_Phase_list(self):
      Amp_mask = np.char.startswith(self.x['labels'],"Amp")
      Amp_vals = self.x['values'][Amp_mask][self.del_mask[Amp_mask]]
      Amp_total = np.sum(Amp_vals)
      del_list = []
      for i,Amp_val in np.ndenumerate(Amp_vals):
         if Amp_val/Amp_total < 1e-5:
            del_list += [False]
            self.del_mask = np.logical_and(self.del_mask, \
               np.invert(np.isin(np.array(range(0,len(RietveldPhases.x))), \
               np.array(range(self.Phase_list[i[0]].phase_0_index, \
                  self.Phase_list[i[0]].phase_0_index \
                     +self.Phase_list[i[0]].num_params)))))
         else:
            del_list += [True]
      if not all(del_list):
         print "\nLength of Phase_list (Before): " + str(len(self.Phase_list))
         self.Phase_list = [self.Phase_list[i] for i in \
            xrange(len(self.Phase_list)) if del_list[i]]
         print "Length of Phase_list (After): " + str(len(self.Phase_list)) \
            + "\n"

   def callback(self,x):

      sys.stdout.write('.')
      sys.stdout.flush()
      if self.show_plots:
         self.update_plot()

   def minimize_Amplitude(self,display_plots = True):
      self.mask = np.char.startswith(self.x['labels'],"Amp")
      self.minimize()

   def minimize_Amplitude_Offset(self,display_plots = True):
      self.mask = np.logical_or( \
         np.char.startswith(self.x['labels'],"Amp"),
         np.char.startswith(self.x['labels'],"two_"))
      self.minimize()
      self.Update_Phase_list()

   def minimize_only_Alite(self,display_plots = True):
      self.mask = np.logical_or( \
                     np.logical_or( \
                        np.char.startswith(self.x['labels'],"Amp"), \
                        np.char.startswith(self.x['labels'],"two_")),
                     np.isin(np.array(range(0,len(RietveldPhases.x))), \
               np.array(range(self.Phase_list[0].phase_0_index, \
                  self.Phase_list[0].phase_0_index \
                     +self.Phase_list[0].num_params))))
      self.minimize()

   def minimize_First_n_Phases(self,n=1,display_plots = True):
      tmp = np.argsort(self.x['values'] 
         [np.char.startswith(self.x['labels'],"Amp")])
      n_phases = [self.Phase_list[i] for i in tmp]

      tmp = np.zeros(len(self.x),dtype=bool)
      for i in xrange(0,n):
         tmp = np.logical_or(tmp, \
            np.isin(np.array(range(0,len(RietveldPhases.x))), \
               np.array(range(n_phases[i].phase_0_index, \
                  n_phases[i].phase_0_index \
                     +n_phases[i].num_params))))
      self.mask = np.logical_or( \
                     np.logical_or( \
                        np.char.startswith(self.x['labels'],"Amp"), \
                        np.char.startswith(self.x['labels'],"two_")),
                     tmp)
      self.minimize()

   def minimize_Amplitude_Offset_W(self,display_plots = True):
      self.mask = np.logical_or( \
         np.char.startswith(self.x['labels'],"Amp"), \
            np.logical_or(np.char.startswith(self.x['labels'],"two_"), \
               self.x['labels'] == "W"))
      self.minimize()
      self.Update_Phase_list()

   def minimize_Amplitude_Bkgd_Offset_W(self,display_plots = True):
      self.mask = np.logical_or( 
         np.char.startswith(self.x['labels'],"Amp"), 
            np.logical_or(np.char.startswith(self.x['labels'],"two_"), 
               np.logical_or(np.char.startswith(self.x['labels'],"Bkgd"), 
                  self.x['labels'] == "W")))
      self.minimize()
      self.Update_Phase_list()

   def minimize_Bkgd(self,display_plots = True):
      self.mask = np.char.startswith(self.x['labels'],"Bkgd")
      self.minimize()

   def minimize_Bkgd_0(self,display_plots = True):
      self.mask = np.char.startswith(self.x['labels'],"Bkgd_0")
      self.minimize()

   def minimize_Amplitude_Bkgd_Offset(self,display_plots = True):
      self.mask = np.logical_or( \
         np.char.startswith(self.x['labels'],"Amp"), \
            np.logical_or(np.char.startswith(self.x['labels'],"Bkgd"), \
               np.char.startswith(self.x['labels'],"two_")))
      self.minimize()

   def minimize_Amplitude_Bkgd(self,display_plots = True):
      self.mask = np.logical_or( \
         np.char.startswith(self.x['labels'],"Amp"), \
         np.char.startswith(self.x['labels'],"Bkgd"))
      self.minimize()

   def minimize_All(self,display_plots = True):
      self.mask = np.ones(len(self.x),dtype=bool)
      self.minimize()

   def minimize_All_except_Amplitude(self,display_plots = True):
      self.mask = np.invert(np.logical_or(np.zeros(len(self.x),dtype=bool), \
         np.char.startswith(self.x['labels'],"Amp")))
      self.minimize()

   def display(self,fn,**kwargs):
      # self.show_multiplot("Before: " + fn.__name__, \
      # two_theta_roi=30, \
      # delta_theta=10, \
      # autohide=False)
      if self.show_plots:
         self.show_multiplot("Progress: " + fn.__name__, \
         two_theta_roi=32.5, \
         delta_theta=3, \
         autohide=False)

      if kwargs is not None:
         fn(**kwargs)
      else:
         fn()

      if (self.store_intermediate_state and self.show_plots):
         self.update_plot()
         self.fig.suptitle(fn.__name__)


      if ((not self.store_intermediate_state) and self.show_plots):
         self.show_multiplot("After: " + fn.__name__, \
         two_theta_roi=32.5, \
         delta_theta=3, \
         autohide=False)
      self.display_parameters(fn)
      self.display_stats(fn)
      if self.show_plots:
         plt.ioff()
         plt.show()

   def update_plot(self):
      self.current_profile.set_ydata(self.TotalProfile_state)
      self.current_profile_masked.set_ydata(self.TotalProfile_state 
         [self.pltmask])
      self.residuals.set_ydata(self.Weighted_Squared_Errors_state)
      self.fig.canvas.draw()

   def display_parameters(self,fn=None):
      if fn is not None:
         print "After " + fn.__name__ + ":"
      for label,value,l,u in zip( \
         RietveldPhases.x['labels'][self.mask], \
         RietveldPhases.x['values'][self.mask], \
         RietveldPhases.x['l_limits'][self.mask], \
         RietveldPhases.x['u_limits'][self.mask]):
         print label + " = " + str(value) + " (" + str(l) + ", " + str(u) + ")"

   def display_stats(self,fn=None):
      self.mask = np.ones(len(self.x),dtype=bool)
      WSS = self.Weighted_Sum_of_Squares(self.x['values'])
      R_wp = np.sqrt(WSS/self.sum_y)
      R_e = np.sqrt((len(self.two_theta)-len(self.x))/self.sum_y)
      if fn is not None:
         print "\nTime taken to run " + fn.__name__ + ": " \
            + str(round(self.t1-self.t0,3)) + " seconds"
      print "R_wp: " + str(R_wp)
      print "R_e: " + str(R_e)
      print "Goodness-of-Fit: " + str(R_wp/R_e)

      if not self.use_bkgd_mask:
         Amplitudes = self.x['values'] \
            [np.char.startswith(self.x['labels'],"Amp")]
         total = np.sum(Amplitudes) 
         print "\n"
         for i,val in np.ndenumerate(Amplitudes):
            print "Phase " + str(i[0]+1) + ": " + str(val/total*100) + " %"
         print "\n"

   def show_plot(self,plottitle,scale_factor=1,autohide=True):
      if autohide:
         plt.ion()
      # fig = plt.figure(figsize=(9,7))
      fig = plt.figure(figsize=(10,8))
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
         delta_theta=0.5,autohide=True,interactive=True):
      if autohide or interactive:
         plt.ion()

      # self.fig = plt.figure(figsize=(7,5))
      self.fig = plt.figure(figsize=(6,4))
      self.fig.suptitle(plottitle)

      self.subplot1 = self.fig.add_subplot(311) #plt.subplot(3,1,1)
      plt.scatter(self.two_theta,self.y,label='Data',s=1, color='red')

      self.current_profile, = self.subplot1.plot(self.two_theta,
         self.TotalProfile(), label=r'$I_{\rm calc}$')
      plt.legend(bbox_to_anchor=(.8,.7))
      plt.ylabel(r"$I$")

      subplot2 = self.fig.add_subplot(312) #plt.subplot(3,1,2)
      self.pltmask = np.abs(self.two_theta - two_theta_roi) \
         < delta_theta
      plt.scatter(self.two_theta[self.pltmask],self.y[self.pltmask],
         label='Data',s=2, color='red')

      self.current_profile_masked, = subplot2.plot(self.two_theta[self.pltmask],
         self.TotalProfile()[self.pltmask], label=r'$I_{\rm calc}$')
      # plt.legend(bbox_to_anchor=(.8,.7))
      plt.ylabel(r"$I$")

      subplot3 = self.fig.add_subplot(313) #plt.subplot(3,1,3)
      self.residuals, = subplot3.plot(self.two_theta,
         self.Weighted_Squared_Errors(),'bo',ms=2)
      plt.ylabel(r"$\frac{1}{I} \, (I-I_{\rm calc})^2$")
      plt.xlabel(r'$2\,\theta$')

      self.fig.canvas.draw()

      if autohide:
         fig.canvas.flush_events()
         time.sleep(1)
         plt.close('all')
