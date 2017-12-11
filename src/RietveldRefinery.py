from __future__ import division
import os, random, math
import time
import sys, subprocess
import numpy as np
import itertools
import matplotlib.pyplot as plt
from scipy.optimize import approx_fprime
from scipy.optimize import fmin_l_bfgs_b as minimizer
import json,codecs

from RietveldPhases import RietveldPhases

default_factr = 1e7
default_iprint = 1
default_maxiter = 150
default_m = 15
default_pgtol = 1e-5
default_epsilon = 1e-8

class RietveldRefinery:
   """
      This class is used to assemble and organize the inputs required to run a
      series of Rietveld refinements, as per some specifications loaded from
      a file/string.
   """


   def __init__(self,phase_list, \
      bkgd_refine=False,
      store_intermediate_state=False,
      show_plots=False,
      factr=default_factr,
      iprint=default_iprint,
      maxiter=default_maxiter,
      m=default_m,
      pgtol=default_pgtol,
      epsilon=default_epsilon,
      ):
      """
         Given some list of phases, an instance of the refinery is initialized
         to readily run a refinement process.
      """
      self.phase_list = phase_list
      self.two_theta = RietveldPhases.two_theta
      self.y = RietveldPhases.I
      self.bkgd_refine = bkgd_refine
      self.store_intermediate_state = store_intermediate_state
      self.show_plots = show_plots
      self.factr = factr
      self.iprint = iprint
      self.maxiter = maxiter
      self.m = 15
      self.pgtol = pgtol
      self.epsilon = epsilon


      # if input_string[-4:] == ".txt":
      #    self.params_from_file(input_string)
      # else:
      #    self.params_from_string(input_string)

      # if self.bkgd_refine:
      #    bkgd_mask = np.invert(np.any( \
      #       [self.Phase_list[i].bkgd_mask(self.two_theta, \
      #          bkgd_delta_theta=self.bkgd_delta_theta) \
      #          for i in xrange(0,len(self.Phase_list))] \
      #          ,axis=0))
      #    self.two_theta = self.two_theta[bkgd_mask]
      #    self.y = self.y[bkgd_mask]

      # self.del_mask = np.ones(len(self.x['values']),dtype = bool)

      self.sum_y = np.sum(self.y)

      if self.bkgd_refine:
         self.raw_profile_state = np.sum( x for x in self.phase_profiles())
      
      RietveldPhases.assemble_global_x()
      for phase in self.phase_list:
         phase.assemble_phase_x()
      self.x = np.hstack((RietveldPhases.global_x,
         np.hstack((x.phase_x for x in self.phase_list))))

      self.mask = np.zeros(len(self.x),dtype=bool)

      self.global_mask = np.isin(np.array(xrange(len(self.x))),
         np.array(xrange(len(RietveldPhases.global_x))))
      self.bkgd_mask = np.char.startswith(self.x['labels'],"bkgd")
      self.global_mask_no_bkgd = np.logical_xor(self.global_mask,self.bkgd_mask)

      n=len(RietveldPhases.global_x)
      self.phase_masks = []
      self.global_and_phase_masks = []
      for phase in self.phase_list:
         self.phase_masks.append(np.isin(np.array(xrange(len(self.x))),
            np.array(xrange(n,n+len(phase.phase_x)))))
         self.global_and_phase_masks.append(np.logical_or(
            self.global_mask_no_bkgd,self.phase_masks[-1]))
         n+=len(phase.phase_x)

      self.update_state()

   def total_profile(self):
      if self.bkgd_refine:
         return RietveldPhases.background_polynomial() \
            + self.raw_profile_state
      else:
         return RietveldPhases.background_polynomial() \
            + np.sum( x for x in self.phase_profiles())

   def update_xs(self):
      if np.any(self.mask[self.global_mask_no_bkgd]):
         RietveldPhases.update_global_x_no_bkgd(
            self.x[self.global_mask_no_bkgd])
      if np.any(self.mask[self.bkgd_mask]):
         RietveldPhases.update_bkgd(self.x[self.bkgd_mask])
      for i in xrange(len(phase_list)):
         if np.any(self.mask[self.phase_masks[i]]):
            self.phase_list[i].update_phase_x(self.x[self.phase_masks[i]])

   def total_profile_x(self,x):
      self.x['values'][self.mask] = x
      self.update_state()
      if self.bkgd_refine:
         return RietveldPhases.background_polynomial() \
            + self.raw_profile_state
      else:
         return RietveldPhases.background_polynomial() \
            + np.sum( x for x in self.phase_profiles())

   def phase_profiles(self):
      for i in xrange(0,len(self.phase_list)):
         yield self.phase_list[i].phase_profile()

   def weighted_squared_errors(self):
         # try:
         #    json.dump(self.TotalProfile().tolist(), 
         #       codecs.open("current_profile.json", 'w', encoding='utf-8'), 
         #       separators=(',', ':'), 
         #       # sort_keys=True, 
         #       indent=4)
         #    json.dump(self.x.tolist(), 
         #       codecs.open("xparams.json", 'w', encoding='utf-8'), 
         #       separators=(',', ':'), 
         #       # sort_keys=True, 
         #       indent=4)
         # except IOError:
         #    pass
      return (self.y - self.total_profile_state)**2/self.y

   def relative_differences(self):
      return (self.total_profile_state-self.y)/self.y

   def update_state(self):
      # if not self.bkgd_refine:
      RietveldPhases.update_global_x(self.x[self.global_mask])
      if not self.bkgd_refine:
         for i in xrange(len(self.phase_list)):
            self.phase_list[i].phase_x = self.x[self.phase_masks[i]]
            self.phase_list[i].update_params()
               # self.x[self.phase_masks[i]][self.mask[self.phase_masks[i]]])
      self.total_profile_state = self.total_profile()
      self.relative_differences_state = self.relative_differences()
      # self.Relative_Differences_state.shape = \
      #    (self.Relative_Differences_state.shape[0],1)
      self.weighted_squared_errors_state = self.weighted_squared_errors()

   def weighted_sum_of_squares(self,x):
      # print self.mask
      self.x['values'][self.mask] = x
      # if not self.bkgd_refine:
      self.update_state()
      return np.sum(self.weighted_squared_errors_state)

   def weighted_sum_of_squares_grad(self,x):
      self.x['values'][self.mask] = x

      result = np.zeros(len(self.x[self.mask]))

      bkgd_mask = np.logical_and(self.mask, self.bkgd_mask)[self.mask]
      if np.any(bkgd_mask):
         result[bkgd_mask] = \
            np.sum(2*np.multiply(
               RietveldPhases.two_theta_powers[:len(RietveldPhases.bkgd)],
               self.relative_differences_state),axis=1)
      # print "here"
      # global_mask = np.logical_and(self.mask,
      #    self.global_mask_no_bkgd)[self.mask]
      # if np.any(global_mask):
      #    result[global_mask] = np.sum(2*np.multiply(self.phase_list[i].
      #          phase_profile_grad(mask, epsilon=self.epsilon),
      #          self.relative_differences_state),axis=1)

      for i in xrange(0,len(self.phase_list)):
         mask = np.logical_and(self.mask,self.global_and_phase_masks[i]) 
         if np.any(mask[self.global_and_phase_masks[i]]):
            result[mask[self.mask]] += \
               np.sum(2*np.multiply(self.phase_list[i].
               phase_profile_grad(mask[self.global_and_phase_masks[i]], 
                  epsilon=self.epsilon),
               self.relative_differences_state),axis=1)
      return result

   def minimize(self):
      self.t0 = time.time()
      # self.mask = np.logical_and(self.del_mask,self.mask)
      self.display_parameters()
      # if self.bkgd_refine:
      #    self.raw_profile_state = self.total_profile()
      minimizer(self.weighted_sum_of_squares,self.x['values'][self.mask],
         callback = self.callback, \
         fprime = self.weighted_sum_of_squares_grad, \
         # approx_grad = True,
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

   def update_phase_list(self):
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
      elif self.store_intermediate_state:
         self.update_state()

   def minimize_with_mask(self,mask):
      self.mask = mask
      self.minimize()

   def minimize_Amplitude(self,display_plots = True):
      self.mask = np.char.startswith(self.x['labels'],"Amp")
      self.minimize()

   def minimize_Amplitude_Offset(self,display_plots = True):
      self.mask = np.logical_or( \
         np.char.startswith(self.x['labels'],"Amp"),
         np.char.startswith(self.x['labels'],"two_"))
      self.minimize()
      # self.Update_Phase_list()

   def minimize_Amplitude_Offset_unit_cell(self,display_plots = True):
      self.mask = np.logical_or(np.logical_or( \
         np.char.startswith(self.x['labels'],"Amp"),
         np.char.startswith(self.x['labels'],"two_")),
         np.char.startswith(self.x['labels'],"unit_cell"))
      self.minimize()
      self.Update_Phase_list()

   def minimize_unit_cell(self,display_plots = True):
      self.mask = np.char.startswith(self.x['labels'],"unit_cell")
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
      # self.Update_Phase_list()

   def minimize_Amplitude_Bkgd_Offset_W(self,display_plots = True):
      self.mask = np.logical_or( 
         np.char.startswith(self.x['labels'],"Amp"), 
            np.logical_or(np.char.startswith(self.x['labels'],"two_"), 
               np.logical_or(np.char.startswith(self.x['labels'],"bkgd"), 
                  self.x['labels'] == "W")))
      self.minimize()
      # self.Update_Phase_list()

   def minimize_bkgd(self,display_plots = True):
      self.mask = self.bkgd_mask
      self.minimize()

   def minimize_Bkgd_Offset(self,display_plots = True):
      self.mask = np.logical_or(
         np.char.startswith(self.x['labels'],"bkgd"),
         np.char.startswith(self.x['labels'],"two_"))
      self.minimize()

   def minimize_Bkgd_0(self,display_plots = True):
      self.mask = np.char.startswith(self.x['labels'],"Bkgd_0")
      self.minimize()

   def minimize_Amplitude_Bkgd_Offset(self,display_plots = True):
      self.mask = np.logical_or( \
         np.char.startswith(self.x['labels'],"Amp"), \
            np.logical_or(np.char.startswith(self.x['labels'],"bkgd"), \
               np.char.startswith(self.x['labels'],"two_")))
      self.minimize()

   def minimize_Amplitude_Bkgd(self,display_plots = True):
      self.mask = np.logical_or( \
         np.char.startswith(self.x['labels'],"Amp"), \
         np.char.startswith(self.x['labels'],"bkgd"))
      self.minimize()

   def minimize_eta(self,display_plots = True):
      self.mask = np.char.startswith(self.x['labels'],"eta")
      self.minimize()

   def minimize_W(self,display_plots = True):
      self.mask = np.char.startswith(self.x['labels'],"W")
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
         # self.total_profile_state = self.total_profile()
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
         # self.total_profile_state = self.total_profile()
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
      self.current_profile.set_ydata(self.total_profile_state)
      self.current_profile_masked.set_ydata(self.total_profile_state 
         [self.pltmask])
      self.residuals.set_ydata(self.weighted_squared_errors_state)
      self.fig.canvas.draw()

   def display_parameters(self,fn=None):
      if fn is not None:
         print "After " + fn.__name__ + ":"
      for label,value,l,u in zip( \
         self.x['labels'][self.mask], \
         self.x['values'][self.mask], \
         self.x['l_limits'][self.mask], \
         self.x['u_limits'][self.mask]):
         print label + " = " + str(value) + " (" + str(l) + ", " + str(u) + ")"

   def display_stats(self,fn=None):
      self.mask = np.ones(len(self.x),dtype=bool)
      WSS = self.weighted_sum_of_squares(self.x['values'])
      R_wp = np.sqrt(WSS/self.sum_y)
      R_e = np.sqrt((len(self.two_theta)-len(self.x))/self.sum_y)
      if fn is not None:
         print "\nTime taken to run " + fn.__name__ + ": " \
            + str(round(self.t1-self.t0,3)) + " seconds"
      print "R_wp: " + str(R_wp)
      print "R_e: " + str(R_e)
      print "Goodness-of-Fit: " + str(R_wp/R_e)

      if not self.bkgd_refine:
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
         self.total_profile_state, label=r'$I_{\rm calc}$')
      plt.legend(bbox_to_anchor=(.8,.7))
      plt.ylabel(r"$I$")

      subplot2 = self.fig.add_subplot(312) #plt.subplot(3,1,2)
      self.pltmask = np.abs(self.two_theta - two_theta_roi) \
         < delta_theta
      plt.scatter(self.two_theta[self.pltmask],self.y[self.pltmask],
         label='Data',s=2, color='red')

      self.current_profile_masked, = subplot2.plot(self.two_theta[self.pltmask],
         self.total_profile_state[self.pltmask], label=r'$I_{\rm calc}$')
      # plt.legend(bbox_to_anchor=(.8,.7))
      plt.ylabel(r"$I$")

      subplot3 = self.fig.add_subplot(313) #plt.subplot(3,1,3)
      self.residuals, = subplot3.plot(self.two_theta,
         self.weighted_squared_errors_state,'bo',ms=2)
      plt.ylabel(r"$\frac{1}{I} \, (I-I_{\rm calc})^2$")
      plt.xlabel(r'$2\,\theta$')

      self.fig.canvas.draw()

      if autohide:
         fig.canvas.flush_events()
         time.sleep(1)
         plt.close('all')
