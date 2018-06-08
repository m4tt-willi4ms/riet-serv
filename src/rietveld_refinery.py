from __future__ import division
import os, random, math
import time
import sys, subprocess
import numpy as np
import itertools
from scipy.optimize import approx_fprime
from scipy.optimize import fmin_l_bfgs_b as minimizer
from scipy.optimize import minimize
import json,codecs
import operator

from rietveld_phases import RietveldPhases
from rietveld_plot import RietveldPlot
import src.refinement_parameters as refinement_parameters

default_factr = 1e2
default_iprint = 1
default_maxiter = 150
default_m = 10
default_pgtol = 1e-10
default_epsilon = 1e-11

# default_composition_cutoff = 3

def with_stats(wrappee):
   def fancy_minimize(self):
      print self.display_parameters()
      wrappee(self)
      print self.display_parameters(fn=wrappee)
      print self.display_stats(fn=wrappee)
   return fancy_minimize

class RietveldRefinery:
   """
      This class is used to assemble and organize the inputs required to run (a
      series of) Rietveld refinements.
   """
   def __init__(self,phase_list, rietveld_plot=None, \
      bkgd_refine=False,
      store_intermediate_state=False,
      show_plots=False,
      factr=default_factr,
      mask=None,
      iprint=default_iprint,
      maxiter=default_maxiter,
      m=default_m,
      pgtol=default_pgtol,
      epsilon=default_epsilon,
      # input_weights=None,
      # composition_cutoff=default_composition_cutoff,
      ):
      """
         Given some list of phases, an instance of the refinery is initialized
         to readily run a refinement process.
      """
      self.phase_list = phase_list
      self.rietveld_plot = rietveld_plot
      self.bkgd_refine = bkgd_refine
      self.store_intermediate_state = store_intermediate_state
      self.show_plots = show_plots
      self.factr = factr
      self.iprint = iprint
      self.maxiter = maxiter
      self.m = m
      self.pgtol = pgtol
      self.epsilon = epsilon
      # self.composition_cutoff = composition_cutoff

      self.weighted_sum_of_I_squared = np.sum(
         RietveldPhases.I**2/RietveldPhases.sigma**2)

      for key in refinement_parameters.keys:
         if key in ("labels", "values", "l_limits", "u_limits"):
            if len(self.phase_list) > 0:
               # self.x = np.hstack((RietveldPhases.global_x,
               #    np.hstack((x.phase_x for x in self.phase_list))))
               tmp = np.hstack((RietveldPhases.global_parameters.x[key],
      np.hstack((phase.phase_parameters.x[key] for phase in self.phase_list))))
            else:
               RietveldPhases.assemble_global_x()
               tmp = RietveldPhases.global_parameters.x[key]
         setattr(self, "x_" + key, tmp)
      self.x = self.x_values
         # for key in refinement_parameters.keys:
         #    setattr(self, key, RietveldPhases.global_parameters.x[key])

      if self.bkgd_refine:
         self.raw_profile_state = np.sum(x.phase_profile() for x in
            self.phase_list)

      self.global_mask = np.isin(np.array(xrange(len(self.x))),
         np.array(xrange(len(RietveldPhases.global_x))))
      self.bkgd_mask = self.make_mask(["bkgd"])
      self.global_mask_no_bkgd = np.logical_xor(
         self.global_mask, self.bkgd_mask)
      self.scale_mask = self.make_mask(["sca"])

      if mask is not None:
         self.mask = mask
      elif self.bkgd_refine:
         self.mask = self.bkgd_mask
      else:
         self.mask = np.zeros(len(self.x),dtype=bool)

      self.set_compositions()

      # self.composition_by_volume = np.zeros(len(self.phase_list))
      # if input_weights is not None:
      #    self.composition_by_weight = input_weights
      # else: self.composition_by_weight = 100*np.ones(len(self.phase_list))


      n=len(RietveldPhases.global_x)
      self.phase_masks = []
      self.global_and_phase_masks = []
      for phase in self.phase_list:
         self.phase_masks.append(np.isin(np.array(xrange(len(self.x))),
            np.array(xrange(n,n+len(phase.phase_x)))))
         self.global_and_phase_masks.append(np.logical_or(
            self.global_mask_no_bkgd, self.phase_masks[-1]))
         n+=len(phase.phase_x)

      self.update_state()
      if self.rietveld_plot is not None:
         self.rietveld_plot.setplotdata()
         self.rietveld_plot.updateplotprofile(self.total_profile_state,
            wse=self.relative_differences_state)

   def total_profile(self):
      if self.bkgd_refine:
         return RietveldPhases.background_polynomial() \
            + self.raw_profile_state
      else:
         return RietveldPhases.background_polynomial() \
            + np.sum( x.phase_profile() for x in self.phase_list)

   def weighted_squared_errors(self):
      return (RietveldPhases.I - self.total_profile_state)**2 \
         /RietveldPhases.sigma**2

   def relative_differences(self):
      return (self.total_profile_state-RietveldPhases.I) \
         /RietveldPhases.sigma**2

   def update_state(self):
      # if not self.bkgd_refine:
      RietveldPhases.global_parameters.update_x(self.x[self.global_mask],
         self.mask[self.global_mask])
      if not self.bkgd_refine:
         for phase, phase_mask in zip(self.phase_list, self.phase_masks):
               phase.phase_parameters.update_x(
                  self.x[phase_mask], self.mask[phase_mask])
               # self.x[self.phase_masks[i]][self.mask[self.phase_masks[i]]])
      self.total_profile_state = self.total_profile()
      # print "profile: " + str(self.phase_list[0].phase_x[3])
      self.relative_differences_state = self.relative_differences()
      # self.Relative_Differences_state.shape = \
      #    (self.Relative_Differences_state.shape[0],1)
      self.weighted_squared_errors_state = self.weighted_squared_errors()

   def revert_to_x(self, x):
      self.x = x
      self.update_state()
      self.set_compositions()
      if self.rietveld_plot is not None:
         self.rietveld_plot.updateplotprofile(self.total_profile_state,
               wse=self.relative_differences_state)

   def weighted_sum_of_squares(self,x):
      # print self.mask
      self.x[self.mask] = x
      # if not self.bkgd_refine:
      self.update_state()
      return np.sum(self.weighted_squared_errors_state)

   def weighted_sum_of_squares_grad(self,x):
      self.x[self.mask] = x

      result = np.zeros(len(self.x[self.mask]))
      # print "here"
      bkgd_mask = np.logical_and(self.mask, self.bkgd_mask)[self.mask]
      if np.any(bkgd_mask):
         result[bkgd_mask] = \
            np.sum(2*np.multiply(
               RietveldPhases.two_theta_powers[:len(RietveldPhases.bkgd)],
               self.relative_differences_state),axis=1)

      for i in xrange(0,len(self.phase_list)):
         mask = self.global_and_phase_refine_masks[i]
         if np.any(mask):
            result[mask[self.mask]] += \
               np.sum(2*np.multiply(self.phase_list[i].
               phase_profile_grad(mask[self.global_and_phase_masks[i]],
                  epsilon=self.epsilon),
               self.relative_differences_state),axis=1)
      # print "post-grad: " + str(self.phase_list[0].phase_x[3])
      return result

   def set_mask(self, list_of_parameter_strings=[]):
      mask = self.make_mask(list_of_parameter_strings=list_of_parameter_strings)
      if np.any(mask):
         self.mask = mask

   def make_mask(self, list_of_parameter_strings=[]):
      mask = np.zeros(len(self.x), dtype=bool)
      if len(list_of_parameter_strings) == 1:
         mask = np.char.startswith(self.x_labels, list_of_parameter_strings[0])
      elif len(list_of_parameter_strings) > 1:
         for param in list_of_parameter_strings:
            mask = np.logical_or(mask, np.char.startswith(self.x_labels, param))
      return mask

   @with_stats
   def minimize(self, callback_functions=[]):
      self.callback_functions = callback_functions
      # self.mask = np.logical_and(self.del_mask,self.mask)
      self.display_parameters()
      self.count = 0

      self.global_and_phase_refine_masks = []
      for i in xrange(len(self.phase_list)):
         self.global_and_phase_refine_masks.append(
            np.logical_and(self.mask, self.global_and_phase_masks[i]))

      if not self.bkgd_refine and self.rietveld_plot is not None:
         self.rietveld_plot.fig.suptitle("In Progress...")

      options = {
         'ftol': self.factr*2.2e-16,
         # 'xtol': self.factr*2.2e-16,
         'gtol': self.pgtol,
         'maxcor': self.m,
         'maxiter': self.maxiter,
         'disp': True,
         }

      self.t0 = time.time()

      if np.sum(self.mask) > 0:
         self.result = \
            minimize(self.weighted_sum_of_squares,self.x[self.mask],
            method = 'L-BFGS-B',#'Newton-CG',#
            options=options,
            # jac = False,
            jac = self.weighted_sum_of_squares_grad,
            callback = self.callback,
            # approx_grad = True,
            # epsilon = self.epsilon, \
            bounds = zip(self.x_l_limits[self.mask],
               self.x_u_limits[self.mask]),
            )
      else:
         self.result = None

      self.t1 = time.time()

      if self.rietveld_plot is not None and self.result is not None:
         if self.result['message'][0:4] == "STOP":
            self.rietveld_plot.fig.suptitle("Optimization Ended...")
         elif self.result['message'][0:4] == "CONV":
            self.rietveld_plot.fig.suptitle("Optimization Complete.")
         elif self.result['message'][0:4] == "ABNO":
            self.rietveld_plot.fig.suptitle("Try Again...")
      elif self.rietveld_plot is not None:
         self.rietveld_plot.fig.suptitle("No parameters to be refined.")

      self.set_compositions()

      if self.rietveld_plot is not None:
         self.rietveld_plot.updateplotprofile(self.total_profile_state,
            wse=self.relative_differences_state)

   def set_compositions(self):
      Scales = self.x[self.scale_mask]
      total = np.sum(Scales)
      weight_moments = []
      # print "\n"
      self.composition_by_volume = np.zeros(len(self.phase_list))
      for i,val in np.ndenumerate(Scales):
         density = self.phase_list[i[0]].phase_data["crystal_density"]
         weight_moments.append(val*density)
         self.composition_by_volume[i] = val/total*100
         # print "Phase " + str(i[0]+1) + ": " + str(val/total*100) + " %"

      weight_moments = np.array(weight_moments)
      weight_total = np.sum(weight_moments)
      # print "\nBy weight:"
      self.composition_by_weight = np.zeros(len(self.phase_list))
      for i,val in np.ndenumerate(weight_moments):
         self.composition_by_weight[i] = val/weight_total*100
         # print "Phase " + str(i[0]+1) + ": " + \
         #    str(val/weight_total*100)  + " %"

   def update_phase_list(self):
      Scale_mask = self.scale_mask
      Scale_vals = self.x['values'][Scale_mask][self.del_mask[Scale_mask]]
      Scale_total = np.sum(Scale_vals)
      del_list = []
      for i,Scale_val in np.ndenumerate(Scale_vals):
         if Scale_val/Scale_total < 1e-5:
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

   def callback(self, x):

      # sys.stdout.write('.')
      # sys.stdout.flush()
      # if self.phase_list:
         # print "grad:", self.weighted_sum_of_squares_grad(self.x[self.mask])
         # print "two_theta_0:", RietveldPhases.two_theta_0
         # print "scale:", self.phase_list[0].phase_parameters.scale
      if self.rietveld_plot is not None and self.count % 5 == 0:
         self.rietveld_plot.updateplotprofile(self.total_profile_state,
            wse=self.relative_differences_state)
         map(operator.methodcaller('__call__'), list(self.callback_functions))
      # print self.x[self.mask]
      self.count += 1


   def minimize_with_mask(self,mask):
      self.mask = mask
      self.minimize()

   def minimize_bkgd(self,display_plots = True):
      self.mask = self.bkgd_mask
      self.minimize()

   def minimize_all(self,display_plots = True):
      self.mask = np.ones(len(self.x),dtype=bool)
      self.minimize()

   def display_parameters(self, fn=None):
      if fn is not None:
         param_list = "After " + fn.__name__ + ":\n"
      else:
         param_list = ""
      mask = np.logical_and(self.mask,self.global_mask)
      if np.any(mask):
         param_list += "Global:\n"
      for label,value,l,u in zip(
            self.x_labels[mask],
            self.x[mask],
            self.x_l_limits[mask],
            self.x_u_limits[mask]):
            param_list += \
               "  " + label + " = " \
               + ('%.4g' % value) + " (" \
               + ('%.4g' % l) + ", " \
               + ('%.4g' % u) + ")\n"
      for i,phase_mask in enumerate(self.phase_masks):
         mask = np.logical_and(self.mask,phase_mask)
         if np.any(mask):
            param_list += "Phase " + str(i+1) + ": " + \
               self.phase_list[i].phase_settings["chemical_name"] +"\n"
            for label,value,l,u in zip(
               self.x_labels[mask],
               self.x[mask],
               self.x_l_limits[mask],
               self.x_u_limits[mask]):
               param_list += \
                  "  " + label + " = " \
                  + ('%.4g' % value) + " (" \
                  + ('%.4g' % l) + ", " \
                  + ('%.4g' % u) + ")\n"
      # print param_list
      return param_list

   def display_stats(self, fn=None):
      # self.mask = np.ones(len(self.x),dtype=bool)
      WSS = np.sum(self.weighted_squared_errors_state)
      R_wp = np.sqrt(WSS/self.weighted_sum_of_I_squared)
      R_e = np.sqrt((len(RietveldPhases.two_theta)-len(self.x[self.mask])
         )/self.weighted_sum_of_I_squared)

      self.num_params = np.sum(self.mask)
      self.GoF = R_wp/R_e

      output = ''

      if fn is not None:
         output += "\n" + self.result['message'] +"\n"
         output +=  "\nTime taken to run " + fn.__name__ + " with " \
            + str(self.num_params) + " parameters: " \
            + str(round(self.t1-self.t0,3)) + " seconds (" \
            + str(self.result['nit']) + " iterations, " \
            + str(self.result['nfev']) + " function calls)\n"
      output +=  "R_wp: " + str(R_wp) +"\n"
      output +=  "R_e: " + str(R_e) +"\n"
      output +=  "Goodness-of-Fit: " + str(self.GoF) +"\n"

      if not self.bkgd_refine:
         output +=  "\n"
         for i in xrange(len(self.phase_list)):
            output +=  "Phase " + str(i+1) + ": " \
               + ('%.3g' % self.composition_by_volume[i]) + " %\n"

         output +=  "\nBy weight:\n"
         for i in xrange(len(self.phase_list)):
            output +=  "Phase " + str(i+1) + ": " \
               + ('%.3g' % self.composition_by_weight[i]) + " %\n"
         output +=  "\n"
      # print output
      return output

   def get_plot_data(self):
      intensities = self.total_profile()
      d = {}
      d['two_thetas'] = []
      d['errors'] = []
      d['intensities'] = list(intensities)
      return d