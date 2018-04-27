from __future__ import division
import os, random, math
import time
import sys, subprocess
import numpy as np
import itertools
import matplotlib.pyplot as plt
from scipy.optimize import approx_fprime
from scipy.optimize import fmin_l_bfgs_b as minimizer
from scipy.optimize import minimize
import json,codecs
import operator

from RietveldPhases import RietveldPhases

default_factr = 1e2
default_iprint = 1
default_maxiter = 150
default_m = 10
default_pgtol = 1e-10
default_epsilon = 1e-11

# default_composition_cutoff = 3

class RietveldRefinery:
   """
      This class is used to assemble and organize the inputs required to run (a
      series of) Rietveld refinements.
   """
   def __init__(self,phase_list, rietveld_plot, \
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
      self.m = 15
      self.pgtol = pgtol
      self.epsilon = epsilon
      # self.composition_cutoff = composition_cutoff


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
      #    self.I = self.I[bkgd_mask]

      # self.del_mask = np.ones(len(self.x['values']),dtype = bool)

      self.weighted_sum_of_I_squared = np.sum(
         RietveldPhases.I**2/RietveldPhases.sigma**2)

      RietveldPhases.assemble_global_x()
      for phase in self.phase_list:
         phase.assemble_phase_x()

      self.x = np.hstack((RietveldPhases.global_x,
         np.hstack((x.phase_x for x in self.phase_list))))

      if self.bkgd_refine:
         self.raw_profile_state = np.sum( x for x in self.phase_profiles())

      if mask is not None:
         self.mask = mask
      else:
         self.mask = np.zeros(len(self.x),dtype=bool)

      self.set_compositions()

      # self.composition_by_volume = np.zeros(len(self.phase_list))
      # if input_weights is not None:
      #    self.composition_by_weight = input_weights
      # else: self.composition_by_weight = 100*np.ones(len(self.phase_list))

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

   # def update_xs(self):
   #    if np.any(self.mask[self.global_mask_no_bkgd]):
   #       RietveldPhases.update_global_x_no_bkgd(
   #          self.x[self.global_mask_no_bkgd])
   #    if np.any(self.mask[self.bkgd_mask]):
   #       RietveldPhases.update_bkgd(self.x[self.bkgd_mask])
   #    for i in xrange(len(phase_list)):
   #       if np.any(self.mask[self.phase_masks[i]]):
   #          self.phase_list[i].update_phase_x(self.x[self.phase_masks[i]])

   # def total_profile_x(self,x):
   #    self.x['values'][self.mask] = x
   #    self.update_state()
   #    if self.bkgd_refine:
   #       return RietveldPhases.background_polynomial() \
   #          + self.raw_profile_state
   #    else:
   #       return RietveldPhases.background_polynomial() \
   #          + np.sum( x for x in self.phase_profiles())

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
      return (RietveldPhases.I - self.total_profile_state)**2 \
         /RietveldPhases.sigma**2

   def relative_differences(self):
      return (self.total_profile_state-RietveldPhases.I) \
         /RietveldPhases.sigma**2

   def update_state(self):
      # if not self.bkgd_refine:
      RietveldPhases.update_global_x(self.x[self.global_mask],
         self.mask[self.global_mask])
      if not self.bkgd_refine:
         for i in xrange(len(self.phase_list)):
            self.phase_list[i].update_params(self.x[self.phase_masks[i]],
               self.mask[self.phase_masks[i]])
               # self.x[self.phase_masks[i]][self.mask[self.phase_masks[i]]])
      self.total_profile_state = self.total_profile()
      # print "profile: " + str(self.phase_list[0].phase_x[3])
      self.relative_differences_state = self.relative_differences()
      # self.Relative_Differences_state.shape = \
      #    (self.Relative_Differences_state.shape[0],1)
      self.weighted_squared_errors_state = self.weighted_squared_errors()

   def revert_to_x(self,x):
      self.x = x
      self.update_state()
      self.set_compositions()
      self.rietveld_plot.updateplotprofile(self.total_profile_state,
            wse=self.relative_differences_state)

   def weighted_sum_of_squares(self,x):
      # print self.mask
      self.x['values'][self.mask] = x
      # if not self.bkgd_refine:
      self.update_state()
      return np.sum(self.weighted_squared_errors_state)

   def weighted_sum_of_squares_grad(self,x):
      self.x['values'][self.mask] = x

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

   def minimize(self,callback_functions=None):
      if callback_functions is not None:
         self.callback_functions = callback_functions
      self.t0 = time.time()
      # self.mask = np.logical_and(self.del_mask,self.mask)
      self.display_parameters()
      self.count = 0

      self.global_and_phase_refine_masks = []
      for i in xrange(len(self.phase_list)):
         self.global_and_phase_refine_masks.append(
            np.logical_and(self.mask,self.global_and_phase_masks[i]))
      # if self.bkgd_refine:
      #    self.raw_profile_state = self.total_profile()
      # self.result = \
      #    minimizer(self.weighted_sum_of_squares,self.x['values'][self.mask],
      #    callback = self.callback,
      #    fprime = self.weighted_sum_of_squares_grad,
      #    # approx_grad = True,
      #    factr = self.factr,
      #    iprint = self.iprint,
      #    m = self.m,
      #    maxiter = self.maxiter,
      #    pgtol = self.pgtol,
      #    # epsilon = self.epsilon, \
      #    bounds = zip(self.x['l_limits'][self.mask],
      #       self.x['u_limits'][self.mask]))

      if not self.bkgd_refine:
         self.rietveld_plot.fig.suptitle("In Progress...")

      options = {
         'ftol': self.factr*2.2e-16,
         # 'xtol': self.factr*2.2e-16,
         'gtol': self.pgtol,
         'maxcor': self.m,
         'maxiter': self.maxiter,
         'disp': True,
         }

      self.result = \
         minimize(self.weighted_sum_of_squares,self.x['values'][self.mask],
         method = 'L-BFGS-B',#'Newton-CG',#
         options=options,
         # jac = False,
         jac = self.weighted_sum_of_squares_grad,
         callback = self.callback,
         # approx_grad = True,
         # epsilon = self.epsilon, \
         bounds = zip(self.x['l_limits'][self.mask],
            self.x['u_limits'][self.mask]),
         )

      self.t1 = time.time()

      if not self.bkgd_refine:
         if self.result['message'][0:4] == "STOP":
            self.rietveld_plot.fig.suptitle("Optimization Ended...")
         elif self.result['message'][0:4] == "CONV":
            self.rietveld_plot.fig.suptitle("Optimization Complete.")
         elif self.result['message'][0:4] == "ABNO":
            self.rietveld_plot.fig.suptitle("Try Again...")

      self.set_compositions()

      self.rietveld_plot.updateplotprofile(self.total_profile_state,
         wse=self.relative_differences_state,update_view=True)
      print "\n"

   def set_compositions(self):
      Scales = self.x['values'] \
         [np.char.startswith(self.x['labels'],"Sca")]
      total = np.sum(Scales)
      weight_moments = []
      # print "\n"
      self.composition_by_volume = np.zeros(len(self.phase_list))
      for i,val in np.ndenumerate(Scales):
         weight_moments.append(val*self.phase_list[i[0]].crystal_density)
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
      Scale_mask = np.char.startswith(self.x['labels'],"Sca")
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

   def callback(self,x):

      sys.stdout.write('.')
      sys.stdout.flush()
      # if self.show_plots:
      #    self.update_plot()
      # elif self.store_intermediate_state:
      #    self.update_state()
      if not self.bkgd_refine and self.count % 5 == 0:
         self.rietveld_plot.updateplotprofile(self.total_profile_state,
            wse=self.relative_differences_state)
         map(operator.methodcaller('__call__'),list(self.callback_functions))
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

      # if ((not self.store_intermediate_state) and self.show_plots):
      #    # self.total_profile_state = self.total_profile()
      #    self.show_multiplot("After: " + fn.__name__, \
      #    two_theta_roi=32.5, \
      #    delta_theta=3, \
      #    autohide=False)
      self.display_parameters(fn)
      self.display_stats(fn)
      if self.show_plots:
         self.fig.suptitle("After: " + fn.__name__)
         plt.ioff()
         # plt.show()

   def update_plot(self):
      self.current_profile.set_ydata(self.total_profile_state)
      self.current_profile_masked.set_ydata(self.total_profile_state
         [self.pltmask])
      self.residuals.set_ydata(self.weighted_squared_errors_state)
      self.fig.canvas.draw_idle()

   def display_parameters(self,fn=None):
      if fn is not None:
         param_list = "After " + fn.__name__ + ":\n"
      else:
         param_list = ""
      mask = np.logical_and(self.mask,self.global_mask)
      if np.any(mask):
         param_list += "Global:\n"
      for label,value,l,u in zip(
            self.x['labels'][mask],
            self.x['values'][mask],
            self.x['l_limits'][mask],
            self.x['u_limits'][mask]):
            param_list += \
               "  " + label + " = " \
               + ('%.4g' % value) + " (" \
               + ('%.4g' % l) + ", " \
               + ('%.4g' % u) + ")\n"
      for i,phase_mask in enumerate(self.phase_masks):
         mask = np.logical_and(self.mask,phase_mask)
         if np.any(mask):
            param_list += "Phase " + str(i+1) + ": " + \
               self.phase_list[i].chemical_name +"\n"
            for label,value,l,u in zip(
               self.x['labels'][mask],
               self.x['values'][mask],
               self.x['l_limits'][mask],
               self.x['u_limits'][mask]):
               param_list += \
                  "  " + label + " = " \
                  + ('%.4g' % value) + " (" \
                  + ('%.4g' % l) + ", " \
                  + ('%.4g' % u) + ")\n"
      # print param_list
      return param_list

   def display_stats(self,fn=None):
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
            + str(round(self.t1-self.t0,3)) + " seconds" +"\n"
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

   def show_plot(self,plottitle,Scale_factor=1,autohide=True):
      if autohide:
         plt.ion()
      # fig = plt.figure(figsize=(9,7))
      fig = plt.figure(figsize=(10,8))
      plt.scatter(self.two_theta,self.I,label='Data',s=1, color='red')
      plt.title(plottitle)

      plt.plot(self.two_theta,Scale_factor*self.TotalProfile(), \
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
      plt.scatter(RietveldPhases.two_theta,
         RietveldPhases.I,label='Data',s=1, color='red')

      self.current_profile, = self.subplot1.plot(RietveldPhases.two_theta,
         self.total_profile_state, label=r'$I_{\rm calc}$')
      plt.legend(bbox_to_anchor=(.8,.7))
      plt.ylabel(r"$I$")

      subplot2 = self.fig.add_subplot(312) #plt.subplot(3,1,2)
      self.pltmask = np.abs(RietveldPhases.two_theta - two_theta_roi) \
         < delta_theta
      plt.scatter(RietveldPhases.two_theta[self.pltmask],
         RietveldPhases.I[self.pltmask],
         label='Data',s=2, color='red')

      self.current_profile_masked, = subplot2.plot(
         RietveldPhases.two_theta[self.pltmask],
         self.total_profile_state[self.pltmask], label=r'$I_{\rm calc}$')
      # plt.legend(bbox_to_anchor=(.8,.7))
      plt.ylabel(r"$I$")

      subplot3 = self.fig.add_subplot(313) #plt.subplot(3,1,3)
      self.residuals, = subplot3.plot(RietveldPhases.two_theta,
         self.weighted_squared_errors_state,'bo',ms=2)
      plt.ylabel(r"$(\Delta I/\sigma)^2$")
      plt.xlabel(r'$2\,\theta$')

      self.fig.canvas.draw()

      if autohide:
         fig.canvas.flush_events()
         time.sleep(1)
         plt.close('all')

# import matplotlib
# matplotlib.use("TkAgg")
# matplotlib.style.use("ggplot")
# from matplotlib.backends.backend_tkagg import \
#    FigureCanvasTkAgg, NavigationToolbar2TkAgg
# from  matplotlib.figure import Figure
# import matplotlib.animation as animation
# from matplotlib.widgets import SpanSelector

# class RietveldPlot:
#    """This class is used to group together objects associated with the
#    matplotlib plots used in displaying the results of RietveldRefinery
#    instances.
#    """
#    def __init__(self,
#       width=4,height=3, #the figure width and height, in inches
#       ):
#       self.fig = Figure(figsize=(width,height),dpi=100)
#       # self.subplot1 = self.fig.subplot2grid((3,1),(0,0),rowspan=2)
#       # self.subplot1 = self.fig.subplot2grid((3,1),(0,2),rowspan=2)
#       self.subplot1 = self.fig.add_subplot(2,1,1)
#       self.subplot2 = self.fig.add_subplot(2,1,2)

#       # def format_coord(x,y):
#       #    return 'test'

#       # self.subplot1.format_coord = format_coord
#       # self.subplot2.format_coord = format_coord

#       # self.subplot3 = self.fig.add_subplot(313)
#       # self.canvas = canvas
#       # self.span = None

#    def setplotdata(self,ROI_center=None,delta_theta=3):

#       two_theta = RietveldPhases.two_theta
#       I = RietveldPhases.I
#       if ROI_center is None:
#          ROI_center = two_theta[np.argmax(I)]
#       ROI_mask = np.abs(two_theta-ROI_center) < delta_theta
#       self.I_max = I.max()
#       self.subplot1.clear()
#       self.subplot1.scatter(two_theta,I,label='Data',s=1, color='blue')
#       # self.subplot1.legend(bbox_to_anchor=(.8,.7))
#       self.subplot1.set_ylabel(r"$I$",rotation=0)
#       self.profile1, = self.subplot1.plot(two_theta,np.zeros(len(two_theta)),
#          label=r'$I_{\rm calc}$',alpha=0.7,color='red')
#       self.subplot1.axhline(y=-self.I_max/2,
#          # xmin=two_theta[0],
#          # xmax=two_theta[-1],
#          color='black')
#       self.diff1, = self.subplot1.plot(
#          two_theta,-self.I_max/2*np.ones(len(two_theta)),
#          label=r'$\Delta I$',color='green')
#       self.subplot1.legend(bbox_to_anchor=(.8,.7))

#       self.subplot2.clear()
#       self.subplot2.scatter(two_theta,I,label='Data',s=1, color='blue')
#       self.subplot2.set_ylabel(r"$I$",rotation=0)
#       self.profile2, = self.subplot2.plot(two_theta,np.zeros(len(two_theta)),
#          alpha=0.7,color='red')
#       self.subplot2.axhline(y=-self.I_max/2,
#          # xmin=two_theta[0],
#          # xmax=two_theta[-1],
#          color='black')
#       self.diff2, = self.subplot2.plot(two_theta,
#          -self.I_max/2*np.ones(len(two_theta)),
#          label=r'$\Delta I$',color='green')

#       self.subplot2.set_xlim(two_theta[ROI_mask][0],two_theta[ROI_mask][-1])
#       self.subplot2.axes.set_xlabel(r'$2\,\theta$')
#       # self.subplot2.set_ylim(top=I[ROI_mask].max())

#       # self.subplot3.set_ylabel(r"$\Delta I$")
#       # self.subplot3.set_ylim(-self.I_max/2,self.I_max/2)


#       self.fig.canvas.show()

#    def updateplotprofile(self,profile,wse=None,update_view=False):
#       # if len(self.subplot1.axes.lines) is not 0:
#       #    for x in (self.subplot1,self.subplot2):#,self.subplot3):
#       #       for line in x.axes.lines:
#       #          line.remove()
#       self.profile1.set_ydata(profile)
#       self.profile2.set_ydata(profile)
#       if wse is not None:
#          sigma = RietveldPhases.sigma
#          self.diff1.set_ydata(-self.I_max/2+wse*sigma**2)
#          self.diff2.set_ydata(-self.I_max/2+wse*sigma**2)
#       # self.subplot1.plot(two_theta,profile,label=r'$I_{\rm calc}$',
#       #    alpha=0.7,color='red')
#       #    self.subplot1.plot(two_theta,-self.I_max/2+wse*RietveldPhases.sigma**2,
#       #       label=r'$\Delta I$',color='green')
#       # self.subplot1.legend(bbox_to_anchor=(.8,.7))

#       # self.subplot2.plot(two_theta,profile,
#       #    alpha=0.7,color='red')
#       # if wse is not None:
#       #    self.subplot2.plot(two_theta,-self.I_max/2+wse*RietveldPhases.sigma**2,
#       #       color='green')

#       # if wse is not None:
#       #    self.subplot3.plot(two_theta,wse*RietveldPhases.sigma**2,
#       #       label=r'$\Delta I$',color='green')

#       def onselect(xmin, xmax):
#          x = RietveldPhases.two_theta
#          indmin, indmax = np.searchsorted(x,
#             (xmin, xmax))
#          indmax = min(len(x) - 1, indmax)

#          thisx = RietveldPhases.two_theta[indmin:indmax]
#          # thisy = RietveldPhases.I[indmin:indmax]
#          # line2.set_data
#          # subplot2.axes.lines[0].set_data(thisx, thisy)
#          self.subplot2.set_xlim(thisx[0], thisx[-1])
#          self.subplot2.axes.set_ylim(
#             top=1.07*max(np.max(profile[indmin:indmax]),
#                np.max(RietveldPhases.I[indmin:indmax])))
#          self.fig.canvas.draw_idle()

#       self.span = SpanSelector(self.subplot1, onselect, 'horizontal',
#                     rectprops=dict(alpha=0.5, facecolor='green'))

#       if update_view:
#          self.subplot1.axes.relim()
#          self.subplot1.axes.autoscale_view()
#          self.subplot2.axes.relim()
#          self.subplot2.axes.autoscale_view()
#          # self.subplot3.axes.relim()
#          # self.subplot3.axes.autoScale_view()

#       self.fig.canvas.draw_idle()

#    def reset_plot_profile(self):
#       if len(self.subplot1.axes.lines) is not 0:
#          for x in (self.subplot1,self.subplot2):#,self.subplot3):
#             for line in x.axes.lines:
#                line.remove()
#       self.fig.suptitle("")
#       self.subplot1.axes.legend_.remove()
#       self.fig.canvas.draw_idle()