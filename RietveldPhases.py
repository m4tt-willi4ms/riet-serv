from __future__ import division
import os, random, math
import time
import sys, subprocess
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing

import paths

import iotbx.cif, cctbx.miller
from cctbx import xray
from cctbx import crystal
from cctbx.array_family import flex
from cctbx.eltbx import wavelengths
import jsonpickle
from libtbx import easy_pickle
from scitbx import lbfgsb

class RietveldPhases:
   r"""
      Used to group together methods for calculating powder profiles
      for Rietveld phases.  

      Parameters
      -----------
      fn_cif : string
         This string stores the file name of the CIF card (including the .cif 
         xtension) for the corresponding phase
      input_string_or_file_name : string
         This string contains either: 
            -  The name of a file (ending in .txt) in which the initial 
               parameters, limits will be stored
            -  The input string itself (for debugging purposes)

      Other Parameters
      -----------------
      two_theta_0 : float
         :math:`2\theta_0` is a refinable offset, used 
         to adjust for any miscalibration of detector motor positions when
         collecting data
      Bkgd : np.array
         This array contains the coefficients of a polynomial
         used to model the profile background

      :cvar float two_theta_0: :math:`2\theta_0` is a refinable offset, used 
         to adjust for any miscalibration of detector motor positions when
         collecting data
      :cvar np.array Bkgd: This array contains the coefficients of a polynomial
         used to model the profile background
      :cvar float K_alpha_2_factor: The value by which the :math:`K_{\alpha 2}`
         peaks are rescaled relative to the primary (:math:`K_{\alpha 1}`) peaks
         (not refinable at present)
      :cvar float delta_theta: :math:`\delta\theta` describes the window over
         which a given peak profile is evaluated. (In practice it's inefficient 
         to evaluate profile tails many standard deviations away from the peak, 
         so some limit is set hard-coded here.)

      Notes
      -----
      Testing out the notes functionality in numpydocs
   """

   two_theta_0 = 0 #:
   Bkgd = np.zeros(0) #: Initialized to a zero-dimensional zero array

   #: Specify a custom d-type for the numpy array which will contain the 
   #: parameters used in refinement
   custom_dtype = np.dtype([ \
      ('labels','S12'), \
      ('values','f8'), \
      ('l_limits','f8'), \
      ('u_limits','f8') \
      ])

   x = np.empty(0,dtype=custom_dtype)

   lambdas=["CUA1","CUA2"] #: Default values
   K_alpha_factors = [1,0.48]  #: Default values

   @classmethod
   def read_param_line(cls,line):
      return np.array([( \
                  line.split()[0], \
                  float(line.split()[1]), \
                  float(line.split()[2]), \
                  float(line.split()[3]) \
                  )],dtype=cls.custom_dtype)

   @classmethod
   def global_params_from_file(cls,filename):
      """
         Reads in a set of global refinement parameters from a file
         
         :param str filename: the name of the file from which to read parameters
      """
      with open(filename) as file:
         cls.global_params_from_string(file.read())

   @classmethod
   def global_params_from_string(cls,input_string,two_theta,I):
      cls.two_theta = two_theta
      cls.I = I

      cls.num_global_params = 0
      for line in input_string.splitlines():
         if line.split()[0][-1] != ':':
            if line.split()[0] == "two_theta_0":
               cls.two_theta_0_index = cls.x['values'].shape[0]
               cls.x = np.append(cls.x,cls.read_param_line(line))
               cls.num_global_params += 1
            # if line.split()[0] == "K_alpha_2_factor":
            #    cls.K_alpha_2_factor = float(line.split()[1])
         else:
            if line.split()[0] == "Bkgd:":
               assert int(line.split()[1]) > 0
               cls.Bkgd_rank = int(line.split()[1])
               cls.Bkgd_0_index = cls.x.shape[0]
               for p in xrange(0,cls.Bkgd_rank):
                  cls.x = np.append(cls.x, \
                     cls.read_param_line("Bkgd_" + str(p) +" 0.0 -inf inf"))
               cls.num_global_params += 1

   @classmethod
   def Background_Polynomial(cls, two_theta):
      r""" Returns a numpy array populated by the values of a background 
      polynomial, :math:`P(2\theta)`, with input parameters :math:`c_i` stored
      in the class variable ``RietveldPhases.x`` with label ``Bkgd_i``:

      .. math:: P(2\theta) = \sum_{i=0}^{N} c_i (2\theta)^i

      where *N* is the length of the numpy array ``RietveldPhases.Bkgd``.

      :param np.array two_theta: a list of :math:`2\theta` values at which to 
         evaluate :math:`P(2\theta)`
      :return: the values of the background polynomial at points in 
         ``two_theta``
      :rtype: np.array

      """
      mask = np.char.startswith(cls.x['labels'],"Bkgd")
      Bkgd = cls.x['values'][mask]
      dim = len(Bkgd)
      powers = np.array(range(dim))
      powers.shape = (dim,1)
      return np.dot(Bkgd,np.power(two_theta,powers))

   @classmethod
   def empty_x(self):
      RietveldPhases.x = np.empty(0,dtype=RietveldPhases.custom_dtype)

   def __init__(self,fn_cif,input_string_or_file_name,d_min,d_max, 
      I_max=None,delta_theta = 0.5,Intensity_Cutoff=0.01):

      self.load_cif(fn_cif)
      self.d_min = d_min
      self.d_max = d_max
      self.Intensity_Cutoff = Intensity_Cutoff
      self.delta_theta = delta_theta

      if I_max is not None:
         self.I_max= I_max
      else:
         self.I_max = np.amax(RietveldPhases.I)

      if input_string_or_file_name[-4:] == ".txt":
         self.params_from_file(input_string_or_file_name)
      else:
         self.params_from_string(input_string_or_file_name)

      self.LP_factors = self.LP_Intensity_Scaling(RietveldPhases.two_theta)
      self.Compute_Relative_Intensities()

      # print "two_theta_peaks_max: " + str(self.two_theta_peaks \
      #    [self.weighted_intensities.argmax()])
      # print "I_max: " + str(self.I_max)

      # print "Phase Profile Max (Before): " + \
      #    str(np.amax(self.Phase_Profile(self.two_theta_peaks[:,0])))
      # print "x (Before): " + str(RietveldPhases.x['values'] \
      #    [self.Amplitude_index]) 
      RietveldPhases.x['values'][self.Amplitude_index] =  \
         RietveldPhases.x['values'][self.Amplitude_index] * \
            self.I_max/np.amax(self.Phase_Profile())
      # print "Phase Profile Max (After): " + \
      #    str(np.amax(self.Phase_Profile(self.two_theta_peaks[:,0])))
      # print "x (After): " + str(RietveldPhases.x['values'] \
      #    [self.Amplitude_index]) 

   def params_from_file(self,filename):
      """
         Reads in a set of refinement parameters from a file
         
         :param str filename: the name of the file from which to read parameters

      """
      with open(filename) as file:
         self.params_from_string(file.read())

   def params_from_string(self,input_string):
      """
         Reads in a set of refinement parameters from a string
      
         :param str input_string: a string containing input parameters
      """
      self.num_params = 0
      self.phase_0_index = RietveldPhases.x.shape[0]
      for l,line in enumerate(input_string.splitlines()):
         if line.split()[0][-1] != ':':
            if line.split()[0] == "U":
               self.U_index = RietveldPhases.x['values'].shape[0]
               RietveldPhases.x = np.append(RietveldPhases.x, \
                  self.read_param_line(line))
               self.num_params += 1
            if line.split()[0] == "V":
               self.V_index = RietveldPhases.x['values'].shape[0]
               RietveldPhases.x = np.append(RietveldPhases.x, \
                  self.read_param_line(line))
               self.num_params += 1
            if line.split()[0] == "W":
               self.W_index = RietveldPhases.x['values'].shape[0]
               RietveldPhases.x = np.append(RietveldPhases.x, \
                  self.read_param_line(line))
               self.num_params += 1
            if line.split()[0] == "Amplitude":
               self.Amplitude_index = RietveldPhases.x['values'].shape[0]
               RietveldPhases.x = np.append(RietveldPhases.x, \
                  self.read_param_line(line))
               self.num_params += 1
            # if line.split()[0] == "K_alpha_2_factor":
            #    cls.K_alpha_2_factor = float(line.split()[1])
         else:
            if line.split()[0] == "eta:":
               assert int(line.split()[1]) > 0
               self.eta_rank = int(line.split()[1])
               self.eta_0_index = RietveldPhases.x.shape[0]
               for p in xrange(0,self.eta_rank):
                  if p == 0:
                     RietveldPhases.x = np.append(RietveldPhases.x, \
                        RietveldPhases.read_param_line( \
                           "eta_" + str(p) +" 0.5 0.0 1.0"))
                     self.num_params += 1
                  else:
                     RietveldPhases.x = np.append(RietveldPhases.x, \
                        RietveldPhases.read_param_line( \
                           "eta_" + str(p) +" 0.0 0.0 "+ str(0.001**p)))
                     self.num_params += 1

   def load_cif(self,fn,d_min = 1.0,lammbda = "CUA1"):
      """Reads in a crystal structure, unit cell from iotbx
      
         :param str fn: The file name
         :param float d_min: Minimum d-spacing
         :param str lammbda: Wavelength label; one of: {"CrA1", 2.28970}, 
            {"CrA2", 2.29361}, {"Cr", 2.2909}, {"FeA1", 1.93604}, 
            {"FeA2", 1.93998}, {"Fe", 1.9373}, {"CuA1", 1.54056}, 
            {"CuA2", 1.54439}, {"Cu", 1.5418}, {"MoA1", 0.70930}, 
            {"MoA2", 0.71359}, {"Mo", 0.7107}, {"AgA1", 0.55941}, 
            {"AgA2", 0.56380}, {"Ag", 0.5608}, {"", 0} 
            (c.f. eltbx/wavelengths_ext.cpp)

      """
      with open(fn, 'r') as file:
         as_cif = file.read()
      self.structure = iotbx.cif.reader( \
         input_string=as_cif).build_crystal_structures() \
            [os.path.split(fn)[1][0:7]]
      self.unit_cell = self.structure.unit_cell()

      for scatterer in self.structure.scatterers():
         if (scatterer.scattering_type == "O-2"):
            scatterer.scattering_type = "O2-"
         if (scatterer.scattering_type == "Ca+2"):
            scatterer.scattering_type = "Ca2+"
         if (scatterer.scattering_type == "Si+4"):
            scatterer.scattering_type = "Si4+"

   def Compute_Relative_Intensities(self,lammbda="CUA1"):
      r"""Returns squared structure factors, weighted by the multiplicity of 
         each reflection.
         
         :param float d_min: The minimum resolution needed
         :param str lammbda: Wavelength label
         
         :return: *d*-spacings and :math:`m\times|F|^2` for each reflection,
             where *m* is the multiplicity and *F* is the structure factor.
         :rtype: numpy array
      """
      anomalous_flag = True
      f_miller_set = self.structure.build_miller_set(anomalous_flag, \
         d_min=self.d_min).sort()
      # Let's use scattering factors from the International Tables
      self.structure.scattering_type_registry(table="it1992") #,  "it1992",  
         # "wk1995" "n_gaussian"\
      self.structure.set_inelastic_form_factors( \
         photon=wavelengths.characteristic(lammbda),table="sasaki")

      f_calc =  self.structure.structure_factors(d_min=self.d_min, \
         anomalous_flag=anomalous_flag).f_calc().sort()

      unit_cell_volume = self.unit_cell.volume()
      f_calc_sq = f_calc.as_intensity_array().sort().data() \
         /unit_cell_volume/unit_cell_volume
      f_calc_mult = f_calc.multiplicities().sort().data()

      self.d_spacings = f_miller_set.d_spacings().data().as_numpy_array() 
      self.relative_intensities = f_calc_sq * f_calc_mult.as_double() \
         .as_numpy_array() #: weight intensities by the corresponding
         #: multiplicity
      # print "Volume: " + str(self.unit_cell.volume())

      # Drop any peaks below the Intensity Cutoff
      rel_I_max_calc = np.amax(self.relative_intensities)
      if self.relative_intensities.shape[0] != 0:
         mask = np.logical_and( \
            self.relative_intensities > self.Intensity_Cutoff*rel_I_max_calc, \
            self.d_spacings < self.d_max)
         self.d_spacings = self.d_spacings \
            [mask]
         self.relative_intensities = self.relative_intensities \
            [mask]

      two_thetas = np.zeros((2,len(self.d_spacings)))
      factors = np.zeros((2,len(self.d_spacings)))
      for i in xrange(0,len(self.lambdas),1):
         # read wavelength
         lambda_i = wavelengths.characteristic(self.lambdas[i]).as_angstrom()
         # Compute two_theta for each d-spacing
         two_thetas[i] = 360/math.pi*np.arcsin(lambda_i/2/self.d_spacings)
      
      # Assemble into a single array
      self.two_theta_peaks = np.concatenate((two_thetas[0],two_thetas[1]))
      #: list of peak positions as a function of :math:`2\theta`
      self.two_theta_peaks.shape = (self.two_theta_peaks.shape[0],1)
      self.weighted_intensities = np.concatenate( \
         (self.K_alpha_factors[0]*self.relative_intensities, \
          self.K_alpha_factors[1]*self.relative_intensities))
      self.weighted_intensities.shape = (self.weighted_intensities.shape[0],1)
      self.tan_two_theta_peaks = np.tan(math.pi/360.0*self.two_theta_peaks)
      self.tan_two_theta_peaks.shape = (self.tan_two_theta_peaks.shape[0],1)
      self.masks = self.peak_masks()

      self.two_theta_masked = \
         np.broadcast_to(self.two_theta,
            (len(self.two_theta_peaks),len(self.two_theta)))[self.masks]
      self.LP_factors_masked = \
         np.broadcast_to(self.LP_factors,
            (len(self.two_theta_peaks),len(self.two_theta)))[self.masks]
      self.weighted_intensities_masked = \
         np.broadcast_to(self.weighted_intensities,
            (len(self.two_theta_peaks),len(self.two_theta)))[self.masks]
      self.two_theta_peaks_masked = \
         np.broadcast_to(self.two_theta_peaks,
            (len(self.two_theta_peaks),len(self.two_theta)))[self.masks]
      self.tan_two_theta_peaks_masked = \
         np.broadcast_to(self.tan_two_theta_peaks,
            (len(self.two_theta_peaks),len(self.two_theta)))[self.masks]
      self.tan_two_theta_peaks_sq_masked = self.tan_two_theta_peaks_masked**2 

   def LP_Intensity_Scaling(self,two_theta):
      r"""
         Computes the Lorentz-Polarization intensity scaling factors for a 
         set of two-theta values listed in ``two_theta``, via the equation

         .. math:: LP(2\theta_{\rm peak}) = \frac{1+\cos^2(2\theta)}{\sin\theta
            \,\sin(2\theta)} \,.
            :label: LPDefn

         :param two_theta: list of :math:`2\theta` positions

      """
      return (1+np.cos(math.pi/180*two_theta)**2) \
         /np.sin(math.pi/360*two_theta) \
         /np.sin(math.pi/180*two_theta)

   def eta_Polynomial(self, two_theta):
      r""" Returns a numpy array populated by the values of the eta 
      polynomial, :math:`\eta(2\theta)`, with input parameters :math:`\eta_i` 
      stored in the class variable ``RietveldPhases.eta``:

      .. math:: P(2\theta) = \sum_{i=0}^{M} \eta_i (2\theta)^i

      where *M* is the length of the numpy array ``RietveldPhases.eta``.

      :param np.array two_theta: a list of :math:`(2\theta)` values
      :return: the values of the eta polynomial at points in 
         ``two_theta``
      :rtype: np.array

      """
      mask = np.isin(np.array(range(0,len(RietveldPhases.x))), \
         np.array(range(self.eta_0_index,self.eta_0_index+self.eta_rank)))
      eta = RietveldPhases.x['values'][mask]
      dim = len(eta)
      powers = np.array(range(dim))
      powers.shape = (dim,1)
      # print "two_theta_shape: " + str(two_theta.shape)
      # powers = np.tile(powers,(1,two_theta.shape[1]))
      return np.dot(eta,np.power(two_theta,powers))

   # def PseudoVoigtProfile(self, two_theta,two_theta_peak,weighted_intensity,
   #       tan_two_theta_peak,eta,LP_factor):
   #    r"""Computes the *Pseudo-Voigt* profile using the function 
   #    in eq. :eq:`PVDefn`:
      
   #    .. math:: PV(2\theta) = \frac{\eta}{1+\overline{\Delta\theta}^2}
   #       +\left(1-\eta\right)2^{-\overline{\Delta\theta}^2}\,, \quad{\rm where}
   #       \quad
   #       \overline{\Delta\theta}^2 
   #       := \frac{(2\theta-2\theta_0-2\theta_{\rm peak})^2}{\omega^2}
   #       :label: PVDefn

   #    and where

   #    .. math:: \omega := \left| U\,\tan^2\theta_{\rm peak}
   #       +V\,\tan\theta_{\rm peak}+W\right|

   #    is the Caglioti equation, describing the variation of peak width as a
   #    function of the Bragg angle, :math:`\theta_{\rm peak}`. (In eq. 
   #    :eq:`PVDefn`, :math:`2\theta_0` describes a refinable offset.)
   #    """

   #    return Amplitude*weighted_intensity*LP_factor* \
   #       (eta/(1 +two_thetabar_squared) +(1-eta) \
   #       *np.exp(-np.log(2)*two_thetabar_squared))

   def Phase_Profile(self):
      result = np.zeros((len(self.two_theta_peaks),len(self.two_theta)))
      two_theta_0 = RietveldPhases.x['values'][RietveldPhases.two_theta_0_index]
      Amplitude = RietveldPhases.x['values'][self.Amplitude_index]
      eta_vals = np.broadcast_to(self.eta_Polynomial(self.two_theta), \
         (len(self.two_theta_peaks),len(self.two_theta)))[self.masks]
      omegaUVW_squareds = np.abs(
         RietveldPhases.x['values'][self.U_index] 
            *self.tan_two_theta_peaks_sq_masked
         +RietveldPhases.x['values'][self.V_index] 
            *self.tan_two_theta_peaks_masked 
         +RietveldPhases.x['values'][self.W_index])
      two_thetabar_squared = (self.two_theta_masked -two_theta_0 
            -self.two_theta_peaks_masked)**2 \
            /omegaUVW_squareds
      result[self.masks] = Amplitude*self.weighted_intensities_masked \
         *self.LP_factors_masked*(eta_vals/(1+two_thetabar_squared) \
            +(1-eta_vals)*np.exp(-np.log(2)*two_thetabar_squared))
      # for i in xrange(0,len(self.two_theta_peaks)):
      #    two_thetabar_squared = (self.two_theta[self.masks[i]] -two_theta_0 
      #       -self.two_theta_peaks[i])**2 \
      #       /omegaUVW_squareds[i]
      #    eta = eta_vals[self.masks[i]]
      #    result[self.masks[i]] += Amplitude*self.weighted_intensities[i]* \
      #       self.LP_factors[self.masks[i]]* \
      #    (eta/(1 +two_thetabar_squared) +(1-eta) \
      #    *np.exp(-np.log(2)*two_thetabar_squared))
      return np.sum(result,axis=0)

   def PseudoVoigtProfile(self,index):
      result = np.zeros(len(self.two_theta[self.masks[index]]))
      two_theta_0 = RietveldPhases.x['values'][RietveldPhases.two_theta_0_index]
      Amplitude = RietveldPhases.x['values'][self.Amplitude_index]
      eta_vals = self.eta_Polynomial(self.two_theta[self.masks[index]])
      omegaUVW_squareds = np.abs(
         RietveldPhases.x['values'][self.U_index] 
            *self.tan_two_theta_peaks_sq_masked[index]
         +RietveldPhases.x['values'][self.V_index] 
            *self.tan_two_theta_peaks_masked[index]
         +RietveldPhases.x['values'][self.W_index])
      two_thetabar_squared = (self.two_theta[self.masks[index]] -two_theta_0 
            -self.two_theta_peaks_masked[index])**2 \
            /omegaUVW_squareds
      result = Amplitude*self.weighted_intensities_masked[index] \
         *self.LP_factors[self.masks[index]]*(eta_vals/(1+two_thetabar_squared) \
            +(1-eta_vals)*np.exp(-np.log(2)*two_thetabar_squared))
      return result

   def peak_masks(self,delta_theta=None):
      if delta_theta is not None:
         return np.abs(self.two_theta-self.two_theta_peaks) < delta_theta
      else:
         return np.abs(self.two_theta-self.two_theta_peaks) < self.delta_theta

   def bkgd_mask(self,two_theta,bkgd_delta_theta):
      return np.any(self.peak_masks(bkgd_delta_theta) \
         ,axis=0)

   def showPVProfilePlot(self,plottitle,index,autohide=True):
      # if autohide:
      plt.ion()
      fig = plt.figure(figsize=(6,4))
      # plt.subplot(3,1,1)
      plt.scatter(self.two_theta[self.masks[index]],self.I[self.masks[index]],
         label='Data',s=1,color='red')
      plt.title(plottitle)

      plt.plot(self.two_theta[self.masks[index]],self.Phase_Profile() \
         [self.masks[index]],label=r'Total $I_{\rm calc}$')
      plt.plot(self.two_theta[self.masks[index]],self.PseudoVoigtProfile(index), \
         label=r'$I_{\rm calc}$')

      plt.legend(bbox_to_anchor=(.8,.7))
      plt.ylabel(r"$I$")

      # fig, ax = plt.subplots()
      plt.ioff()
      plt.show()
      if autohide:
         fig.canvas.flush_events()
         time.sleep(1)
         plt.close('all')