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

class RietveldPhases:
   """
      Used to group together methods for calculating powder profiles
      for Rietveld phases.  
   """
   U = 0
   V = 0
   W = 0
   Amplitude = 0
   two_theta_0 = 0
   Bkgd = np.zeros(1)
   eta = np.zeros(1)
   K_alpha_2_factor = 0.48
   delta_theta = 0.5

   # def __init__(self,U,V,W,eta,theta_0,Amplitude,K_alpha_2_factor):
   #    self.U = U
   #    self.V = V
   #    self.W = W
   #    self.theta_0 = theta_0
   #    self.Amplitude = Amplitude
   #    self.K_alpha_2_factor = K_alpha_2_factor

   @classmethod
   def fromfile(cls,filename):
      """
         Reads in a set of global refinement parameters from a file
         
         :param str filename: the name of the file from which to read parameters

      """
      with open(filename) as file:
         cls.fromstring(file.read())

   @classmethod
   def fromstring(cls,input_string):
      """
         Reads in a set of global refinement parameters from a string
      
         :param str input_string: a string containing input parameters
      """
      # print input_string
      # print input_string.splitlines()[0].split()
      for line in input_string.splitlines():
         if line.split()[0][-1] != ':':
            if line.split()[0] == "U":
               cls.U = float(line.split()[1])
            if line.split()[0] == "V":
               cls.V = float(line.split()[1])
            if line.split()[0] == "W":
               cls.W = float(line.split()[1])
            if line.split()[0] == "Amplitude":
               cls.Amplitude = float(line.split()[1])
            if line.split()[0] == "two_theta_0":
               cls.two_theta_0 = float(line.split()[1])
            # if line.split()[0] == "K_alpha_2_factor":
            #    cls.K_alpha_2_factor = float(line.split()[1])
         else:
            if line.split()[0] == "Bkgd:":
               cls.Bkgd = np.zeros(int(line.split()[1]))
            if line.split()[0] == "eta:":
               cls.eta = np.zeros(int(line.split()[1]))
      # for line in inp
      # return input_string.splitlines()[0]

   @classmethod
   def Background_Polynomial(cls, two_theta):
      dim = len(cls.Bkgd)
      powers = np.array(range(dim))
      powers.shape = (dim,1)
      # print str(two_theta)
      # print str(np.dot(x_bkgd,np.power(two_theta,powers)))
      return np.dot(cls.Bkgd,np.power(two_theta,powers))

   def __init__(self,fn_cif,common_params_fn_or_string=""):
      if common_params_fn_or_string != "":
         if common_params_fn_or_string[-4:] == ".txt":
            RietveldPhases.fromfile(common_params_fn_or_string)
         else:
            RietveldPhases.fromstring(common_params_fn_or_string)
      self.load_cif(fn_cif)

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
      self.structure = iotbx.cif.reader(
       input_string=as_cif).build_crystal_structures()[fn[0:-4]]
      self.unit_cell = self.structure.unit_cell()

   def Compute_Relative_Intensities(self,d_min=1.0,lammbda="CUA1"):
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
         d_min=d_min).sort()
      # Let's use scattering factors from the International Tables
      self.structure.scattering_type_registry(table="it1992"), # "it1992",  
         # "wk1995" "n_gaussian"\
      self.structure.set_inelastic_form_factors( \
         photon=wavelengths.characteristic(lammbda),table="sasaki")

      f_calc =  self.structure.structure_factors(d_min=d_min, \
         anomalous_flag=anomalous_flag).f_calc().sort()

      f_calc_sq = f_calc.as_intensity_array().sort().data() 
      f_calc_mult = f_calc.multiplicities().sort().data()

      self.d_spacings = f_miller_set.d_spacings().data().as_numpy_array()
      self.relative_intensities = f_calc_sq * f_calc_mult.as_double() \
         .as_numpy_array()

      return np.vstack((self.d_spacings, self.relative_intensities))

   def Compile_Weighted_Peak_Intensities(self,d_spacings = np.zeros(1), \
      relative_intensities = np.zeros(1),lambdas=["CUA1","CUA2"]):
      r"""
         Takes a list of d-spacings and relative intensities (by default, those
         belonging to the instance) and compiles a list of peak positions, 
         labeled by their two_theta value

         :param np.array d_spacings: a list of *d*-spacings
         :param np.array relative_intensities: a list of relative 
            intensities
      """
      # Set to the instance values if none are passed
      if d_spacings == np.zeros(1):
         d_spacings = self.d_spacings
      if relative_intensities == np.zeros(1):
         relative_intensities = self.relative_intensities

      K_alpha_factors = [1,self.K_alpha_2_factor]
      two_thetas = np.zeros((2,len(d_spacings)))
      factors = np.zeros((2,len(d_spacings)))
      for i in xrange(0,len(lambdas),1):
         # read wavelength
         lambda_i = wavelengths.characteristic(lambdas[i]).as_angstrom()
         # Compute two_theta for each d-spacing
         two_thetas[i] = 180/math.pi*np.arcsin(lambda_i/2/d_spacings)
         # Compute the corresponding LP intensity weights for each peak
         factors[i] = self.LP_Intensity_Scaling(two_thetas[i], \
            K_alpha_factors[i])
      
      # Assemble into a single array
      self.two_theta_peaks = np.concatenate((two_thetas[0],two_thetas[1]))
      self.weighted_intensities = np.concatenate( \
         (factors[0]*relative_intensities,factors[1]*relative_intensities))

   def LP_Intensity_Scaling(self,two_theta_peaks,K_alpha_factor):
      r"""
         Computes the intensity scaling factors for a set of peaks at locations 
         listed in ``two_theta_peaks``, via the equation in :eq:`LPDefn`:

         .. math:: LP(2\theta) = K_{\alpha}\,\frac{1+\cos^2(2\theta)}
            {\sin\theta\,\sin(2\theta)}
            :label: LPDefn

         where :math:`K_{\alpha} \in \{K_{\alpha 1},K_{\alpha 2}\}`, as 
         appropriate.

         :param np.array two_theta_peaks: list of :math:`2\theta` peak locations
         :param float K_alpha_factor: the :math:`K_\alpha` factor (either
            1 or 0.48, by default)

      """
      return K_alpha_factor*abs((1+np.cos(math.pi/180*two_theta_peaks)**2) \
         /np.sin(math.pi/360*two_theta_peaks) \
         /np.sin(math.pi/180*two_theta_peaks))

   def PseudoVoigtProfile(self, two_theta):
      r"""
         Computes the *Pseudo-Voigt* profile using the function 
         in eq. :eq:`PVDefn`:
         
         .. math:: PV(2\theta) = \frac{\eta}{1+\overline{\Delta\theta}^2}
            +\left(1-\eta\right)2^{-\overline{\Delta\theta}^2}\,, \quad{\rm where}
            \quad
            \overline{\Delta\theta}^2 
            := \frac{(2\theta-2\theta_0-2\theta_{\rm peak})^2}{\omega^2}
            :label: PVDefn

         and where

         .. math:: \omega := \left| U\,\tan^2\theta_{\rm peak}
            +V\,\tan\theta_{\rm peak}+W\right|

         is the Cagliotti equation, describing the variation of peak width as a
         function of the Bragg angle. (In eq. :eq:`PVDefn`, :math:`2\theta_0`
         describes a refinable offset.)

      """
      tan_thetapeak = np.tan(math.pi/360.0*two_theta_calc_peak)
      omegaUVW_squared = abs(self.U*tan_thetapeak**2+self.V*tan_thetapeak \
         +self.W)
      two_thetabar_squared = (two_theta-two_theta0-two_theta_calc_peak)**2 \
         /omegaUVW_squared
      eta = eta_0
      return I_res_calc*Amp*(eta/(1 +two_thetabar_squared) \
         +(1-eta)*np.exp(-np.log(2)*two_thetabar_squared))

   def Profile_Calc(self, two_theta):
      # print Rel_Peak_Intensity
      two_theta_peaks = Rel_Peak_Intensity[0,:]
      # print two_theta_peaks
      two_theta_peaks.shape = (Rel_Peak_Intensity.shape[1],1)
      # print two_theta_peaks
      Intensities = Rel_Peak_Intensity[1,:]
      # print Intensities
      Intensities.shape = (Rel_Peak_Intensity.shape[1],1)
      # print Intensities
      mask = np.abs(two_theta - two_theta_peaks) < delta_theta
      # print delta_theta
      # print two_theta
      # print two_theta[mask[0]]
      # print two_theta_peaks[0]       
      result = np.zeros(len(two_theta))
      for i in xrange(0,Rel_Peak_Intensity.shape[1],1):
         result[mask[i]] += PseudoVoigtProfile(
             x[0], # eta
            x[1], # two_theta_0
            x[2], # U
            x[3], # V
            x[4], # W
            x[5], # Amplitude
            x[10], # eta_1
            x[11], # eta_2
            two_theta[mask[i]],two_theta_peaks[i],Intensities[i])
      return result+Background_Polynomial(two_theta,np.array([x[6],x[7],x[8]]))