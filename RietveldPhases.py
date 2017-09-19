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

   # def __init__(self,U,V,W,eta,theta_0,Amplitude,K_alpha_2_factor):
   #    self.U = U
   #    self.V = V
   #    self.W = W
   #    self.theta_0 = theta_0
   #    self.Amplitude = Amplitude
   #    self.K_alpha_2_factor = K_alpha_2_factor

   @classmethod
   def fromfile(cls,filename):
      """Reads in a set of refinement parameters from a file
      Args:
          filename (string): the name of the file from which to read parameters
      """
      with open(filename) as file:
         cls.fromstring(file.read())

   @classmethod
   def fromstring(cls,input_string):
      """Reads in a set of refinement parameters from a string
      Args:
         input_string (string): a string containing input parameters
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
         else:
            if line.split()[0] == "Bkgd:":
               cls.Bkgd = np.zeros(int(line.split()[1]))
            if line.split()[0] == "eta:":
               cls.eta = np.zeros(int(line.split()[1]))
      # for line in inp
      # return input_string.splitlines()[0]

   def __init__(self,fn_cif,fn_or_string=""):
      if fn_or_string != "":
         if fn_or_string[-4:] == ".txt":
            RietveldPhases.fromfile(fn_or_string)
         else:
            RietveldPhases.fromstring(fn_or_string)
      self.load_cif(fn_cif)

   def load_cif(self,fn,d_min = 1.0,lammbda = "CUA1"):
      """Reads in a crystal structure, unit cell from iotbx
      
      Args:
         fn (string): The file name
         d_min (float, optional): Minimum d-spacing
         lammbda (str, optional): Wavelength label; one of:
         {"CrA1", 2.28970}, {"CrA2", 2.29361}, {"Cr", 2.2909},
         {"FeA1", 1.93604}, {"FeA2", 1.93998}, {"Fe", 1.9373},
         {"CuA1", 1.54056}, {"CuA2", 1.54439}, {"Cu", 1.5418},
         {"MoA1", 0.70930}, {"MoA2", 0.71359}, {"Mo", 0.7107},
         {"AgA1", 0.55941}, {"AgA2", 0.56380}, {"Ag", 0.5608},
         {"", 0} (c.f. eltbx/wavelengths_ext.cpp)
      
      Returns:
          None
      """
      with open(fn, 'r') as file:
         as_cif = file.read()
      self.structure = iotbx.cif.reader(
       input_string=as_cif).build_crystal_structures()[fn[0:-4]]
      self.unit_cell = self.structure.unit_cell()
      return

   def Compute_Relative_Intensities(self,d_min=1.0,lammbda="CUA1"):
      """Returns squared structure factors, weighted by the multiplicity of each reflection.
      
      Args:
          d_min (float, optional): Description
          lammbda (str, optional): Wavelength label
      
      Returns:
          Numpy array: d-spacings and :math:`m\times|F|^2` for each reflection,
          where *m* is the multiplicity and F is the structure factor.
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

      d_spacings = f_miller_set.d_spacings().data().as_numpy_array()

      return np.vstack((d_spacings,(f_calc_sq * f_calc_mult.as_double()) \
         .as_numpy_array()))

   def Background_Polynomial(self, two_theta):
      dim = len(self.Bkgd)
      powers = np.array(range(dim))
      powers.shape = (dim,1)
      # print str(two_theta)
      # print str(np.dot(x_bkgd,np.power(two_theta,powers)))
      return np.dot(x_bkgd,np.power(two_theta,powers))

   def PseudoVoigtProfile(self, two_theta):
      r"""
         Computes the *Pseudo-Voigt* profile using the function in eq. 
         :eq:`PVDef`
         
         .. math:: PV(2\theta) = \frac{\eta}{1+\overline{\Delta\theta}^2}
            +\left(1-\eta\right)2^{-\overline{\Delta\theta}^2} \,,{\rm where}
            \quad
            \overline{\Delta\theta}^2 = \frac{(2\theta-2\theta_0)^2}{\omega^2}
            :label: PVDef

      """
      tan_thetapeak = np.tan(math.pi/360.0*two_theta_calc_peak)
      omegaUVW_squared = abs(U*tan_thetapeak**2+V*tan_thetapeak+W)
      two_thetabar_squared = (two_theta-two_theta0-two_theta_calc_peak)**2 \
         /omegaUVW_squared
      eta = eta_0
      return I_res_calc*Amp*(eta/(1 \
         +two_thetabar_squared) +(1-eta)*np.exp(-np.log(2)*two_thetabar_squared))

   def Profile_Calc(x,two_theta,Rel_Peak_Intensity,delta_theta):
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