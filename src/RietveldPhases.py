from __future__ import division
import os, random, math
import time
# import inspect
import sys
import numpy as np
from scipy.optimize import approx_fprime
import matplotlib.pyplot as plt

import paths

import iotbx.cif, cctbx.miller
from cctbx import xray
from cctbx import crystal
from cctbx import uctbx
from cctbx.array_family import flex
from cctbx.eltbx import wavelengths
import jsonpickle
from libtbx import easy_pickle
from scitbx import lbfgsb

custom_dtype = np.dtype([ \
      ('labels','S12'), \
      ('values','f8'), \
      ('l_limits','f8'), \
      ('u_limits','f8') \
      # ('round','i4')
      ])

default_bkgd_order = 3
default_two_theta_0 = np.array([('two_theta_0',0.0,-0.2,0.2)],
   dtype=custom_dtype)
default_vertical_offset = False #:False = angular offset; True = Vertical Offset

default_U = np.array([('U',0.00,0,0.1)],dtype=custom_dtype)
default_V = np.array([('V',-0.00,-0.1,0)],dtype=custom_dtype)
default_W = np.array([('W',0.001,0.0001,1)],dtype=custom_dtype)
default_Amplitude = np.array([('Amplitude',0.1,0,float('inf'))],
   dtype=custom_dtype)
default_eta_order = 2
default_lattice_dev = 0.01

uc_labels = ["a","b","c","alpha","beta","gamma"]


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

   custom_dtype = custom_dtype

   two_theta_0 = default_two_theta_0
   vertical_offset = default_vertical_offset
   bkgd = None #: Initialized to None

   max_polynom_order = 5

   lambdas=["CUA1","CUA2"] #: Default values
   K_alpha_factors = [1,0.48]  #: Default values

   @classmethod
   def set_bkgd_order(cls,order):
      assert isinstance(order,int) == True
      if order > cls.max_polynom_order:
         order = cls.max_polynom_order
      elif order < 1:
         order = 1
      cls.bkgd = np.hstack((x for x in cls.bkgd_param_gen(order)))
      return cls.bkgd

   @classmethod
   def bkgd_param_gen(cls,order=default_bkgd_order):
      n=0
      while n < order:
         # if cls.bkgd == None:
         yield np.array([('bkgd_'+str(n),0.0,-float('inf'),float('inf'))],
            dtype=custom_dtype)
         # else:
         #    yield cls.bkgd[n]
         n+=1

   @classmethod
   def set_profile(cls,filename,
      number_of_columns=3,
      min_two_theta=0,
      max_two_theta=180):
      two_theta = []
      I = []
      with open(filename) as file:
         for line in file.readlines():#[4:]:
            if number_of_columns == 2:
               two_thetatmp, Itmp = line.split()
            # if float(two_thetatmp) < 15:
               I.append(Itmp)
            elif number_of_columns == 3:
               two_thetatmp, Itmp, sigmatmp = line.split()
               I.append(float(sigmatmp)**2)
            two_theta.append(float(two_thetatmp))
      cls.two_theta = np.array(two_theta) 
      cls.I = np.array(I)

      cls.I = cls.I \
         [np.logical_and(cls.two_theta > min_two_theta,
            cls.two_theta < max_two_theta)]
      cls.two_theta = cls.two_theta \
         [np.logical_and(cls.two_theta > min_two_theta,
            cls.two_theta < max_two_theta)]

      CU_wavelength = wavelengths.characteristic(RietveldPhases.lambdas[0]) \
         .as_angstrom()
      cls.d_min = CU_wavelength/2/np.sin(np.pi/360*cls.two_theta[-1])
      cls.d_max = CU_wavelength/2/np.sin(np.pi/360*cls.two_theta[0])

      cls.two_theta_powers = np.power(cls.two_theta,np.array(
         xrange(0,cls.max_polynom_order)).reshape(cls.max_polynom_order,1))

      cls.set_bkgd_order(default_bkgd_order)

   @classmethod
   def background_polynomial(cls):
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
      dim = len(cls.bkgd)
      return np.dot(cls.bkgd['values'],cls.two_theta_powers[:dim,:])

   @classmethod
   def LP_intensity_scaling(self):
      r"""
         Computes the Lorentz-Polarization intensity scaling factors for a 
         set of two-theta values listed in ``two_theta``, via the equation

         .. math:: LP(2\theta) = \frac{1+\cos^2(2\theta)}{\sin\theta
            \,\sin(2\theta)} \,.
            :label: LPDefn

         :param two_theta: list of :math:`2\theta` positions

      """
      return (1+np.cos(math.pi/180*RietveldPhases.two_theta)**2) \
         /np.sin(math.pi/360*RietveldPhases.two_theta) \
         /np.sin(math.pi/180*RietveldPhases.two_theta)

   @classmethod
   def assemble_global_x(cls):
      cls.global_x = np.hstack((x for x in cls.global_param_gen()))
      cls.global_x_no_bkgd_mask = np.invert(
         np.char.startswith(cls.global_x['labels'],'bkgd'))
      
      cls.two_theta_0 = cls.global_x[0]
      cls.bkgd = cls.global_x[1:1+cls.bkgd.shape[0]]

   @classmethod
   def update_global_x(cls,global_x,mask=None):
      if mask is None:
         cls.global_x = global_x
      else:
         cls.global_x[mask] = global_x[mask]


   @classmethod
   def global_param_gen(cls):
      yield cls.two_theta_0
      yield cls.bkgd

   def __init__(self,fn_cif,
      I_max=None,
      delta_theta = 0.5,
      Intensity_Cutoff=0.01,
      Amplitude=default_Amplitude,
      U=default_U,
      V=default_V,
      W=default_W,
      eta_order=default_eta_order,
      lattice_dev=default_lattice_dev,
      recompute_peak_positions=False,
      ):

      if I_max is not None:
         self.I_max= I_max
      else:
         self.I_max = np.amax(RietveldPhases.I)

      self.Intensity_Cutoff = Intensity_Cutoff
      self.delta_theta = delta_theta
      self.U = U
      self.V = V
      self.W = W
      self.Amplitude = Amplitude
      self.eta = self.set_eta_order(eta_order)
      self.lattice_dev = lattice_dev
      self.recompute_peak_positions = recompute_peak_positions

      self.load_cif(fn_cif)
      self.lattice_parameters = self.set_lattice_parameters()

      self.compute_relative_intensities()

      self.Amplitude['values'] = self.Amplitude['values']* \
         self.I_max/np.amax(self.phase_profile())

   def phase_param_gen(self):
      yield self.U
      yield self.V
      yield self.W
      yield self.Amplitude
      yield self.eta
      yield self.lattice_parameters

   def assemble_lattice_params(self):
      crystal_system = self.structure.space_group().crystal_system()
      uc_params = self.unit_cell.parameters()
      if crystal_system == "Triclinic":
         mask = [True,True,True,True,True,True]
      elif crystal_system == "Monoclinic":
         mask = [True,True,True]
         for i in xrange(3,6):
            if np.isclose(uc_params[i], 90):
               mask.append(False)
            else: mask.append(True)
      elif crystal_system == "Orthorhombic":
         mask = [True,True,True,False,False,False]
      elif crystal_system == "Tetragonal":
         mask = [True]
         for i in xrange(1,3):
            if np.isclose(uc_params[i],uc_params[0]):
               mask.append(False)
            else: mask.append(True)
         mask.append([False,False,False])
      elif crystal_system == "Trigonal":
         if np.isclose(uc_params[3],uc_params[4]) and \
            np.isclose(uc_params[3],uc_params[5]):
            mask = [True,False,False,True,False,False]
         else: mask = [True,False,True,False,False,False]
      elif crystal_system == "Hexagonal":
         mask = [True,False,True,False,False,False]
      elif crystal_system == "Cubic":
         mask = [True,False,False,False,False,False]

      assert len(mask) == 6
      return np.array([x for x in self.unit_cell_parameter_gen(mask)],
            dtype=custom_dtype)

   def assemble_phase_x(self):
      self.phase_x = np.hstack((x for x in self.phase_param_gen()))
      global_x_no_bkgd = RietveldPhases.global_x \
         [RietveldPhases.global_x_no_bkgd_mask]
      self.global_and_phase_x = np.hstack((global_x_no_bkgd,self.phase_x))
      
      self.global_mask_no_bkgd = np.hstack((
         np.ones(len(global_x_no_bkgd),dtype=bool),
         np.zeros(len(self.phase_x),dtype=bool)))
      self.phase_mask = np.hstack((
         np.zeros(len(global_x_no_bkgd),dtype=bool),
         np.ones(len(self.phase_x),dtype=bool)))

      self.U = self.phase_x[0]
      self.V = self.phase_x[1]
      self.W = self.phase_x[2]
      self.Amplitude = self.phase_x[3]
      self.eta = self.phase_x[4:4+self.eta.shape[0]]
      self.lattice_parameters = self.phase_x[4+self.eta.shape[0]:
         4+self.eta.shape[0]+self.lattice_parameters.shape[0]]

   def update_params(self,phase_x,mask=None):
      if mask is None:
         self.phase_x = phase_x
      else:
         self.phase_x[mask] = phase_x[mask] 
      
      
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
      self.unit_cell_parameters = np.array(self.unit_cell.parameters())
      # print self.unit_cell_parameters

      for scatterer in self.structure.scatterers():
         if (scatterer.scattering_type == "O-2"):
            scatterer.scattering_type = "O2-"
         if (scatterer.scattering_type == "Ca+2"):
            scatterer.scattering_type = "Ca2+"
         if (scatterer.scattering_type == "Si+4"):
            scatterer.scattering_type = "Si4+"

   def unit_cell_parameter_gen(self,mask):
      """returns all unit cell parameters specified by the mask 
      (a list of six booleans)
      """
      global uc_labels
      for i in xrange(6):
         if mask[i]:
            yield ('uc_'+uc_labels[i],
               self.unit_cell_parameters[i],
               self.unit_cell_parameters[i]*(1-self.lattice_dev),
               self.unit_cell_parameters[i]*(1+self.lattice_dev))

   def set_lattice_parameters(self):
      """Returns a numpy array consisting of the lattice parameters
         
      """
      mask = [True,True,True,True,True,True]
      # mask = [False,False,False,False,False,False]
      return  np.array(
         [x for x in self.unit_cell_parameter_gen(mask)],dtype=custom_dtype)

   # def update_unit_cell(self):
   #    if self.refine_unit_cell:
   #       unit_cell_parameters = list(self.unit_cell.parameters()) #self.get_unit_cell_parameters()
   #       for (i,param) in enumerate(unit_cell_parameters):
   #          if self.unit_cell_indices[i] != 0:
   #             unit_cell_parameters[i] = \
   #                self.x['values'][self.unit_cell_indices[i]]

   #       # print unit_cell_parameters
   #       self.unit_cell = uctbx.unit_cell(tuple(unit_cell_parameters))
   #       # self.Update_two_thetas()
   #       self.Compute_Relative_Intensities()

   def compute_relative_intensities(self,lammbda="CUA1",anomalous_flag = True):
      r"""Returns squared structure factors, weighted by the multiplicity of 
         each reflection.
         
         :param float d_min: The minimum resolution needed
         :param str lammbda: Wavelength label
         
         :return: *d*-spacings and :math:`m\times|F|^2` for each reflection,
             where *m* is the multiplicity and *F* is the structure factor.
         :rtype: numpy array

      """
      self.f_miller_set = self.structure.build_miller_set(anomalous_flag, \
         d_min=RietveldPhases.d_min).sort()
      # Let's use scattering factors from the International Tables
      self.structure.scattering_type_registry(table="it1992") #,  "it1992",  
         # "wk1995" "n_gaussian"\
      self.structure.set_inelastic_form_factors( \
         photon=wavelengths.characteristic(lammbda),table="sasaki")

      f_calc =  self.structure.structure_factors(d_min=RietveldPhases.d_min, \
         anomalous_flag=anomalous_flag).f_calc().sort()

      unit_cell_volume = self.unit_cell.volume()
      f_calc_sq = f_calc.as_intensity_array().sort().data() \
         /unit_cell_volume/unit_cell_volume
      f_calc_mult = f_calc.multiplicities().sort().data()

      self.d_spacings = self.f_miller_set.d_spacings().data().as_numpy_array() 
      self.relative_intensities = f_calc_sq * f_calc_mult.as_double() \
         .as_numpy_array() #: weight intensities by the corresponding
         #: multiplicity
      # print "Volume: " + str(self.unit_cell.volume())

      # Drop any peaks below the Intensity Cutoff
      self.rel_I_max_calc = np.amax(self.relative_intensities)
      if self.relative_intensities.shape[0] != 0:
         self.d_mask = np.logical_and( \
            self.relative_intensities > self.Intensity_Cutoff
            *self.rel_I_max_calc, self.d_spacings < RietveldPhases.d_max)
         self.d_spacings = self.d_spacings \
            [self.d_mask]
         self.relative_intensities = self.relative_intensities \
            [self.d_mask]

      two_thetas = np.zeros((2,len(self.d_spacings)))
      # factors = np.zeros((2,len(self.d_spacings)))
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
         np.broadcast_to(RietveldPhases.LP_intensity_scaling(),
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

   def update_two_thetas(self):
      self.d_spacings = self.f_miller_set.d_spacings().data().as_numpy_array() \
         [self.d_mask] 

      two_thetas = np.zeros((2,len(self.d_spacings)))
      # factors = np.zeros((2,len(self.d_spacings)))
      for i in xrange(0,len(self.lambdas),1):
         # read wavelength
         lambda_i = wavelengths.characteristic(self.lambdas[i]).as_angstrom()
         # Compute two_theta for each d-spacing
         two_thetas[i] = 360/math.pi*np.arcsin(lambda_i/2/self.d_spacings)
      
      # Assemble into a single array
      self.two_theta_peaks = np.concatenate((two_thetas[0],two_thetas[1]))
      #: list of peak positions as a function of :math:`2\theta`
      self.two_theta_peaks.shape = (self.two_theta_peaks.shape[0],1)
      # self.weighted_intensities = np.concatenate( \
      #    (self.K_alpha_factors[0]*self.relative_intensities, \
      #     self.K_alpha_factors[1]*self.relative_intensities))
      # self.weighted_intensities.shape = (self.weighted_intensities.shape[0],1)
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

   def set_eta_order(self,order):
      if order > RietveldPhases.max_polynom_order:
         order = RietveldPhases.max_polynom_order
      elif order <1:
         order = 1
      return np.hstack((x for x in self.eta_param_gen(order)))

   def eta_param_gen(self,order):
      n=0
      while n < order:
         limit = np.power(0.001,n)
         if n == 0:
            yield np.array([('eta_'+str(n),0.5,0,1)],
               dtype=custom_dtype)
         else: 
            yield np.array([('eta_'+str(n),0.0,-limit,limit)],
               dtype=custom_dtype)
         n+=1

   def eta_polynomial(self):
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
      dim = self.eta.shape[0]
      return np.dot(self.eta['values'],
         RietveldPhases.two_theta_powers[:dim,:])

   def phase_profile(self):
      if self.recompute_peak_positions:
         self.update_two_thetas()
      # print "called phase_profile()", inspect.stack()[1][3]
      result = np.zeros((len(self.two_theta_peaks),len(self.two_theta)))
      # two_theta_0 = RietveldPhases.x['values'][RietveldPhases.two_theta_0_index]
      # Amplitude = RietveldPhases.x['values'][self.Amplitude_index]
      eta_vals = np.broadcast_to(self.eta_polynomial(), \
         (len(self.two_theta_peaks),len(self.two_theta)))[self.masks]

      if RietveldPhases.vertical_offset:
         two_theta_0_vals = None
      else:
         two_theta_0_vals = np.broadcast_to(
            RietveldPhases.two_theta_0['values'],
            (len(self.two_theta_peaks),len(self.two_theta)))[self.masks]

      omegaUVW_squareds = np.abs(
          self.U['values']*self.tan_two_theta_peaks_sq_masked
         +self.V['values']*self.tan_two_theta_peaks_masked 
         +self.W['values'])
      two_thetabar_squared = (self.two_theta_masked 
            -two_theta_0_vals-self.two_theta_peaks_masked)**2 \
            /omegaUVW_squareds
      result[self.masks] = self.Amplitude['values'] \
         *self.weighted_intensities_masked \
         *self.LP_factors_masked \
         *(eta_vals/(1+two_thetabar_squared) \
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

   def phase_profile_x(self,x,mask):
      # print "x: " + str(x)
      self.global_and_phase_x['values'][mask] = x
      RietveldPhases.global_x['values'][RietveldPhases.global_x_no_bkgd_mask] \
         =self.global_and_phase_x['values'][self.global_mask_no_bkgd] 
         # [mask[self.global_mask_no_bkgd]] = \
            # [mask[self.global_mask_no_bkgd]]
      # RietveldPhases.global_x['values'][RietveldPhases.global_x_no_bkgd_mask] = 3
      # print "mask: " + str(mask[self.global_mask_no_bkgd])
      # print "RHS: " +str(self.global_and_phase_x['values'][self.global_mask_no_bkgd])
      #       # [mask[self.global_mask_no_bkgd]])
      # print str(RietveldPhases.global_x['values'][RietveldPhases.global_x_no_bkgd_mask])
      self.phase_x['values'] = \
         self.global_and_phase_x['values'][self.phase_mask]
      # self.global_and_phase_x['values'][mask] = x
      # self.update_params(global_and_phase=True)
      # print "eps: " + str(self.phase_x[3])
      if self.recompute_peak_positions:
         self.update_unit_cell()
      # print sys._getframe(1).f_code.co_name
      return self.phase_profile()

   def phase_profile_grad(self, mask, epsilon=1e-6):
      num_params = np.sum(mask)
      result = np.zeros((num_params,len(RietveldPhases.two_theta)))
      epsilons = epsilon*np.identity(num_params)
      self.global_and_phase_x = np.hstack((RietveldPhases.global_x
         [RietveldPhases.global_x_no_bkgd_mask], self.phase_x))

      for i,eps in enumerate(epsilons):
         # print eps
         result[i,:] = (
             self.phase_profile_x(self.global_and_phase_x['values'][mask]+eps,
               mask)
            -self.phase_profile_x(self.global_and_phase_x['values'][mask]-eps,
               mask)
            )/epsilon
      # print np.sum(result,axis=1)
      return result

   # def PseudoVoigtProfile(self,index):
   #    result = np.zeros(len(self.two_theta[self.masks[index]]))
   #    two_theta_0 = RietveldPhases.x['values'][RietveldPhases.two_theta_0_index]
   #    Amplitude = RietveldPhases.x['values'][self.Amplitude_index]
   #    eta_vals = self.eta_Polynomial(self.two_theta[self.masks[index]])
   #    omegaUVW_squareds = np.abs(
   #       RietveldPhases.x['values'][self.U_index] 
   #          *self.tan_two_theta_peaks_sq_masked[index]
   #       +RietveldPhases.x['values'][self.V_index] 
   #          *self.tan_two_theta_peaks_masked[index]
   #       +RietveldPhases.x['values'][self.W_index])
   #    two_thetabar_squared = (self.two_theta[self.masks[index]] -two_theta_0 
   #          -self.two_theta_peaks_masked[index])**2 \
   #          /omegaUVW_squareds
   #    result = Amplitude*self.weighted_intensities_masked[index] \
   #       *self.LP_factors[self.masks[index]]*(eta_vals/(1+two_thetabar_squared) \
   #          +(1-eta_vals)*np.exp(-np.log(2)*two_thetabar_squared))
   #    return result

   def peak_masks(self,delta_theta=None):
      if delta_theta is not None:
         return np.abs(self.two_theta-RietveldPhases.two_theta_0['values']
            -self.two_theta_peaks) < delta_theta
      else:
         return np.abs(self.two_theta-RietveldPhases.two_theta_0['values']
            - self.two_theta_peaks) < self.delta_theta

   # def bkgd_mask(self,two_theta,bkgd_delta_theta):
   #    return np.any(self.peak_masks(bkgd_delta_theta) \
   #       ,axis=0)

   # def showPVProfilePlot(self,plottitle,index,autohide=True):
   #    # if autohide:
   #    plt.ion()
   #    fig = plt.figure(figsize=(6,4))
   #    # plt.subplot(3,1,1)
   #    plt.scatter(self.two_theta[self.masks[index]],self.I[self.masks[index]],
   #       label='Data',s=1,color='red')
   #    plt.title(plottitle)

   #    plt.plot(self.two_theta[self.masks[index]],self.Phase_Profile() \
   #       [self.masks[index]],label=r'Total $I_{\rm calc}$')
   #    plt.plot(self.two_theta[self.masks[index]],self.PseudoVoigtProfile(index), \
   #       label=r'$I_{\rm calc}$')

   #    plt.legend(bbox_to_anchor=(.8,.7))
   #    plt.ylabel(r"$I$")

   #    # fig, ax = plt.subplots()
   #    # plt.ioff()
   #    fig.canvas.draw()
   #    plt.show()
   #    if autohide:
   #       fig.canvas.flush_events()
   #       time.sleep(1)
   #       plt.close('all')