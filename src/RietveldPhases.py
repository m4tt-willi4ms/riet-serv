from __future__ import division
import os
import random
import math
import time
import inspect
import sys
import numpy as np
from scipy.optimize import approx_fprime
import matplotlib.pyplot as plt
import jsonpickle

import src.paths

import iotbx.cif
from cctbx import xray, miller, crystal, uctbx
from cctbx.eltbx import wavelengths
from libtbx import easy_pickle
from scitbx import lbfgsb

custom_dtype = np.dtype([
   ('labels', 'S12'),
   ('values', 'f8'),
   ('l_limits', 'f8'),
   ('u_limits', 'f8'),
   # ('round','i4')
   ])

default_bkgd_order = 3
default_two_theta_0 = np.array([('two_theta_0', 0.0, -0.1, 0.1)],
                               dtype=custom_dtype)
default_vertical_offset = False #:False = angular offset; True = Vertical Offset

default_U = np.array([('U', 0.00, -0.1, 0.1)], dtype=custom_dtype)
default_V = np.array([('V', 0.00, -0.1, 0.1)], dtype=custom_dtype)
default_W = np.array([('W', 0.003, 0.000001, 1)], dtype=custom_dtype)
default_Amplitude = np.array([('Amplitude', 0.1, 0, float('inf'))],
                             dtype=custom_dtype)
default_delta_theta = 0.5
default_intensity_cutoff = 0.01
default_eta_order = 2
default_lattice_dev = 0.01
default_recompute_peak_positions = True

class RietveldPhases:
   r"""
      Used to collect parameters and methods for calculating powder profiles
      for Rietveld phases.

      Parameters
      -----------
      fn_cif : string
         This string stores the location of the CIF card (including the .cif
         xtension) for the corresponding phase, either as an absolute path
         or a relative one (relative to the root directory)

      I_max : float, optional
         Can be used to specify the maximum intensity relative to which the
         computed intensities should be scaled. If unspecified, the maximum
         intensity is determined from profile data (which can be loaded
         using the :func:`~src.RietveldPhases.RietveldPhases.set_profile()`
         class method, described later).

      delta_theta : float, optional
         :math:`\Delta\theta` specifies the region around which each peak
         profile is generated, using the formula

         .. math:: |\theta-\theta_{\rm peak}| < \Delta \theta \,.

         The default value is

         .. literalinclude:: ../src/RietveldPhases.py
            :lines: 40

      intensity_cutoff : float, optional
         The relative intensity, below which peaks are not generated. (In
         practice this is implemented when computing the squares of structure
         factors, and before any Lorentz, polarization rescaling is applied.)
         The default value of :math:`|F|^2_{\rm cutoff}` is

         .. literalinclude:: ../src/RietveldPhases.py
            :lines: 41


      lattice_dev : float, optional
         This parameter specifices the maximum allowed relative deviation of
         any lattice parameters (assuming lattice parameters will be refined).
         The default value of `lattice_dev` is 0.01.

         .. literalinclude:: ../src/RietveldPhases.py
            :lines: 43

      recompute_peak_positions : bool, optional
         This boolean variable determines whether or not peak positions are
         recomputed before determining the phase profile. (This is only
         necessary if lattice parameters are being refined.) The default value
         is

         .. literalinclude:: ../src/RietveldPhases.py
            :lines: 44

      Attributes
      ----------
      Amplitude : np.array (custom dtype)
         The initial input parameters for the phase scale factor (Amplitude).
         Its default label, value, lower- and upper-limit are set to be

         .. literalinclude:: ../src/RietveldPhases.py
            :lines: 38,39

         respectively.

      U : np.array (custom dtype)
         The initial input parameters for the Caglioti `U` parameter.
         Its default label, value, lower- and upper-limit are set to be

         .. literalinclude:: ../src/RietveldPhases.py
            :lines: 35

         respectively.

      V : np.array (custom dtype)
         The initial input parameters for the Caglioti `V` parameter.
         Its default label, value, lower- and upper-limit are set to be

         .. literalinclude:: ../src/RietveldPhases.py
            :lines: 36

         respectively.

      W : np.array (custom dtype)
         The initial input parameters for the Caglioti `W` parameter.
         Its default label, value, lower- and upper-limit are set to be

         .. literalinclude:: ../src/RietveldPhases.py
            :lines: 37

         respectively.

      eta_order : int
         `eta_order` is used to specify the number of parameters to appear in
         the corresponding eta polynomial (for more information, see
         :func:`~src.RietveldPhases.RietveldPhases.eta_polynomial()`).
         The default value is

         .. literalinclude:: ../src/RietveldPhases.py
            :lines: 42

      Notes
      -----
      This class contains many classmethods, which can be used by any of the
      RietveldPhases instances.

   """
   custom_dtype = custom_dtype

   bkgd_order = default_bkgd_order

   Amplitude = default_Amplitude
   U = default_U
   V = default_V
   W = default_W
   eta_order = default_eta_order

   two_theta_0 = default_two_theta_0
   vertical_offset = default_vertical_offset

   max_polynom_order = 5
   '''The maximum number of parameters allowed in any parameter represented
      as a polynomial (e.g. bkgd, eta)'''

   lambdas = ["CUA1", "CUA2"] #: Default values
   K_alpha_factors = [1, 0.48]  #: Default values

   @classmethod
   def set_bkgd_order(cls, order):
      assert isinstance(order, int) == True
      if order > cls.max_polynom_order:
         order = cls.max_polynom_order
      elif order < 1:
         order = 1
      cls.bkgd_order = order
      cls.bkgd = np.hstack((x for x in cls.bkgd_param_gen(order)))
      return cls.bkgd

   @classmethod
   def set_vertical_offset(cls, value):
      assert type(value) == bool
      cls.vertical_offset = value
      if value:
         cls.cos_theta = -360/np.pi*np.cos(np.pi/360*cls.two_theta)

   @classmethod
   def bkgd_param_gen(cls, order=default_bkgd_order):
      n = 0
      while n < order:
         # if cls.bkgd == None:
         yield np.array([('bkgd_'+str(n), 0.0, -float('inf'), float('inf'))],
                        dtype=custom_dtype)
         # else:
         #    yield cls.bkgd[n]
         n += 1

   @classmethod
   def set_profile(cls, filename,
                   number_of_columns=3,
                   min_two_theta=0,
                   max_two_theta=180,
                   wavelength=None):
      two_theta = []
      I = []
      sigma = []
      with open(filename) as file:
         for line in file.readlines()[1:]:#[4:]:
            if number_of_columns == 2:
               two_thetatmp, Itmp = line.split()
            # if float(two_thetatmp) < 15:
               I.append(float(Itmp))
               sigma.append(np.sqrt(float(Itmp)))
            elif number_of_columns == 3:
               two_thetatmp, Itmp, sigmatmp = line.split()
               # I.append(float(sigmatmp)**2)
               I.append(float(Itmp))
               sigma.append(float(sigmatmp))
            two_theta.append(float(two_thetatmp))
      cls.two_theta = np.array(two_theta)
      cls.I = np.array(I)
      cls.sigma = np.array(sigma)

      min_max_mask = np.logical_and(cls.two_theta > min_two_theta,
                                    cls.two_theta < max_two_theta)
      cls.I = cls.I[min_max_mask]
      cls.sigma = cls.sigma[min_max_mask]
      cls.two_theta = cls.two_theta[min_max_mask]

      CU_wavelength = wavelengths.characteristic(RietveldPhases.lambdas[0]) \
         .as_angstrom()
      cls.d_min = CU_wavelength/2/np.sin(np.pi/360*cls.two_theta[-1])
      cls.d_max = CU_wavelength/2/np.sin(np.pi/360*cls.two_theta[0])

      cls.two_theta_powers = np.power(cls.two_theta, np.array(
         xrange(0, cls.max_polynom_order)).reshape(cls.max_polynom_order, 1))

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
      return np.dot(cls.bkgd['values'], cls.two_theta_powers[:dim, :])

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
      two_theta = RietveldPhases.two_theta \
         - RietveldPhases.two_theta_0['values']
      return (1+np.cos(math.pi/180*two_theta)**2) \
         /np.sin(math.pi/360*two_theta) \
         /np.sin(math.pi/180*two_theta)

   @classmethod
   def assemble_global_x(cls):
      cls.global_x = np.hstack((x for x in cls.global_param_gen()))
      cls.global_x_no_bkgd_mask = np.invert(
         np.char.startswith(cls.global_x['labels'], 'bkgd'))

      cls.two_theta_0 = cls.global_x[0]
      cls.bkgd = cls.global_x[1:1+cls.bkgd.shape[0]]

   @classmethod
   def update_global_x(cls, global_x, mask=None):
      if mask is None:
         cls.global_x = global_x
      else:
         cls.global_x[mask] = global_x[mask]


   @classmethod
   def global_param_gen(cls):
      yield cls.two_theta_0
      yield cls.bkgd

   def __init__(self, fn_cif,
                I_max=None,
                delta_theta=default_delta_theta,
                intensity_cutoff=default_intensity_cutoff,
                lattice_dev=default_lattice_dev,
                recompute_peak_positions=default_recompute_peak_positions,
               ):

      self.fn_cif = fn_cif

      if I_max is not None:
         self.I_max = I_max
      else:
         self.I_max = np.amax(RietveldPhases.I)

      self.intensity_cutoff = intensity_cutoff
      self.delta_theta = delta_theta
      self.eta = self.set_eta_order(self.eta_order)
      self.lattice_dev = lattice_dev

      self.recompute_peak_positions = recompute_peak_positions

      self.load_cif(fn_cif, d_min=RietveldPhases.d_min)

      # self.uc_mask = [True,True,True,True,True,True,]
      # self.lattice_parameters = self.set_lattice_parameters()
      RietveldPhases.assemble_global_x()
      self.assemble_phase_x()

      self.compute_relative_intensities()

      self.Amplitude['values'] = self.Amplitude['values']* \
         self.I_max/np.amax(self.phase_profile())

   def phase_param_gen(self):
      yield self.U
      yield self.V
      yield self.W
      yield self.Amplitude
      yield self.eta
      if self.recompute_peak_positions:
         yield self.lattice_parameters

   def assemble_phase_x(self):
      if self.recompute_peak_positions:
         self.assemble_lattice_params()
      self.phase_x = np.hstack((x for x in self.phase_param_gen()))

      global_x_no_bkgd = RietveldPhases.global_x \
         [RietveldPhases.global_x_no_bkgd_mask]
      self.global_and_phase_x = np.hstack((global_x_no_bkgd, self.phase_x))

      self.global_mask_no_bkgd = np.hstack((
         np.ones(len(global_x_no_bkgd), dtype=bool),
         np.zeros(len(self.phase_x), dtype=bool)))
      self.phase_mask = np.hstack((
         np.zeros(len(global_x_no_bkgd), dtype=bool),
         np.ones(len(self.phase_x), dtype=bool)))

      self.U = self.phase_x[0]
      self.V = self.phase_x[1]
      self.W = self.phase_x[2]
      self.Amplitude = self.phase_x[3]
      self.eta = self.phase_x[4:4+self.eta.shape[0]]
      self.lattice_parameters = self.phase_x[4+self.eta.shape[0]:
                                             4+self.eta.shape[0]
                                             +self.lattice_parameters.shape[0]]

   def update_params(self, phase_x, mask=None):
      if mask is None:
         self.phase_x = phase_x
      else:
         self.phase_x[mask] = phase_x[mask]
      if self.recompute_peak_positions:
         self.update_unit_cell()

   def update_unit_cell(self):
      if np.char.startswith(self.crystal_system, "Tri"):
         self.unit_cell = uctbx.unit_cell(
            (float(x) for x in np.nditer(self.lattice_parameters['values'])))
      elif np.char.startswith(self.crystal_system, "M"):
         a = self.lattice_parameters['values'][0]
         b = self.lattice_parameters['values'][1]
         c = self.lattice_parameters['values'][2]
         if self.uc_mask[3]:
            alpha = self.lattice_parameters['values'][3]
         else: alpha = 90
         if self.uc_mask[4]:
            beta = self.lattice_parameters['values'][3]
         else: beta = 90
         if self.uc_mask[5]:
            gamma = self.lattice_parameters['values'][3]
         else: gamma = 90
         self.unit_cell = uctbx.unit_cell(
            (a, b, c, alpha, beta, gamma))
      elif np.char.startswith(self.crystal_system, "O"):
         a = self.lattice_parameters['values'][0]
         b = self.lattice_parameters['values'][1]
         c = self.lattice_parameters['values'][2]
         self.unit_cell = uctbx.unit_cell((a, b, c, 90, 90, 90))
      elif np.char.startswith(self.crystal_system, "Te"):
         a = self.lattice_parameters['values'][0]
         if self.uc_mask[1]:
            b = self.lattice_parameters['values'][1]
         else: b = a
         if self.uc_mask[2]:
            c = self.lattice_parameters['values'][1]
         else: c = a
         self.unit_cell = uctbx.unit_cell(
            (a, b, c, 90, 90, 90))
      elif np.char.startswith(self.crystal_system, "RT"):
         a = self.lattice_parameters['values'][0]
         alpha = self.lattice_parameters['values'][1]
         self.unit_cell = uctbx.unit_cell(
            (a, a, a, alpha, alpha, alpha))
      elif np.char.startswith(self.crystal_system, "HT") \
         or np.char.startswith(self.crystal_system, "He"):
         a = self.lattice_parameters['values'][0]
         c = self.lattice_parameters['values'][1]
         self.unit_cell = uctbx.unit_cell(
            (a, a, c, 90, 90, 120))
      elif np.char.startswith(self.crystal_system, "C"):
         a = self.lattice_parameters['values'][0]
         self.unit_cell = uctbx.unit_cell((a, a, a, 90, 90, 90))

   def load_cif(self, fn, d_min=1.0, lammbda="CUA1"):
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
      cif_reader = iotbx.cif.reader(input_string=as_cif)
      self.structure = cif_reader.build_crystal_structures() \
            [os.path.split(fn)[1][0:7]]
      self.unit_cell = self.structure.unit_cell()
      self.crystal_system = self.structure.space_group().crystal_system()

      cif_model = cif_reader.model()
      try:
         # print cif_model[os.path.split(fn)[1][0:7]]['_chemical_name_mineral']
         self.chemical_name = \
            cif_model[os.path.split(fn)[1][0:7]]['_chemical_name_mineral']
      except KeyError:
         try:
            self.chemical_name = \
               cif_model[os.path.split(fn)[1][0:7]]['_chemical_name_systematic']
         except KeyError:
            self.chemical_name = os.path.split(fn)[1]

      for scatterer in self.structure.scatterers():
         if scatterer.scattering_type == "O-2":
            scatterer.scattering_type = "O2-"
         if scatterer.scattering_type == "Ca+2":
            scatterer.scattering_type = "Ca2+"
         if scatterer.scattering_type == "Si+4":
            scatterer.scattering_type = "Si4+"

   def assemble_lattice_params(self):
      """Sets the independent self.lattice_parameters attributes,
      according to the crystal system.

      """
      uc_params = self.unit_cell.parameters()
      if self.crystal_system == "Triclinic":
         self.uc_mask = [True, True, True, True, True, True]
      elif self.crystal_system == "Monoclinic":
         self.uc_mask = [True, True, True]
         for i in xrange(3, 6):
            if np.isclose(uc_params[i], 90):
               self.uc_mask.append(False)
            else: self.uc_mask.append(True)
      elif self.crystal_system == "Orthorhombic":
         self.uc_mask = [True, True, True, False, False, False]
      elif self.crystal_system == "Tetragonal":
         self.uc_mask = [True]
         if np.isclose(uc_params[1], uc_params[0]):
            self.uc_mask += [False, True]
         else: self.uc_mask += [True, False]
         self.uc_mask += [False, False, False]
      elif self.crystal_system == "Trigonal":
         if np.isclose(uc_params[3], uc_params[4]) and \
            np.isclose(uc_params[3], uc_params[5]):
            self.uc_mask = [True, False, False, True, False, False]
            self.crystal_system = "RTrigonal"
         else:
            self.uc_mask = [True, False, True, False, False, False]
            self.crystal_system = "HTrigonal"
      elif self.crystal_system == "Hexagonal":
         self.uc_mask = [True, False, True, False, False, False]
      elif self.crystal_system == "Cubic":
         self.uc_mask = [True, False, False, False, False, False]

      assert len(self.uc_mask) == 6
      self.lattice_parameters = np.array(
         [x for x in self.unit_cell_parameter_gen()], dtype=custom_dtype)

   def unit_cell_parameter_gen(self):
      """returns all unit cell parameters specified by the mask
      (a list of six booleans)
      """
      uc_labels = ["a", "b", "c", "alpha", "beta", "gamma"]
      uc_params = self.unit_cell.parameters()
      for i in xrange(6):
         if self.uc_mask[i]:
            yield ('uc_'+uc_labels[i],
                   uc_params[i],
                   uc_params[i]*(1-self.lattice_dev),
                   uc_params[i]*(1+self.lattice_dev)
                  )

   def set_lattice_parameters(self):
      """Returns a numpy array consisting of the lattice parameters

      """
      self.lattice_parameters = np.array(
         [x for x in self.unit_cell_parameter_gen()], dtype=custom_dtype)
      return self.lattice_parameters

   def compute_relative_intensities(self, anomalous_flag=True):
      r"""Returns squared structure factors, weighted by the multiplicity of
         each reflection.

         :param float d_min: The minimum resolution needed

         :return: *d*-spacings and :math:`m\times|F|^2` for each reflection,
             where *m* is the multiplicity and *F* is the structure factor.
         :rtype: numpy array

      """
      self.f_miller_set = \
         self.structure.build_miller_set(anomalous_flag,
                                         d_min=RietveldPhases.d_min
                                        ).sort()
      # Let's use scattering factors from the International Tables
      self.structure.scattering_type_registry(table="it1992") #,  "it1992",
         # "wk1995" "n_gaussian"\
      self.structure.set_inelastic_form_factors( \
         photon=wavelengths.characteristic(self.lambdas[0]), table="sasaki")

      f_calc = self.structure.structure_factors(d_min=RietveldPhases.d_min,
                                                anomalous_flag=anomalous_flag
                                               ).f_calc().sort()

      self.crystal_density = self.structure.crystal_density()
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
      rel_I_max_calc = np.amax(self.relative_intensities)
      if self.relative_intensities.shape[0] != 0:
         self.d_mask = np.logical_and( \
            self.relative_intensities > self.intensity_cutoff
            *rel_I_max_calc, self.d_spacings < RietveldPhases.d_max)
         self.d_spacings = self.d_spacings \
            [self.d_mask]
         self.relative_intensities = self.relative_intensities \
            [self.d_mask]

      # two_thetas = np.zeros((2,len(self.d_spacings)))
      # # factors = np.zeros((2,len(self.d_spacings)))
      # for i in xrange(0,len(self.lambdas),1):
      #    # read wavelength
      #    lambda_i = wavelengths.characteristic(self.lambdas[i]).as_angstrom()
      #    # Compute two_theta for each d-spacing
      #    two_thetas[i] = 360/math.pi*np.arcsin(lambda_i/2/self.d_spacings)

      # Assemble into a single array
      # self.two_theta_peaks = np.concatenate((two_thetas[0],two_thetas[1]))
      self.weighted_intensities = np.concatenate( \
         (self.K_alpha_factors[0]*self.relative_intensities, \
          self.K_alpha_factors[1]*self.relative_intensities))
      self.weighted_intensities.shape = (self.weighted_intensities.shape[0], 1)

      self.set_two_theta_peaks()
      self.masks = self.peak_masks()
      self.set_masked_arrays()

   def set_two_theta_peaks(self):
      self.two_theta_peaks = np.concatenate((
         self.unit_cell.two_theta(self.f_miller_set.indices(),
                                  wavelengths.characteristic(self.lambdas[0]
                                                            ).as_angstrom(),
                                  deg=True
                                 ).as_numpy_array()[self.d_mask],
         self.unit_cell.two_theta(self.f_miller_set.indices(),
                                  wavelengths.characteristic(self.lambdas[1]
                                                            ).as_angstrom(),
                                  deg=True
                                 ).as_numpy_array()[self.d_mask]
         ))
      self.two_theta_peaks.shape = (self.two_theta_peaks.shape[0], 1)

      self.tan_two_theta_peaks = np.tan(math.pi/360.0*self.two_theta_peaks)
      self.tan_two_theta_peaks.shape = (self.tan_two_theta_peaks.shape[0], 1)

   def set_masked_arrays(self):
      self.two_theta_masked = \
         np.broadcast_to(self.two_theta,
                         (len(self.two_theta_peaks), len(self.two_theta))
                        )[self.masks]
      self.two_theta_peaks_masked = \
         np.broadcast_to(self.two_theta_peaks,
                         (len(self.two_theta_peaks), len(self.two_theta))
                        )[self.masks]
      self.tan_two_theta_peaks_masked = \
         np.broadcast_to(self.tan_two_theta_peaks,
                         (len(self.two_theta_peaks), len(self.two_theta))
                        )[self.masks]
      self.tan_two_theta_peaks_sq_masked = self.tan_two_theta_peaks_masked**2

      self.LP_factors_masked = \
         np.broadcast_to(RietveldPhases.LP_intensity_scaling(),
                         (len(self.two_theta_peaks), len(self.two_theta))
                        )[self.masks]
      self.weighted_intensities_masked = \
         np.broadcast_to(self.weighted_intensities,
                         (len(self.two_theta_peaks), len(self.two_theta))
                        )[self.masks]

   def update_two_thetas(self, anomalous_flag=False):

      self.update_unit_cell()

      self.set_two_theta_peaks()
      # self.masks = self.peak_masks()
      self.set_masked_arrays()

   def set_eta_order(self, order):
      if order > RietveldPhases.max_polynom_order:
         order = RietveldPhases.max_polynom_order
      elif order < 1:
         order = 1
      self.eta_order = order
      self.eta = np.hstack((x for x in self.eta_param_gen(order)))
      return self.eta

   def eta_param_gen(self, order):
      n = 0
      while n < order:
         limit = np.power(0.001, n)
         if n == 0:
            yield np.array([('eta_'+str(n), 0.5, 0, 1)],
                           dtype=custom_dtype)
         else:
            yield np.array([('eta_'+str(n), 0.0, -limit, limit)],
                           dtype=custom_dtype)
         n += 1

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
                    RietveldPhases.two_theta_powers[:dim, :])

   def phase_profile(self):
      if self.recompute_peak_positions:
         self.update_two_thetas()
      # print "called phase_profile()", inspect.stack()[1][3]
      result = np.zeros((len(self.two_theta_peaks), len(self.two_theta)))
      eta_vals = np.broadcast_to(self.eta_polynomial(), \
                                 (len(self.two_theta_peaks),
                                  len(self.two_theta))
                                )[self.masks]

      if RietveldPhases.vertical_offset:
         two_theta_0_vals = \
            np.broadcast_to(RietveldPhases.cos_theta*
                            RietveldPhases.two_theta_0['values'],
                            (len(self.two_theta_peaks), len(self.two_theta))
                           )[self.masks]
      else:
         two_theta_0_vals = \
            np.broadcast_to(RietveldPhases.two_theta_0['values'],
                            (len(self.two_theta_peaks), len(self.two_theta))
                           )[self.masks]

      omegaUVW_squareds = \
         np.abs(self.U['values']*self.tan_two_theta_peaks_sq_masked
                +self.V['values']*self.tan_two_theta_peaks_masked
                +self.W['values'])
      two_thetabar_squared = (self.two_theta_masked-two_theta_0_vals
                              -self.two_theta_peaks_masked)**2 \
                             /omegaUVW_squareds

      result[self.masks] = self.Amplitude['values'] \
         *self.weighted_intensities_masked \
         *self.LP_factors_masked \
         *(eta_vals/(1+two_thetabar_squared) \
            +(1-eta_vals)*np.exp(-np.log(2)*two_thetabar_squared))
      # print result[self.masks]
      # print np.sum(result,axis=0)
      # for i in xrange(0,len(self.two_theta_peaks)):
      #    two_thetabar_squared = (self.two_theta[self.masks[i]] -two_theta_0
      #       -self.two_theta_peaks[i])**2 \
      #       /omegaUVW_squareds[i]
      #    eta = eta_vals[self.masks[i]]
      #    result[self.masks[i]] += Amplitude*self.weighted_intensities[i]* \
      #       self.LP_factors[self.masks[i]]* \
      #    (eta/(1 +two_thetabar_squared) +(1-eta) \
      #    *np.exp(-np.log(2)*two_thetabar_squared))
      # print "here"
      self.phase_profile_state = np.sum(result, axis=0)
      return self.phase_profile_state

   def update_global_and_phase_x(self, x, mask):
      self.global_and_phase_x['values'][mask] = x
      RietveldPhases.global_x['values'][RietveldPhases.global_x_no_bkgd_mask] \
         = self.global_and_phase_x['values'][self.global_mask_no_bkgd]
      self.phase_x['values'] = \
         self.global_and_phase_x['values'][self.phase_mask]

   def phase_profile_x(self, x, mask):
      self.update_global_and_phase_x(x, mask)
      return self.phase_profile()

   def phase_profile_grad(self, mask, epsilon=1e-6):
      num_params = np.sum(mask)
      result = np.zeros((num_params, len(RietveldPhases.two_theta)))
      epsilons = epsilon*np.identity(num_params)
      self.global_and_phase_x = np.hstack(
         (RietveldPhases.global_x[RietveldPhases.global_x_no_bkgd_mask],
          self.phase_x))
      self.prev_state = np.copy(self.phase_profile_state)

      for i, eps in enumerate(epsilons):
         # print eps
         self.update_global_and_phase_x(
            self.global_and_phase_x['values'][mask]+eps, mask)
         result[i, :] = (
            # self.phase_profile_x(self.global_and_phase_x['values'][mask]+eps,
            #   mask)
            self.phase_profile()-self.prev_state
            # -self.phase_profile_x(self.global_and_phase_x['values'][mask]-eps,
            #    mask)
            )/epsilon
         self.update_global_and_phase_x(
            self.global_and_phase_x['values'][mask]-eps, mask)

      # print np.sum(result, axis=1)
      return result

   def peak_masks(self, delta_theta=None):
      # print "called peak_masks()", inspect.stack()[1][3]
      if delta_theta is not None:
         return np.abs(
            self.two_theta-RietveldPhases.two_theta_0['values']
            -self.two_theta_peaks) < delta_theta
      return np.abs(
         self.two_theta-RietveldPhases.two_theta_0['values']
         - self.two_theta_peaks) < self.delta_theta
