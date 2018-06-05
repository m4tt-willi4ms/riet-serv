"""
   This module gathers together code used to populate the phase_settings
   dictionary using cctbx's iotbx.cif module.
"""
import os
import numpy as np

import iotbx.cif

def load_cif(phase_settings):
   """Reads in a unit cell, crystal structure, crystal system from iotbx

      :param str file_path: The file name

   """
   file_path = phase_settings["cif_path"]
   with open(file_path, 'r') as opened_file:
      as_cif = opened_file.read()
   cif_reader = iotbx.cif.reader(input_string=as_cif)
   structure = cif_reader.build_crystal_structures() \
      [os.path.split(file_path)[1][0:7]]
   phase_settings["structure"] = structure
   phase_settings["unit_cell"] = structure.unit_cell()
   phase_settings["crystal_system"] = \
      structure.space_group().crystal_system()

   cif_model = cif_reader.model()
   try:
      phase_settings["chemical_name"] = \
         cif_model[os.path.split(file_path)[1][0:7]]['_chemical_name_mineral']
   except KeyError:
      try:
         phase_settings["chemical_name"] = cif_model[os.path.split(
            file_path)[1][0:7]]['_chemical_name_systematic']
      except KeyError:
         phase_settings["chemical_name"] = os.path.split(file_path)[1]

   for scatterer in structure.scatterers():
      if scatterer.scattering_type == "O-2":
         scatterer.scattering_type = "O2-"
      if scatterer.scattering_type == "Ca+2":
         scatterer.scattering_type = "Ca2+"
      if scatterer.scattering_type == "Si+4":
         scatterer.scattering_type = "Si4+"

   return phase_settings

def compute_relative_intensities(phase_settings, anomalous_flag=True):
   r"""Returns squared structure factors, weighted by the multiplicity of
      each reflection.

      :param float d_min: The minimum resolution needed

      :return: *d*-spacings and :math:`m\times|F|^2` for each reflection,
          where *m* is the multiplicity and *F* is the structure factor.
      :rtype: numpy array

   """
   phase_data = {}

   structure = phase_settings["structure"]
   unit_cell = phase_settings["unit_cell"]

   wavelengths = phase_settings["wavelengths"]
   d_min = phase_settings["d_min"]

   f_miller_set = structure.build_miller_set(anomalous_flag, d_min=d_min).sort()
   # Let's use scattering factors from the International Tables
   structure.scattering_type_registry(table="it1992") #,  "it1992",
      # "wk1995" "n_gaussian"\
   structure.set_inelastic_form_factors(wavelengths[0], table="sasaki")

   f_calc = structure.structure_factors(d_min=d_min,
                                        anomalous_flag=anomalous_flag
                                       ).f_calc().sort()

   phase_data["crystal_density"] = structure.crystal_density()
   unit_cell_volume = unit_cell.volume()

   f_calc_sq = f_calc.as_intensity_array().sort().data() \
      /unit_cell_volume/unit_cell_volume
   f_calc_mult = f_calc.multiplicities().sort().data()

   d_spacings = f_miller_set.d_spacings().data().as_numpy_array()
   relative_intensities = f_calc_sq * f_calc_mult.as_double().as_numpy_array()
   #: weight intensities by the corresponding multiplicity
   # print "Volume: " + str(self.unit_cell.volume())

   intensity_cutoff = phase_settings["intensity_cutoff"]
   d_max = phase_settings["d_max"]

   # Drop any peaks below the Intensity Cutoff
   rel_I_max_calc = np.amax(relative_intensities)
   if relative_intensities.shape[0] != 0:
      d_mask = np.logical_and(relative_intensities > intensity_cutoff \
         *rel_I_max_calc, d_spacings < d_max)
      d_spacings = d_spacings[d_mask]
      relative_intensities = relative_intensities[d_mask]

   phase_data["f_miller_set"] = f_miller_set
   phase_data["d_spacings"] = d_spacings
   phase_data["relative_intensities"] = relative_intensities
   phase_data["d_mask"] = d_mask

   # two_thetas = np.zeros((2,len(self.d_spacings)))
   # # factors = np.zeros((2,len(self.d_spacings)))
   # for i in xrange(0,len(self.lambdas),1):
   #    # read wavelength
   #    lambda_i = wavelengths.characteristic(self.lambdas[i]).as_angstrom()
   #    # Compute two_theta for each d-spacing
   #    two_thetas[i] = 360/math.pi*np.arcsin(lambda_i/2/self.d_spacings)

   K_alpha_factors = phase_settings["K_alpha_factors"]
   # Assemble into a single array
   # self.two_theta_peaks = np.concatenate((two_thetas[0],two_thetas[1]))
   weighted_intensities = np.concatenate(
      (K_alpha_factors[0]*relative_intensities, \
       K_alpha_factors[1]*relative_intensities))
   weighted_intensities.shape = (weighted_intensities.shape[0], 1)

   phase_data["weighted_intensities"] = weighted_intensities

   set_two_theta_peaks(phase_settings, phase_data)
   # self.masks = self.peak_masks()
   # self.set_masked_arrays()

   return phase_data

def set_two_theta_peaks(phase_settings, phase_data):
   unit_cell = phase_settings["unit_cell"]
   wavelengths = phase_settings["wavelengths"]
   f_miller_set = phase_data["f_miller_set"]
   d_mask = phase_data["d_mask"]
   two_theta_peaks = np.concatenate((
      unit_cell.two_theta(f_miller_set.indices(), wavelengths[0], deg=True
                         ).as_numpy_array()[d_mask],
      unit_cell.two_theta(f_miller_set.indices(), wavelengths[1], deg=True,
                         ).as_numpy_array()[d_mask]
      ))
   two_theta_peaks.shape = (two_theta_peaks.shape[0], 1)

   tan_two_theta_peaks = np.tan(np.pi/360.0*two_theta_peaks)
   tan_two_theta_peaks.shape = (tan_two_theta_peaks.shape[0], 1)
   #TODO: use [:,np.newaxis] instead

   phase_data["two_theta_peaks"] = two_theta_peaks
   phase_data["tan_two_theta_peaks"] = tan_two_theta_peaks
