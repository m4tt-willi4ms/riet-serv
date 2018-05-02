import numpy as np

from cctbx.eltbx import wavelengths

wavelengths_dict = {
   "CoA1": 1.788965,
   "CoA2": 1.792850,
}
for w in wavelengths.characteristic_iterator():
   wavelengths_dict[w.label()] = w.as_angstrom()

targets = ('Ag','Mo','Cu','Cr','Fe','Co')

K_alpha_factors = [1.0, 0.48]

def set_wavelength(phase_settings, target='Cu', wavelength_mode=2):
   if type(target) == str:
      assert target in targets
      if wavelength_mode == 1:
         wavelengths = [wavelengths_dict[target+"A1"]]
      if wavelength_mode == 2:
         wavelengths = [
            wavelengths_dict[target+"A1"],
            wavelengths_dict[target+"A2"],
            ]
   elif type(target) == float:
      assert target > 0
      wavelengths = [target]
   _set_wavelength_only(phase_settings, wavelengths)

   phase_settings["K_alpha_factors"] = K_alpha_factors

   epsilon = 0.05
   for bound, default_val in zip(("min_two_theta", "max_two_theta"),
                                 (0 + epsilon, 180 - epsilon)):
      if bound in phase_settings:
         exec(bound + " = {}".format(phase_settings[bound]))
      else:
         phase_settings[bound] = default_val
         exec(bound + " = {}".format(default_val))

   phase_settings["d_max"] = wavelengths[0]/2/np.sin(np.pi/360*min_two_theta)
   phase_settings["d_min"] = wavelengths[-1]/2/np.sin(np.pi/360*max_two_theta)

def _set_wavelength_only(phase_settings, wavelength_list):
   if "wavelengths" in phase_settings:
      phase_settings["wavelengths"] = []
   phase_settings["wavelengths"] = wavelength_list
