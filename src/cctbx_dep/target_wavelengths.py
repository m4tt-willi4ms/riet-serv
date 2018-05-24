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

def set_wavelength(phase_settings, target='Cu', wavelength_model=0,
   custom_wavelength=None):
   if isinstance(target, basestring):
      assert target in targets
      assert wavelength_model in (0, 1, 2)
      wavelengths = []
      if wavelength_model == 0:
         wavelengths = [
            wavelengths_dict[target+"A1"],
            wavelengths_dict[target+"A2"],
            ]
      elif wavelength_model == 1:
         wavelengths = [wavelengths_dict[target+"A1"]]
      elif wavelength_model == 2 and custom_wavelength is not None:
         wavelengths == [float(custom_wavelength)]

   phase_settings["wavelengths"] = wavelengths

   phase_settings["K_alpha_factors"] = K_alpha_factors

   epsilon = 0.05
   for bound, default_val in zip(("min_two_theta", "max_two_theta"),
                                 (0 + epsilon, 180 - epsilon)):
      if bound in phase_settings:
         value = phase_settings[bound]
      else:
         value = default_val
      exec(bound + " = {}".format(value))

   print len(wavelengths), wavelength_model, target

   phase_settings["d_max"] = wavelengths[0]/2/np.sin(np.pi/360*min_two_theta)
   phase_settings["d_min"] = wavelengths[-1]/2/np.sin(np.pi/360*max_two_theta)
