import src.cctbx_dep.target_wavelengths as tw
import numpy as np
import pytest

def test_set_wavelength():
   test_phase_settings = {}

   tw.set_wavelength(test_phase_settings, 'Co')
   assert np.isclose(test_phase_settings["wavelengths"][0],
      tw.wavelengths_dict['CoA1'])

   tw.set_wavelength(test_phase_settings, 'Cu', wavelength_mode=1)
   assert np.isclose(test_phase_settings["wavelengths"][-1],
      tw.wavelengths_dict['CuA1'])

   assert pytest.raises(AssertionError,tw.set_wavelength,
      test_phase_settings, 'P')

   test_val = 1.2
   tw.set_wavelength(test_phase_settings, test_val)
   assert np.isclose(test_phase_settings["wavelengths"][-1], test_val)