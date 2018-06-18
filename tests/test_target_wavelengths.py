import src.cctbx_dep.target_wavelengths as tw
import numpy as np
import pytest

@pytest.fixture(scope="function")
def test_phase_settings():
   return {}

def test_set_wavelength_Co_0(test_phase_settings):
   tw.set_wavelength(test_phase_settings, 'Co')
   assert np.isclose(test_phase_settings["wavelengths"][0],
      tw.wavelengths_dict['CoA1'])

def test_set_wavelength_Cu_1(test_phase_settings):
   tw.set_wavelength(test_phase_settings, 'Cu', wavelength_model=1)
   assert np.isclose(test_phase_settings["wavelengths"][-1],
      tw.wavelengths_dict['CuA1'])

def test_error_wrong_target(test_phase_settings):
   assert pytest.raises(AssertionError,tw.set_wavelength,
      test_phase_settings, 'P')

def test_custom_wavelength(test_phase_settings):
   test_val = 1.2
   tw.set_wavelength(test_phase_settings, target=None,
      wavelength_model=2, custom_wavelength=test_val)
   assert np.isclose(test_phase_settings["wavelengths"][-1], test_val)