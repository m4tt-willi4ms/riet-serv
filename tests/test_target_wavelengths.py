import src.cctbx_dep.target_wavelengths as tw
import numpy as np
import pytest

def test_set_wavelength_Co_0():
   wavelengths = tw.set_wavelength('Co')
   assert np.isclose(wavelengths[0], tw.WAVELENGTHS_DICT['CoA1'])

def test_set_wavelength_Cu_1():
   wavelengths = tw.set_wavelength('Cu', wavelength_model=1)
   assert np.isclose(wavelengths, tw.WAVELENGTHS_DICT['CuA1'])

def test_error_wrong_target():
   assert pytest.raises(AssertionError,tw.set_wavelength, 'P')

def test_custom_wavelength():
   test_val = 1.2
   wavelengths = tw.set_wavelength(target=None,
      wavelength_model=2, custom_wavelength=test_val)
   assert np.isclose(wavelengths, test_val)