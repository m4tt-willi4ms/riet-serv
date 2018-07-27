import pytest
import numpy as np
import src.cctbx_dep.preferred_orientation as po
import src.cctbx_dep.target_wavelengths as tw
from src.cctbx_dep.phase_from_cif import load_cif, compute_relative_intensities

@pytest.fixture(scope='module')
def phase_settings():
    wavelengths = tw.set_wavelength('Cu', 0)
    d = {
        'preferred_orientation': True,
        'pref_orient_method': 'march_dollase',
        'pref_orient_hkl': [1, 0, 0],
        'cif_path': '.\\data\\cifs\\Gypsum.cif',
        'wavelengths': wavelengths,
        'd_min': wavelengths[1]/2,
        'd_max': wavelengths[0]/2/np.sin(np.pi/180*5),
        'intensity_cutoff': 0.01,
        'K_alpha_factors': tw.K_ALPHA_FACTORS
    }
    return load_cif(d)

@pytest.fixture(scope='module')
def phase_data(phase_settings):
    return compute_relative_intensities(phase_settings)

@pytest.fixture(scope='module')
def test_params():
    return [('md_param', 1.00, [True], 0.000001, 2)]

def test_get_pref_orient_params(phase_settings, test_params):
    params = po.get_pref_orient_params(phase_settings, test_params)
    assert test_params == params

def test_get_sym_equiv_angles(phase_settings, phase_data):
    angles = po.get_sym_equiv_angles(phase_settings, phase_data)
    assert np.all(np.isclose(
        phase_data['f_miller_indices'][0,:], np.array([0, 2, 0])))
    assert np.all(np.isclose(angles[0], [90, 90]))
    multiplicities = phase_data['multiplicities']
    for i, angle in enumerate(angles):
        assert len(angle) == multiplicities[i]

def test_pref_orient_function():
    val1 = po.pref_orient_function(1., [0.])
    assert np.isclose(val1, 1.0)
    val2 = po.pref_orient_function(1.5, [0.])
    assert np.isclose(val2, np.power(1.5**2, -1.5))
    val3 = po.pref_orient_function(.7, np.array([20., 30.]))
    assert np.isclose(val3, (
        np.power(1/.7+(.7**2-1/.7)*np.cos(np.pi/180*20.)**2, -1.5)
        + np.power(1/.7+(.7**2-1/.7)*np.cos(np.pi/180*30.)**2, -1.5))/2)

def test_update_pref_orient_factors(phase_settings, phase_data):
    r1 = [1.2]
    comp = po.update_pref_orient_factors(phase_settings, phase_data, r1)
    for index, sea in enumerate(phase_data['sym_equiv_angles']):
        assert np.isclose(comp[index], po.pref_orient_function(r1[0], sea))


# def test_get_pref_orient_function(phase_settings, phase_data):
#     result = pr.get_pref_orient_function(
#         phase_settings, phase_data)
#     assert 0
