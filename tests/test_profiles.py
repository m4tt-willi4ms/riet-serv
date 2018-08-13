import pytest
import numpy as np
import src.cctbx_dep.preferred_orientation as po
import src.cctbx_dep.target_wavelengths as tw
from src.cctbx_dep.phase_from_cif import load_cif, compute_relative_intensities
from src.profiles import ProfileBase, PseudoVoigtProfile, PROFILE_CLASSES

@pytest.fixture(scope='module')
def two_theta():
    return np.linspace(5, 100, num=1000)

@pytest.fixture(scope='module')
def phase_settings():
    wavelengths = tw.set_wavelength('Cu', 0)
    d = {
        'cif_path': '.\\data\\cifs\\1000032.cif',
        'preferred_orientation': False,
        'profile': 'PV',
        'delta_theta': 2.0,
        'wavelengths': wavelengths,
        'd_min': wavelengths[1]/2,
        'd_max': wavelengths[0]/2/np.sin(np.pi/180*5),
        'intensity_cutoff': 0.01,
        'K_alpha_factors': tw.K_ALPHA_FACTORS,
        'max_polynom_order': 3,
        'vertical_offset': False,
    }
    return load_cif(d)

@pytest.fixture(scope='module')
def phase_data(phase_settings):
    return compute_relative_intensities(phase_settings)

def test_profile_base_init(phase_settings, phase_data, two_theta):
    prof = ProfileBase(phase_settings, phase_data, two_theta)
    assert prof.masks.shape == \
        (len(phase_data['two_theta_peaks']), len(two_theta))
    num_mask_true = np.sum(prof.masks)
    for item in (
            'two_theta_peaks',
            'tan_two_theta_peaks',
            'tan_two_theta_peaks_sq'):
        itemval = getattr(prof, item + '_masked')
        assert len(itemval) == num_mask_true

@pytest.fixture(scope='module')
def PVprofile(phase_settings, phase_data, two_theta):
    return PseudoVoigtProfile(
        phase_settings, phase_data, two_theta, eta_order=4)

def test_pseudo_voigt_profile(PVprofile):
    params = [
        'pp_U', 'pp_V', 'pp_W', 'pp_P',
        'pp_X', 'pp_Y',
        'eta_0', 'eta_1', 'eta_2', 
    ]
    for x in PVprofile:
        assert x[0] == params.pop(0)
    assert params == []
    gauss_params = np.array([1., 2., 3., 4.])
    lorentz_params = np.array([1., 2.])
    phase_params = np.hstack((gauss_params, lorentz_params))
    global_params = np.array([0.03])
    pv_profile = PVprofile.profile(phase_params, global_params)
    assert len(pv_profile) == np.sum(PVprofile.masks)

def test_set_eta_order(PVprofile):
    PVprofile.phase_settings["max_polynom_order"] = 3
    PVprofile.set_eta_order(7)
    assert PVprofile.eta_order == 3
    PVprofile.set_eta_order(-5)
    assert PVprofile.eta_order == 1
    assert pytest.raises(AssertionError, PVprofile.set_eta_order, -.3)
    assert PVprofile.eta_order == 1
    PVprofile.set_eta_order(3)
    n = 0
    for item in PVprofile:
        if item[0][:3] == 'eta':
            if n == 0:
                assert item[1] == 0.5
            else:
                assert item[1] == 0
            n += 1
    assert n == 3

def test_set_vertical_offset(PVprofile, phase_data, two_theta):
    assert PVprofile.phase_settings["vertical_offset"] == False
    assert 'cos_theta' not in PVprofile.__dict__
    PVprofile.phase_settings["vertical_offset"] = True
    PVprofile = PseudoVoigtProfile(
        PVprofile.phase_settings, phase_data, two_theta)
    assert 'cos_theta_masked' in PVprofile.__dict__
    assert np.isclose(PVprofile.cos_theta_masked[-1],
        np.cos(np.pi/360*PVprofile.two_theta_masked[-1]))


def test_eta_polynomial(PVprofile):
    assert 'eta_polynomial' in dir(PVprofile)
    PVprofile.set_eta_order(2)
    test_eta = np.array([0.5, 0.005])
    tmp = PVprofile.eta_polynomial(test_eta)
    assert np.isclose(tmp[-1], test_eta[0]+PVprofile.two_theta[-1]*test_eta[1])
