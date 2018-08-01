from __future__ import division, print_function, absolute_import
import pytest
import os
import numpy as np
import copy
import json

from src.phase_parameters import PhaseParameters
import src.refinement_parameters as rp
import src.cctbx_dep.phase_from_cif as phase_from_cif
import src.cctbx_dep.unit_cell as unit_cell
import src.cctbx_dep.target_wavelengths as target_wavelengths
import src.profiles as profiles

@pytest.fixture(scope="module")
def phase_settings():
    phase_settings = {
        'cif_path': './data/cifs/1000032.cif'
    }
    phase_settings['wavelengths'] = target_wavelengths.set_wavelength('Cu', 0)
    phase_settings = phase_from_cif.load_cif(phase_settings)
    unit_cell.set_unit_cell_mask(phase_settings)
    phase_settings["max_polynom_order"] = 5
    phase_settings["lattice_dev"] = [0.01]*6
    phase_settings["recompute_peak_positions"] = True
    phase_settings["preferred_orientation"] = True
    return phase_settings

@pytest.fixture
def pp(phase_settings):
    return PhaseParameters(phase_settings)

def test_set_eta_order(pp):
    pp.set_eta_order(3)
    assert len(pp.eta) == 3
    n=0
    for tup in pp.eta_param_gen():
        assert int(tup[0][-1]) == n
        if n == 0:
            assert tup[1] == 0.5
        else:
            assert tup[1] == 0.0
        n += 1
    assert n == 3

def test_set_U(pp):
    pp.reset_x()
    new_U_val = 0.01
    pp.assemble_x()
    U_mask = np.where(np.char.startswith(pp.x['labels'], 'cagliotti_u'))
    pp.update_x(new_U_val, U_mask, apply_mask_to_input=False)
    assert pp.x['values'][U_mask] == new_U_val
    assert pp.cagliotti_u == new_U_val

def test_assemble_x(pp):
    pp.assemble_x()
    for key in rp.keys:
        assert key in pp.x

def test_validate_order_func(phase_settings):
    phase_settings["max_polynom_order"] = 3
    pp = PhaseParameters(phase_settings)
    pp.set_eta_order(7)
    assert pp.eta_order == 3
    pp.set_eta_order(-5)
    assert pp.eta_order == 1
    assert pytest.raises(AssertionError, pp.set_eta_order, -.3)

def test_custom_scale(phase_settings):
    pp = PhaseParameters(phase_settings,
        scale=('blah', 10.0, [False], 0, float('inf')))
    assert pp.scale[1] == 10.0
    pp.assemble_x()
    assert pp.scale == 10.0
    assert isinstance(pp.scale, np.ndarray)
    scale_mask = np.where(np.char.startswith(pp.x['labels'], 'blah'))
    x_new = np.zeros(len(scale_mask), dtype=float)
    x_new[scale_mask] = 100.0
    pp.update_x(x_new, scale_mask)
    assert pp.scale == 100.0
    pp.reset_x()
    assert pp.scale[1] == 10.0
    assert pp.scale[0] == 'blah'

@pytest.fixture(scope="module")
def pp_assembled(phase_settings):
    pp = PhaseParameters(phase_settings)
    pp.assemble_x()
    return pp

@pytest.fixture(scope="module")
def scale_mask(pp_assembled):
    return np.where(np.char.startswith(pp_assembled.x['labels'], 'sca'))

@pytest.fixture(scope="module")
def eta_mask(pp_assembled):
    return np.where(np.char.startswith(pp_assembled.x['labels'], 'eta'))

def test_assemble_x(pp_assembled, scale_mask, eta_mask):
    scale_val = pp_assembled.x['values'][scale_mask]
    assert pp_assembled.scale == scale_val
    assert np.sum(pp_assembled.x['uround'][scale_mask]) == 1
    eta_vals = pp_assembled.x['values'][eta_mask]
    assert np.all(np.isclose(pp_assembled.eta, eta_vals))

@pytest.fixture(scope="module")
def phase_dict():
    with open(os.path.join(os.path.dirname(__file__),
        '../data/server_input/phase_parameters_sample.json')) as f:
        return json.load(f)

@pytest.fixture(scope="module")
def pp_from_json(phase_settings, phase_dict):
    pp_from_json = PhaseParameters(phase_settings, param_dict=phase_dict)
    pp_from_json.assemble_x()
    return pp_from_json

@pytest.fixture(scope="module")
def lp_mask(pp_from_json):
    return np.where(np.char.startswith(pp_from_json.x['labels'], 'uc'))

def test_lattice_parameters(pp_from_json, lp_mask, phase_dict):
    json_lps = phase_dict['lattice_parameters']
    assert list(pp_from_json.x['labels'][lp_mask]) == ['uc_a', 'uc_c']

def test_pref_or(pp_from_json):
    for key in ('pref_orient_hkl', 'pref_orient_method', 'pref_orient_ell'):
        assert key in pp_from_json.phase_settings

def test_as_dict(pp_from_json):
    phase_dict = pp_from_json.as_dict()
    assert 'scale' in phase_dict
    assert isinstance(phase_dict['pref_orient'], list)
# def _set_new_two_theta_0_val(gp_assembled, two_theta_0_mask, val):
#     new_x = copy.deepcopy(gp_assembled.x['values'])
#     new_x[two_theta_0_mask] = val
#     gp_assembled.update_x(new_x, two_theta_0_mask)
#     return new_x

# def test_update_x(gp_assembled, two_theta_0_mask):
#     new_two_theta_0_val = 0.1
#     new_x = _set_new_two_theta_0_val(gp_assembled, two_theta_0_mask,
#         new_two_theta_0_val)
#     assert np.all(np.isclose(gp_assembled.x['values'], new_x))
#     assert new_two_theta_0_val == gp_assembled.two_theta_0

# def test_assembled_param_gen(gp_assembled, two_theta_0_mask, bkgd_mask):
#     gp_assembled.x['values'][two_theta_0_mask] = 0.5
#     gp_assembled.x['values'][bkgd_mask] = np.array([0, 0, 1000])
#     for name, value in gp_assembled.param_gen():
#         print name, value
#         if name == 'two_theta_0':
#             assert value[1] == 0.0
#             assert getattr(gp_assembled, name) == 0.5
#         if name == 'bkgd':
#             assert value[2][1] == 0.0
#             assert getattr(gp_assembled, name)[2] == 1000

# def test_reset_x(phase_settings, bkgd_mask, two_theta_0_mask):
#     gp = GlobalParameters(phase_settings)
#     old_two_theta_0 = copy.deepcopy(gp.two_theta_0)
#     old_bkgd = copy.deepcopy(gp.bkgd)
#     gp.assemble_x()
#     gp.x['values'][bkgd_mask] = np.array([1., 2., 3.])
#     gp.x['values'][two_theta_0_mask] = .4
#     gp.reset_x()
#     assert old_two_theta_0 == gp.two_theta_0
#     assert old_bkgd == gp.bkgd
#     assert not gp.x
#     new_bkgd_0 = ('bkgd_0', 1.0, [False], -float(np.inf), float(np.inf))
#     gp.bkgd[0] = new_bkgd_0
#     for name, value in gp.param_gen():
#         if name == 'bkgd':
#             assert value[0] == new_bkgd_0
#     gp.assemble_x()
#     print gp.x
#     assert gp.x['values'][bkgd_mask][0] == 1.0
