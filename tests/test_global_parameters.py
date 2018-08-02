from __future__ import division, print_function, absolute_import
import pytest
import numpy as np
import copy

from src.global_parameters import GlobalParameters
import src.refinement_parameters as rp

@pytest.fixture(scope="module")
def phase_settings():
    phase_settings = {}
    phase_settings["max_polynom_order"] = 5
    phase_settings["vertical_offset"] = False
    return phase_settings

@pytest.fixture(scope="module")
def gp(phase_settings):
    return GlobalParameters(phase_settings)

def test_set_bkgd_order(gp):
    gp.set_bkgd_order(3)
    assert len(gp.bkgd) == 3
    n=0
    for tup in gp.bkgd_param_gen():
        assert int(tup[0][-1]) == n
        assert tup[1] == 0.0
        n += 1
    assert n == 3

def test_assemble_x(gp):
    gp.assemble_x()
    for key in rp.keys:
        assert key in gp.x

def test_validate_order_func(phase_settings):
    phase_settings["max_polynom_order"] = 3
    gp = GlobalParameters(phase_settings)
    gp.set_bkgd_order(7)
    assert gp.bkgd_order == 3
    gp.set_bkgd_order(-5)
    assert gp.bkgd_order == 1
    assert pytest.raises(AssertionError, gp.set_bkgd_order, -.3)

def test_custom_two_theta_offset(gp, phase_settings):
    gp = GlobalParameters(phase_settings,
        two_theta_offset=('two_theta_offset', 0.5, [True, False], -1, 1))
    assert gp.two_theta_offset[1] == 0.5
    gp = GlobalParameters(phase_settings)

@pytest.fixture(scope="module")
def gp_assembled(gp):
    gp.assemble_x()
    return gp

@pytest.fixture(scope="module")
def two_theta_offset_mask(gp_assembled):
    return np.where(np.char.startswith(gp_assembled.x['labels'], 'two_t'))

@pytest.fixture(scope="module")
def bkgd_mask(gp_assembled):
    return np.where(np.char.startswith(gp_assembled.x['labels'], 'bkgd'))

def test_assemble_x(gp_assembled, two_theta_offset_mask, bkgd_mask):
    two_theta_val = gp_assembled.x['values'][two_theta_offset_mask]
    assert gp_assembled.two_theta_offset == two_theta_val
    assert np.sum(gp_assembled.x['uround'][bkgd_mask]) == gp_assembled.bkgd_order
    bkgd_vals = gp_assembled.x['values'][bkgd_mask]
    assert np.all(np.isclose(gp_assembled.bkgd, bkgd_vals))

def _set_new_two_theta_offset_val(gp_assembled, two_theta_offset_mask, val):
    new_x = copy.deepcopy(gp_assembled.x['values'])
    new_x[two_theta_offset_mask] = val
    gp_assembled.update_x(new_x, two_theta_offset_mask)
    return new_x

def test_update_x(gp_assembled, two_theta_offset_mask):
    new_two_theta_offset_val = 0.1
    new_x = _set_new_two_theta_offset_val(gp_assembled, two_theta_offset_mask,
        new_two_theta_offset_val)
    assert np.all(np.isclose(gp_assembled.x['values'], new_x))
    assert new_two_theta_offset_val == gp_assembled.two_theta_offset

def test_assembled_param_gen(gp_assembled, two_theta_offset_mask, bkgd_mask):
    gp_assembled.x['values'][two_theta_offset_mask] = 0.5
    gp_assembled.x['values'][bkgd_mask] = np.array([0, 0, 1000])
    for name, value in gp_assembled.param_gen():
        print(name, value)
        if name == 'two_theta_offset':
            assert value[1] == 0.0
            assert getattr(gp_assembled, name) == 0.5
        if name == 'bkgd':
            assert value[2][1] == 0.0
            assert getattr(gp_assembled, name)[2] == 1000

def test_reset_x(phase_settings, bkgd_mask, two_theta_offset_mask):
    gp = GlobalParameters(phase_settings)
    old_two_theta_offset = copy.deepcopy(gp.two_theta_offset)
    old_bkgd = copy.deepcopy(gp.bkgd)
    gp.assemble_x()
    gp.x['values'][bkgd_mask] = np.array([1., 2., 3.])
    gp.x['values'][two_theta_offset_mask] = .4
    gp.reset_x()
    assert old_two_theta_offset == gp.two_theta_offset
    assert old_bkgd == gp.bkgd
    assert not gp.x
    new_bkgd_0 = ('bkgd_0', 1.0, [False], -float(np.inf), float(np.inf))
    gp.bkgd[0] = new_bkgd_0
    for name, value in gp.param_gen():
        if name == 'bkgd':
            assert value[0] == new_bkgd_0
    gp.assemble_x()
    assert gp.x['values'][bkgd_mask][0] == 1.0

def test_get_dict(gp, gp_assembled):
    d = gp.as_dict()
    assert 'bkgd' in d.keys()
    assert 'two_theta_offset' in d.keys()
    assert d['bkgd'][0]['value'] == 0.0
    assert 'scale' not in d.keys()
    assert gp_assembled.x
    d = gp_assembled.as_dict()
    assert len(d['bkgd']) == gp_assembled.bkgd_order


@pytest.fixture(scope='module')
def gp_from_json(phase_settings):
    import json
    with open('./data/server_input/global_parameters_sample.json') as file:
        gp_dict = json.load(file)
    return GlobalParameters(phase_settings, param_dict=gp_dict)

def test_from_dict(phase_settings, gp_from_json):
    import json
    with open('./data/server_input/global_parameters_sample.json') as file:
        gp_dict = json.load(file)
    gp_dict_filtered = dict((k, v) for k, v in gp_dict.iteritems()
        if k not in rp.ignored_keys)
    rp.check_dict(gp_from_json.as_dict(), gp_dict_filtered)
    gp_from_json.assemble_x()
    rp.check_dict(gp_from_json.as_dict(), gp_dict_filtered)

def test_check_vertical_offset_from_json(phase_settings, gp_from_json):
    assert gp_from_json.phase_settings["vertical_offset"] == True