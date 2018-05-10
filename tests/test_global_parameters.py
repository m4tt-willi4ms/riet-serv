import pytest
import numpy as np
import copy

from src.global_parameters import GlobalParameters
import src.refinement_parameters as rp

@pytest.fixture(scope="module")
def gp():
    return GlobalParameters()

def test_set_bkgd_order(gp):
    gp.set_bkgd_order(3)
    assert len(gp.bkgd) == 3
    n=0
    for tup in gp.bkgd_param_gen():
        assert int(tup[0][-1]) == n
        assert tup[1] == 0.0
        n += 1
    assert n == 3

def test_set_vertical_offset(gp):
    assert gp.vertical_offset == False
    gp.set_vertical_offset(True)
    assert gp.vertical_offset == True
    gp.set_vertical_offset(False)
    assert gp.vertical_offset == False

def test_assemble_x(gp):
    gp.assemble_x()
    for key in rp.keys:
        assert key in gp.x

def test_validate_order_func():
    def val_func(order):
        assert isinstance(order, int)
        if order > 3:
            order = 3
        elif order < 0:
            order = 0
        return order
    rp.DEFAULT_MAX_POLYNOM_ORDER = 3

    gp = GlobalParameters()
    gp.set_bkgd_order(7)
    assert gp.bkgd_order == 3
    gp.set_bkgd_order(-5)
    assert gp.bkgd_order == 0
    assert pytest.raises(AssertionError, gp.set_bkgd_order, -.3)

def test_custom_two_theta_0(gp):
    gp = GlobalParameters(
        two_theta_0=('two_theta_0', 0.5, [True, False], -1, 1))
    assert gp.two_theta_0[1] == 0.5
    gp = GlobalParameters()

@pytest.fixture(scope="module")
def gp_assembled(gp):
    gp.assemble_x()
    return gp

@pytest.fixture(scope="module")
def two_theta_0_mask(gp_assembled):
    return np.where(np.char.startswith(gp_assembled.x['labels'], 'two_t'))

@pytest.fixture(scope="module")
def bkgd_mask(gp_assembled):
    return np.where(np.char.startswith(gp_assembled.x['labels'], 'bkgd'))

def test_assemble_x(gp_assembled, two_theta_0_mask, bkgd_mask):
    two_theta_val = gp_assembled.x['values'][two_theta_0_mask]
    assert gp_assembled.two_theta_0 == two_theta_val
    assert np.sum(bkgd_mask) == gp_assembled.bkgd_order
    bkgd_vals = gp_assembled.x['values'][bkgd_mask]
    assert np.all(np.isclose(gp_assembled.bkgd, bkgd_vals))

def _set_new_two_theta_0_val(gp_assembled, two_theta_0_mask, val):
    new_x = copy.deepcopy(gp_assembled.x['values'])
    new_x[two_theta_0_mask] = val
    gp_assembled.update_x(new_x, two_theta_0_mask)
    return new_x

def test_update_x(gp_assembled, two_theta_0_mask):
    new_two_theta_0_val = 0.1
    new_x = _set_new_two_theta_0_val(gp_assembled, two_theta_0_mask,
        new_two_theta_0_val)
    assert np.all(np.isclose(gp_assembled.x['values'], new_x))
    assert new_two_theta_0_val == gp_assembled.two_theta_0

def test_assembled_param_gen(gp_assembled, two_theta_0_mask, bkgd_mask):
    gp_assembled.x['values'][two_theta_0_mask] = 0.5
    gp_assembled.x['values'][bkgd_mask] = np.array([0, 0, 1000])
    for name, value in gp_assembled.param_gen():
        print name, value
        if name == 'two_theta_0':
            assert value[1] == 0.0
            assert getattr(gp_assembled, name) == 0.5
        if name == 'bkgd':
            assert value[2][1] == 0.0
            assert getattr(gp_assembled, name)[2] == 1000

def test_reset_x(bkgd_mask, two_theta_0_mask):
    gp = GlobalParameters()
    old_two_theta_0 = copy.deepcopy(gp.two_theta_0)
    old_bkgd = copy.deepcopy(gp.bkgd)
    gp.assemble_x()
    gp.x['values'][bkgd_mask] = np.array([1., 2., 3.])
    gp.x['values'][two_theta_0_mask] = .4
    gp.reset_x()
    assert old_two_theta_0 == gp.two_theta_0
    assert old_bkgd == gp.bkgd
    assert not gp.x
    new_bkgd_0 = ('bkgd_0', 1.0, [False], -float(np.inf), float(np.inf))
    gp.bkgd[0] = new_bkgd_0
    for name, value in gp.param_gen():
        if name == 'bkgd':
            assert value[0] == new_bkgd_0
    gp.assemble_x()
    print gp.x
    assert gp.x['values'][bkgd_mask][0] == 1.0
