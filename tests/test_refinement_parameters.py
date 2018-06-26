import pytest
import json

import src.refinement_parameters as rp
from src.rietveld_phases import RietveldPhases

@pytest.fixture(scope="module")
def global_parameters():
   with open('./data/server_input/global_parameters_sample.json') as f:
      return json.load(f)

def test_not_implemented():
    t = rp.RefinementParameters()
    assert not t.x
    assert pytest.raises(NotImplementedError, t.param_gen)

def test_validate_order():
    t = rp.RefinementParameters()
    assert t.validate_order(7) == 5
    assert t.validate_order(7, max_polynom_order=11) == 7
    assert rp.validate_order(-10) == 1
    assert pytest.raises(AssertionError, t.validate_order, -.2)

def test_as_param(global_parameters):
    bkgd_list = rp.as_param(global_parameters['bkgd'])
    assert len(bkgd_list) == 3
    assert bkgd_list[0] == (
        'bkgd_0', 0.0, [False]*12, -3.40282347e+38, 3.40282347e+38)
    two_theta_offset = rp.as_param(global_parameters['two_theta_offset'])
    print(two_theta_offset)
    assert two_theta_offset == ('two_theta_offset', 0.0, [False]*12, -0.1, 0.1)

def test_get_dict_from_param():
    two_theta_offset = ('two_theta_offset', 0.0, [False]*12, -0.1, 0.1)
    param_dict = rp.get_dict_from_param(two_theta_offset)
    for i, key in enumerate(rp.param_keys):
        assert param_dict[key] == two_theta_offset[i]

def test_as_dict(global_parameters):
    bkgd_list = rp.as_param(global_parameters['bkgd'])
    bkgd_list_as_dict = rp.as_dict(bkgd_list)
    for param, param_as_dict in zip(bkgd_list, bkgd_list_as_dict):
        for i, key in enumerate(rp.param_keys):
            assert param_as_dict[key] == param[i]

def test_rietveld_phase_as_dict_with_cubic_crystal_system():
    cub_phase = RietveldPhases('./data/cifs/Cement/1000039-AluminateCubic.cif')
    assert len(cub_phase.phase_parameters.lattice_parameters) == 1
    phase_dict = cub_phase.as_dict()
    assert len(phase_dict['lattice_parameters']) == 6
