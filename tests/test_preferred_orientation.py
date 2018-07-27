import pytest
import src.cctbx_dep.preferred_orientation as pr

@pytest.fixture(scope='module')
def phase_settings():
    return {
        'pref_orient_method': 'march_dollase',
    }

def test_get_pref_orient_params(phase_settings):
    test_params = [('md_param', 1.00, [True], 0.000001, 2)]
    params = pr.get_pref_orient_params(phase_settings, test_params)
    assert test_params == params
