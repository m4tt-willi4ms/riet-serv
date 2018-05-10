import pytest

import src.refinement_parameters as rp

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