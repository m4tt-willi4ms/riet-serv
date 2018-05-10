import numpy as np

from refinement_parameters import RefinementParameters
from profiles import Profile

DEFAULT_U = ('caglioti_u', 0.00, [False], -0.1, 0.1)
DEFAULT_V = ('caglioti_v', 0.00, [False], -0.1, 0.1)
DEFAULT_W = ('caglioti_w', 0.002, [True], 0.000001, 1)
DEFAULT_SCALE = ('scale', 0.1, [True], 0, float('inf'))
DEFAULT_ETA_ORDER = 2
DEFAULT_LATTICE_DEV = 0.01
DEFAULT_PROFILE = 'PV'

class PhaseParameters(RefinementParameters):
    '''
    A class used to keep track of phase-specific parameters used in computing
    powder diffraction profiles.
    '''
    def __init__(self,
        validate_order_func=lambda x: x,
        scale=DEFAULT_SCALE,
        U=DEFAULT_U,
        V=DEFAULT_V,
        W=DEFAULT_W,
        eta_order=DEFAULT_ETA_ORDER,
        lattice_dev=DEFAULT_LATTICE_DEV,
        profile=DEFAULT_PROFILE,
        ):
        self.validate_order_func = validate_order_func
        self.scale = scale
        self.U = U
        self.V = V
        self.W = W
        self.eta_order = eta_order
        self.eta = set_eta_order(self.eta_order)
        self.lattice_dev = lattice_dev
        self.lattice_parameters = unit_cell.
        assert profile in profiles.profiles
        self.profile = profile
        # self.profile = profiles.Profile(DEFAULT_PROFILE)
        # self.profile_parameters = [x for x in profile.]

    def set_eta_order(self, order):
        self.eta_order = self.validate_order(order)
        self.eta = np.hstack((x for x in self.eta_param_gen(order)))
        return self.eta

    def eta_param_gen(self, order):
        n = 0
        while n < order:
            limit = np.power(0.001, n)
            if n == 0:
                yield np.array([('eta_'+str(n), 0.5, 0, 1)],
                                    dtype=CUSTOM_DTYPE)
            else:
                yield np.array([('eta_'+str(n), 0.0, -limit, limit)],
                                    dtype=CUSTOM_DTYPE)
            n += 1

    def param_gen(self):
        yield self.scale
        yield self.U
        yield self.V
        yield self.W
        yield self.eta
        yield self.lattice_parameters

    def set_param(self, param_name, val_tuple):
        assert len(val_tuple) == 5
        setattr(self, param_name, )