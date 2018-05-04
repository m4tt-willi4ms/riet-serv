import numpy as np

from refinement_parameters import RefinementParameters
from profiles import Profile

DEFAULT_U = ('caglioti_u', 0.00, False, -0.1, 0.1)
DEFAULT_V = ('caglioti_v', 0.00, False, -0.1, 0.1)
DEFAULT_W = ('caglioti_w', 0.002, True, 0.000001, 1)
DEFAULT_SCALE = ('scale', 0.1, True, 0, float('inf'))
DEFAULT_ETA_ORDER = 2
DEFAULT_LATTICE_DEV = 0.01
DEFAULT_PROFILE = 'PV'

class PhaseParameters(RefinementParameters):
    '''
    A class used to keep track of phase-specific parameters used in computing
    powder diffraction profiles.
    '''
    def __init__(self, validate_order_func=lambda x: x):
        self.validate_order_func = validate_order_func
        self.scale = DEFAULT_SCALE
        self.U = DEFAULT_U
        self.V = DEFAULT_V
        self.W = DEFAULT_W
        # self.profile = profiles.Profile(DEFAULT_PROFILE)
        # self.profile_parameters = [x for x in profile.]
