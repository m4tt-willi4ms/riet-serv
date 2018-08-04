from __future__ import division, print_function, absolute_import
import numpy as np

LN2 = np.sqrt(np.log(2))

def pseudo_voigt(x_squared, eta):
    return eta/(1+x_squared)/LN2 + (1-eta)*np.exp2(-x_squared)

def gaussian(x_squared, eta=None):
    return np.exp2(-x_squared)

def lorentz(x_squared,eta=None):
    return 1/(1+x_squared)

def pseudo_voigtTCH(two_theta, gamma_l, gamma_g):
    pass

PROFILES = {
    'PV': pseudo_voigt,
    'PV-TCH': pseudo_voigtTCH,
    'Lorentz': lorentz,
    'Gaussian': gaussian,
}

class ProfileBase:
    def profile_param_gen():
        raise NotImplementedError

class PseudoVoigtProfile(ProfileBase):
    pass


class Profile:
    def __init__(profile_name, two_theta):
        assert profile_name in PROFILES
        self.profile_name = profile_name
        self.two_theta = two_theta
        