import numpy as np

def pseudo_voigt(x_squared, eta):
    return (eta/(1+x_squared) \
            + (1-eta)*np.exp2(-x_squared))

def gaussian(x_squared, eta=None):
    return np.exp2(-x_squared)

def lorentz(x_squared,eta=None):
    return 1/(1+x_squared)

PROFILES = {
    'PV': pseudo_voigt,
    'Lorentz': lorentz,
    'Gaussian': gaussian,
}

class ProfileBase:
    def profile_param_gen():
        raise NotImplementedError

class PseudoVoigtProfile(ProfileBase):
    pass


class Profile:
    def __init__(profile_name):
        assert profile_name in PROFILES
        self.profile_name = profile_name