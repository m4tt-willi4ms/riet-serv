from __future__ import division, print_function, absolute_import
import numpy as np
from src.rietveld_phases import RietveldPhases as RP
import src.peak_masking as peak_masking

LN2 = np.sqrt(np.log(2))

def pseudo_voigt(x_squared, eta):
    return eta/(1+x_squared)/LN2 + (1-eta)*np.exp2(-x_squared)

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
    def __init__(self, phase_settings, phase_data):
        self.profile_name = phase_settings['profile']
        assert self.profile_name in PROFILES
        self.two_theta_peaks = phase_data['two_theta_peaks']
        self.delta_theta = phase_settings['delta_theta']
        self.masks = peak_masking.peak_masks(
            RP.two_theta, 0, self.two_theta_peaks, self.delta_theta)

    def set_two_theta_peaks_masked_arrays(self):
        two_theta_peaks = self.phase_data["two_theta_peaks"]
        two_theta = RietveldPhases.two_theta
        masks = self.masks
        dim = masks.shape
        self.two_theta_masked = peak_masking.get_masked_array(
            two_theta, dim, masks)
        self.two_theta_peaks_masked = peak_masking.get_masked_array(
            two_theta_peaks, dim, masks)

        tan_two_theta_peaks = self.phase_data["tan_two_theta_peaks"]
        self.tan_two_theta_peaks_masked = peak_masking.get_masked_array(
            tan_two_theta_peaks, dim, masks)
        self.tan_two_theta_peaks_sq_masked = self.tan_two_theta_peaks_masked**2

    def profile(phase_peak_parameters, global_peak_parameters):
        raise NotImplementedError

class PseudoVoigtProfile(ProfileBase):
    pass
