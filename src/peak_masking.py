'''This module contains methods used in computing peak masks and masking arrays
'''
import numpy as np

def peak_masks(two_theta, two_theta_0, two_theta_peaks, delta_theta):
    '''

    '''
    # print "called peak_masks()", inspect.stack()[1][3]
    return np.abs(two_theta  - two_theta_peaks) < delta_theta #- two_theta_0

def get_masked_array(array, shape, mask):
    return np.broadcast_to(array, shape)[mask]

def set_masked_arrays(phase_data, two_theta):
    two_theta_peaks = phase_data["two_theta_peaks"]
    # two_theta = RietveldPhases.two_theta
    masks = phase_data["masks"]
    # print masks.shape
    # dim = (len(two_theta_peaks), len(two_theta))
    dim = masks.shape
    self.two_theta_masked = get_masked_array(two_theta, dim, masks)
    # two_theta_masked = np.broadcast_to(two_theta,
    #         (len(two_theta_peaks), len(two_theta)))[masks]
    self.two_theta_peaks_masked = get_masked_array(two_theta_peaks, dim, masks)

    tan_two_theta_peaks = phase_data["tan_two_theta_peaks"]
    self.tan_two_theta_peaks_masked = get_masked_array(
        tan_two_theta_peaks, dim, masks)
    self.tan_two_theta_peaks_sq_masked = self.tan_two_theta_peaks_masked**2

    self.LP_factors_masked = get_masked_array(
        phase_data["LP_factors"], dim, masks)

    weighted_intensities = phase_data["weighted_intensities"]
    self.weighted_intensities_masked = get_masked_array(
        weighted_intensities, dim, masks)
    self.eta_masked = peak_masking.get_masked_array(
        self.eta_polynomial(), dim, masks)

    if RietveldPhases.phase_settings["vertical_offset"]:
        vals = -360/np.pi*np.cos(np.pi/360*RietveldPhases.two_theta) \
             * RietveldPhases.two_theta_0
        self.two_theta_0_masked = peak_masking.get_masked_array(
             vals, dim, masks)
    else:
        self.two_theta_0_masked = peak_masking.get_masked_array(
             RietveldPhases.two_theta_0, dim, masks)

    self.two_theta_all_squared = (self.two_theta_masked
                                             -self.two_theta_0_masked
                                             -self.two_theta_peaks_masked)**2
