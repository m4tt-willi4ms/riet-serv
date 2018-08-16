from __future__ import division, print_function, absolute_import
import numpy as np

SQRTLN2 = np.sqrt(np.log(2))

def pseudo_voigt(x_squared, eta):
    return eta/(1+x_squared)/SQRTLN2 + (1-eta)*np.exp2(-x_squared)

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

def get_masked_array(array, shape, mask):
    return np.broadcast_to(array, shape)[mask]


class ProfileBase(object):
    def __init__(self, phase_settings, phase_data, two_theta):
        self.two_theta = two_theta
        self.phase_settings = phase_settings
        self.update_two_theta_peaks(phase_data)
        self.masks = self.peak_masks()
        self.set_two_theta_peaks_masked_arrays()

        max_polynom_order = phase_settings["max_polynom_order"]
        self.two_theta_powers = np.power(self.two_theta, np.array(
            xrange(0, max_polynom_order)).reshape(max_polynom_order, 1))

    def update_two_theta_peaks(self, phase_data):
        self.two_theta_peaks = phase_data['two_theta_peaks']
        self.tan_two_theta_peaks = \
            np.tan(np.pi/360.0*self.two_theta_peaks)
        self.cos_two_theta_inv_peaks = 1/np.cos(np.pi/360.0*self.two_theta_peaks)

    def set_two_theta_peaks_masked_arrays(self, phase_data=None):
        if phase_data is not None:
            self.update_two_theta_peaks(phase_data)
        dim = self.masks.shape
        masks = self.masks
        self.two_theta_masked = get_masked_array(
            self.two_theta, dim, masks)
        if self.phase_settings['vertical_offset']:
            self.cos_theta_masked = get_masked_array(
                np.cos(np.pi/360.0*self.two_theta), dim, masks)

        self.two_theta_peaks_masked = get_masked_array(
            self.two_theta_peaks, dim, masks)

        self.tan_two_theta_peaks_masked = get_masked_array(
            self.tan_two_theta_peaks, dim, masks)
        self.tan_two_theta_peaks_sq_masked = self.tan_two_theta_peaks_masked**2

        self.cos_two_theta_peaks_inv_masked = get_masked_array(
            self.cos_two_theta_inv_peaks, dim, masks)
        self.cos_two_theta_peaks_inv_sq_masked = \
            self.cos_two_theta_peaks_inv_masked**2

    def peak_masks(self):
        return np.abs(self.two_theta - self.two_theta_peaks) \
            < self.phase_settings['delta_theta']

    def calc_offset(self, offset):
        if self.phase_settings["vertical_offset"]:
            return -360/np.pi**self.cos_theta_masked * offset
        return offset

    def profile(self, phase_peak_parameters, global_peak_parameters):
        raise NotImplementedError

    def __iter__(self):
        raise NotImplementedError

class PseudoVoigtProfileBase(ProfileBase):

    def __iter__(self):
        yield ('pp_U', 0.00, [False], -0.1, 0.1)
        yield ('pp_V', 0.00, [False], -0.1, 0.1)
        yield ('pp_W', 0.001, [True], 0.000001, 1)
        yield ('pp_P', 0.00, [False], -0.1, 0.1)
        yield ('pp_X', 0.01, [False], 0.0001, 100)
        yield ('pp_Y', 0.0, [False], -0.1, 0.1)

    def calc_gamma_g(self, gauss_params):
        return np.abs(
            gauss_params[0]*self.tan_two_theta_peaks_sq_masked \
            + gauss_params[1]*self.tan_two_theta_peaks_masked \
            + gauss_params[2] \
            + gauss_params[3]*self.cos_two_theta_peaks_inv_sq_masked
        )

    def calc_gamma_l(self, lorentz_params):
        return np.abs(
            lorentz_params[0]*self.tan_two_theta_peaks_masked \
            + lorentz_params[1]*self.cos_two_theta_peaks_inv_masked
        )

class PseudoVoigtProfile(PseudoVoigtProfileBase):

    def __init__(self, phase_settings, phase_data, two_theta, eta_order=3):
        super(PseudoVoigtProfile, self).__init__(
            phase_settings, phase_data, two_theta)
        self.set_eta_order(eta_order)

    def __iter__(self):
        for param in super(PseudoVoigtProfile, self).__iter__():
            yield param
        limit = 0.1
        n = 0
        while n < self.eta_order:
            if n == 0:
                yield ('eta_'+str(n), 0.5, [True], 0, 1)
            else:
                yield ('eta_'+str(n), 0.0, [True], -limit, limit)
            n += 1

    def profile(self, phase_pp, global_pp):
        gamma_g = self.calc_gamma_g(phase_pp[:4])
        gamma_l = self.calc_gamma_l(phase_pp[4:6])
        eta = get_masked_array(
            self.eta_polynomial(phase_pp[6:]), self.masks.shape, self.masks)
        offset = self.calc_offset(global_pp[0])
        x_squared_g = (
            2*SQRTLN2*(self.two_theta_masked
                       - self.two_theta_peaks_masked - offset))**2/gamma_g
        x_squared_l = (
            2*(self.two_theta_masked
               - self.two_theta_peaks_masked - offset))**2/gamma_l
        return eta/(1+x_squared_l)/SQRTLN2 + (1-eta)*np.exp2(-x_squared_g)

    def eta_polynomial(self, eta):
        return np.clip(np.dot(
            eta, self.two_theta_powers[:len(eta), :]), 0.0, 1.0)

    def set_eta_order(self, order):
        max_polynom_order = self.phase_settings['max_polynom_order']
        assert isinstance(order, int)
        if order > max_polynom_order:
            self.eta_order = max_polynom_order
        elif order < 1:
            self.eta_order = 1
        else:
            self.eta_order = order

class PseudoVoigtTCHProfile(PseudoVoigtProfileBase):
    def calc_gamma_eta(self, gamma_g, gamma_l):
        gamma = np.power(
            np.power(gamma_g, 5)
            + 2.69269 * np.power(gamma_g, 4) * gamma_l
            + 2.42843 * np.power(gamma_g, 3) * np.power(gamma_l, 2)
            + 4.47163 * np.power(gamma_g, 2) * np.power(gamma_l, 3)
            + 0.07842 * gamma_g * np.power(gamma_l, 4)
            + np.power(gamma_l, 5), 0.2)
        gamma_l_ov_gamma = gamma_l/gamma
        eta = (
            1.366 * gamma_l_ov_gamma
            - 0.47719 * np.power(gamma_l_ov_gamma, 2)
            + 0.11116 * np.power(gamma_l_ov_gamma, 3))
        return gamma, eta

    def profile(self, phase_pp, global_pp):
        gamma_g = self.calc_gamma_g(phase_pp[:4])
        gamma_l = self.calc_gamma_l(phase_pp[4:6])
        gamma, eta = self.calc_gamma_eta(gamma_g, gamma_l)
        offset = self.calc_offset(global_pp[0])
        x_squared = (
            2*(self.two_theta_masked
                - self.two_theta_peaks_masked - offset))**2/gamma
        return eta/(1+x_squared)/SQRTLN2 + (1-eta)*np.exp2(-SQRTLN2*x_squared)


PROFILE_CLASSES = {
    'PV': PseudoVoigtProfile,
    'PV_TCH': PseudoVoigtTCHProfile,
}
