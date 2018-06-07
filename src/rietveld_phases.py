from __future__ import division
import os
import time
import inspect
import sys
import numpy as np

import profiles
import peak_masking

import src.cctbx_dep.phase_from_cif as phase_from_cif
import src.cctbx_dep.unit_cell as unit_cell
import src.cctbx_dep.target_wavelengths as target_wavelengths
from src.global_parameters import GlobalParameters
from src.phase_parameters import PhaseParameters

DEFAULT_MAX_POLYNOM_ORDER = 5
DEFAULT_VERTICAL_OFFSET = False #:False = angular offset; True = Vertical Offset

DEFAULT_DELTA_THETA = 0.5
DEFAULT_INTENSITY_CUTOFF = 0.01
DEFAULT_RECOMPUTE_PEAK_POSITIONS = False
DEFAULT_LATTICE_DEV = [0.01]*6

DEFAULT_TWO_THETAS = np.linspace(0.05, 100, num=1000)

class RietveldPhases:
    r"""
        Used to collect parameters and methods for calculating powder profiles
        for Rietveld phases.

        Parameters
        -----------
        file_path_cif : string
            This string stores the location of the CIF card (including the .cif
            xtension) for the corresponding phase, either as an absolute path
            or a relative one (relative to the root directory)

        I_max : float, optional
            Can be used to specify the maximum intensity relative to which the
            computed intensities should be Scaled. If unspecified, the maximum
            intensity is determined from profile data (which can be loaded
            using the :func:`~src.RietveldPhases.RietveldPhases.set_profile()`
            class method, described later).

        delta_theta : float, optional
            :math:`\Delta\theta` specifies the region around which each peak
            profile is generated, using the formula

            .. math:: |\theta-\theta_{\rm peak}| < \Delta \theta \,.

            The default value is

        intensity_cutoff : float, optional
            The relative intensity, below which peaks are not generated. (In
            practice this is implemented when computing the squares of structure
            factors, and before any Lorentz, polarization rescaling is applied.)
            The default value of :math:`|F|^2_{\rm cutoff}` is

        lattice_dev : float, optional
            This parameter specifices the maximum allowed relative deviation of
            any lattice parameters (assuming lattice parameters will be refined).
            The default value of `lattice_dev` is 0.01.

        recompute_peak_positions : bool, optional
            This boolean variable determines whether or not peak positions are
            recomputed before determining the phase profile. (This is only
            necessary if lattice parameters are being refined.) The default value
            is

        Attributes
        ----------
        Scale : np.array (custom dtype)
            The initial input parameters for the phase Scale factor (Scale).
            Its default label, value, lower- and upper-limit are set to be

            respectively.

        U : np.array (custom dtype)
            The initial input parameters for the Caglioti `U` parameter.
            Its default label, value, lower- and upper-limit are set to be

            respectively.

        V : np.array (custom dtype)
            The initial input parameters for the Caglioti `V` parameter.
            Its default label, value, lower- and upper-limit are set to be

            respectively.

        W : np.array (custom dtype)
            The initial input parameters for the Caglioti `W` parameter.
            Its default label, value, lower- and upper-limit are set to be

            respectively.

        eta_order : int
            `eta_order` is used to specify the number of parameters to appear in
            the corresponding eta polynomial (for more information, see
            :func:`~src.RietveldPhases.RietveldPhases.eta_polynomial()`).
            The default value is

        Note
        ----
        This class contains many classmethods, which can be used by any of the
        RietveldPhases instances.

    """
    phase_settings = {}
    phase_settings["max_polynom_order"] = DEFAULT_MAX_POLYNOM_ORDER
    '''The maximum number of parameters allowed in any parameter represented
        as a polynomial (e.g. bkgd, eta)'''
    phase_settings["vertical_offset"] = DEFAULT_VERTICAL_OFFSET

    global_parameters = GlobalParameters(phase_settings)

    target_wavelengths.set_wavelength(phase_settings)

    phase_data = {}

    two_theta = None
    I = None
    sigma = None

    @classmethod
    def set_profile(cls, filename,
                         number_of_columns=3,
                         min_two_theta=0,
                         max_two_theta=180,
                         target='Cu',
                         wavelength_model=0,
                         lines_to_strip_at_TOF=1,
                        ):
        two_theta = []
        I = []
        sigma = []
        with open(filename) as file:
            for _ in xrange(lines_to_strip_at_TOF):
                file.readline()
            n = len(file.readline().split())
        if n == 2 or n == 3:
            number_of_columns = n


        with open(filename) as file:
            for line in file.readlines()[lines_to_strip_at_TOF:]:
                if number_of_columns == 2:
                    two_thetatmp, Itmp = line.split()
                # if float(two_thetatmp) < 15:
                    I.append(float(Itmp))
                    sigma.append(np.sqrt(float(Itmp)))
                elif number_of_columns == 3:
                    two_thetatmp, Itmp, sigmatmp = line.split()
                    # I.append(float(sigmatmp)**2)
                    I.append(float(Itmp))
                    sigma.append(float(sigmatmp))
                two_theta.append(float(two_thetatmp))
        cls.two_theta = np.array(two_theta)
        cls.I = np.array(I)
        cls.sigma = np.array(sigma)

        min_max_mask = np.logical_and(cls.two_theta >= min_two_theta,
                                                cls.two_theta <= max_two_theta)
        cls.I = cls.I[min_max_mask]
        cls.sigma = cls.sigma[min_max_mask]
        cls.two_theta = cls.two_theta[min_max_mask]

        cls.set_two_theta_powers_and_limits()
        cls.set_wavelength(target, wavelength_model)

    @classmethod
    def set_wavelength(cls, target, wavelength_model, custom_wavelength=None):
        target_wavelengths.set_wavelength(cls.phase_settings,
            target=target, wavelength_model=wavelength_model,
            custom_wavelength=custom_wavelength)

    @classmethod
    def set_two_theta_powers_and_limits(cls):
        cls.phase_settings["min_two_theta"] = cls.two_theta[0]
        cls.phase_settings["max_two_theta"] = cls.two_theta[-1]

        max_polynom_order = cls.phase_settings["max_polynom_order"]
        cls.two_theta_powers = np.power(cls.two_theta, np.array(
            xrange(0, max_polynom_order)).reshape(max_polynom_order, 1))
        #TODO: set cls.
        if cls.phase_settings["vertical_offset"]:
            cls.cos_theta = -360/np.pi*np.cos(np.pi/360*cls.two_theta)

    @classmethod
    def get_plot_data(cls, intensities, two_thetas=[], errors=[]):
        d = {}
        d['two_thetas'] = list(two_thetas)
        d['errors'] = list(errors)
        d['intensities'] = list(intensities)
        return d

    @classmethod
    def get_rietveld_plot(cls, profile, compute_differences=False):
        result = {}
        result['input_data'] = cls.get_plot_data(cls.I,
            two_thetas=cls.two_theta,
            errors=cls.sigma)
        result['profile_data'] = cls.get_plot_data(profile)
        if compute_differences:
            result['differences'] = cls.get_plot_data(cls.I - profile)
        else:
            result['differences'] = cls.get_plot_data([])
        return result

    @classmethod
    def background_polynomial(cls):
        r""" Returns a numpy array populated by the values of a background
        polynomial, :math:`P(2\theta)`, with input parameters :math:`c_i` stored
        in the class variable ``RietveldPhases.x`` with label ``bkgd_i``:

        .. math:: P(2\theta) = \sum_{i=0}^{N} c_i (2\theta)^i

        where *N* is the length of the numpy array ``RietveldPhases.Bkgd``.

        :return: the values of the background polynomial at points in
            ``two_theta``
        :rtype: np.array

        """
        dim = len(cls.bkgd)
        return np.dot(cls.bkgd, cls.two_theta_powers[:dim, :])

    @classmethod
    def LP_intensity_scaling(self):
        r"""
            Computes the Lorentz-Polarization intensity scaling factors for a
            set of two-theta values listed in ``two_theta``, via the equation

            .. math:: LP(2\theta) = \frac{1+\cos^2(2\theta)}{\sin\theta
                \,\sin(2\theta)} \,.
                :label: LPDefn

            :param two_theta: list of :math:`2\theta` positions

        """
        two_theta = RietveldPhases.two_theta \
            # - RietveldPhases.global_parameters.two_theta_0[:]
        return (1+np.cos(np.pi/180*two_theta)**2) \
            /np.sin(np.pi/360*two_theta) \
            /np.sin(np.pi/180*two_theta)

    @classmethod
    def assemble_global_x(cls):
        cls.global_parameters.assemble_x()
        cls.global_x = cls.global_parameters.x['values']
        cls.two_theta_0 = cls.global_parameters.two_theta_0[:]
        cls.bkgd = cls.global_parameters.bkgd[:]
        cls.global_x_no_bkgd_mask = np.invert(
            np.char.startswith(cls.global_parameters.x['labels'], 'bkgd'))
        cls.phase_data["LP_factors"] = cls.LP_intensity_scaling()

    def __init__(self, file_path_cif,
                 I_max=None,
                 delta_theta=DEFAULT_DELTA_THETA,
                 intensity_cutoff=DEFAULT_INTENSITY_CUTOFF,
                 lattice_dev=DEFAULT_LATTICE_DEV,
                 recompute_peak_positions=DEFAULT_RECOMPUTE_PEAK_POSITIONS,
                 target='Cu',
                 profile='PV',
                 phase_parameter_dict=None,
                ):

        self.phase_settings = RietveldPhases.phase_settings.copy()
        self.phase_data = RietveldPhases.phase_data.copy()
        self.file_path_cif = file_path_cif

        if I_max is not None:
            self.I_max = I_max
        elif RietveldPhases.I is not None:
            self.I_max = np.amax(RietveldPhases.I)
        else:
            self.I_max = 100

        assert intensity_cutoff <= 1 and intensity_cutoff >= 0
        self.phase_settings["intensity_cutoff"] = intensity_cutoff

        assert len(lattice_dev) == 6
        assert all([dev > 0 and dev < 1 for dev in lattice_dev])
        self.phase_settings["lattice_dev"] = lattice_dev

        assert isinstance(recompute_peak_positions, bool)
        self.phase_settings["recompute_peak_positions"] = \
            recompute_peak_positions

        assert profile in profiles.PROFILES
        self.phase_settings["profile"] = profile
        self.profile = profiles.PROFILES[profile]

        assert delta_theta > 0
        self.phase_settings["delta_theta"] = delta_theta

        self.phase_settings["cif_path"] = file_path_cif
        phase_from_cif.load_cif(self.phase_settings)

        self.phase_parameters = PhaseParameters(self.phase_settings,
            phase_parameter_dict=phase_parameter_dict)

        if RietveldPhases.two_theta is None:
            RietveldPhases.two_theta = DEFAULT_TWO_THETAS
            RietveldPhases.set_two_theta_powers_and_limits()

        RietveldPhases.assemble_global_x()
        self.assemble_phase_x()

        self.phase_data.update(phase_from_cif.compute_relative_intensities(
            self.phase_settings))

        # self.phase_data["masks"] = peak_masking.peak_masks(
        self.masks = peak_masking.peak_masks(
            RietveldPhases.two_theta,
            RietveldPhases.two_theta_0,
            self.phase_data["two_theta_peaks"],
            self.phase_settings["delta_theta"])

        # peak_masking.set_masked_arrays(self.phase_data,
        #     RietveldPhases.two_theta)
        self.set_masked_arrays()

        scale = self.phase_parameters.scale*self.I_max/ \
            np.amax(self.phase_profile())
        scale_mask = np.char.startswith(self.phase_parameters.x['labels'],'sca')
        self.phase_parameters.update_x(scale, scale_mask,
            apply_mask_to_input=False)
        # print self.global_and_phase_x

        # self.phase_parameters.scale[:] = self.phase_parameters.scale*self.I_max/ \
        #     np.amax(self.phase_profile())

    def assemble_phase_x(self):
        self.phase_parameters.assemble_x()
        self.phase_x = self.phase_parameters.x['values']

        global_x_no_bkgd = RietveldPhases.global_x \
            [RietveldPhases.global_x_no_bkgd_mask]
        # self.global_and_phase_x = np.hstack((global_x_no_bkgd,
        #     self.phase_x))

        self.global_mask_no_bkgd = np.hstack((
            np.ones(len(global_x_no_bkgd), dtype=bool),
            np.zeros(len(self.phase_x), dtype=bool)))
        self.phase_mask = np.hstack((
            np.zeros(len(global_x_no_bkgd), dtype=bool),
            np.ones(len(self.phase_x), dtype=bool)))

    def update_params(self, phase_x, mask=None):
        if mask is None:
            self.phase_x = phase_x
        else:
            self.phase_x[mask] = phase_x[mask]
        if self.phase_settings["recompute_peak_positions"]:
            #TODO: only update if lattice parameter is updated
            unit_cell.update_unit_cell(self.phase_settings,
                self.phase_parameters.lattice_parameters)

    def set_masked_arrays(self):
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

        weighted_intensities = self.phase_data["weighted_intensities"]
        self.weighted_intensities_masked = peak_masking.get_masked_array(
            weighted_intensities, dim, masks)
        self.LP_factors_masked = peak_masking.get_masked_array(
            RietveldPhases.LP_intensity_scaling(), dim, masks)
        self.update_param_arrays()

    def update_param_arrays(self):
        masks = self.masks
        dim = masks.shape
        self.eta_masked = peak_masking.get_masked_array(
            self.eta_polynomial(), dim, masks)
        if RietveldPhases.phase_settings["vertical_offset"]:
            vals = -360/np.pi*np.cos(np.pi/360*RietveldPhases.two_theta) \
                * RietveldPhases.two_theta_0
            self.two_theta_0_masked = peak_masking.get_masked_array(
                vals, dim, masks)
        else:
            self.two_theta_0_masked = RietveldPhases.two_theta_0

    def update_two_thetas(self, anomalous_flag=False):
        unit_cell.update_unit_cell(self.phase_settings,
            self.phase_parameters.lattice_parameters)
        phase_from_cif.set_two_theta_peaks(self.phase_settings,
            self.phase_data)
        # peak_masking.set_masked_arrays(self.phase_settings,
        #     RietveldPhases.two_theta)
        self.set_masked_arrays()

    def eta_polynomial(self):
        r""" Returns a numpy array populated by the values of the eta
        polynomial, :math:`\eta(2\theta)`, with input parameters :math:`\eta_i`
        stored in the class variable ``RietveldPhases.eta``:

        .. math:: P(2\theta) = \sum_{i=0}^{M} \eta_i (2\theta)^i

        where *M* is the length of the numpy array ``RietveldPhases.eta``.

        :param np.array two_theta: a list of :math:`(2\theta)` values
        :return: the values of the eta polynomial at points in
            ``two_theta``
        :rtype: np.array

        """
        eta = self.phase_parameters.eta
        return np.dot(eta, RietveldPhases.two_theta_powers[:len(eta), :])

    def phase_profile(self):
        self.update_param_arrays()
        if self.phase_settings["recompute_peak_positions"]:
            self.update_two_thetas()
        # print "called phase_profile()", inspect.stack()[1][3]
        result = np.zeros(self.masks.shape)

        omegaUVW_squareds = np.abs(
            self.phase_parameters.caglioti_u*self.tan_two_theta_peaks_sq_masked
            +self.phase_parameters.caglioti_v*self.tan_two_theta_peaks_masked
            +self.phase_parameters.caglioti_w)
        two_theta_all_squared = (self.two_theta_masked - self.two_theta_0_masked
                                        - self.two_theta_peaks_masked)**2
        two_thetabar_squared = two_theta_all_squared/omegaUVW_squareds

        # tmp = self.phase_parameters.scale \
        #     *self.weighted_intensities_masked \
        #     *self.LP_factors_masked \
        #     *self.profile(two_thetabar_squared, self.eta_masked)
        # np.putmask(result, self.masks, tmp)
        result[self.masks] = self.phase_parameters.scale \
            *self.weighted_intensities_masked \
            *self.LP_factors_masked \
            *self.profile(two_thetabar_squared, self.eta_masked)

        self.phase_profile_state = np.sum(result, axis=0)
        return self.phase_profile_state

    def increment_global_and_phase_x(self, eps, mask):
        # self.global_and_phase_x[np.where(mask)] = x

        global_x_no_bkgd = RietveldPhases.global_x \
            [RietveldPhases.global_x_no_bkgd_mask]
        global_and_phase_x = np.hstack((global_x_no_bkgd, self.phase_x))

        global_and_phase_x[mask] += eps
        # print '____RP_UPDATE_GLOBAL_AND_PHASE_X_____'
        # print 'x', x
        # print 'mask', mask
        # print 'mask, phase_mask', mask[self.phase_mask]
        # print 'phase_x', self.phase_x
        RietveldPhases.global_parameters.update_x(
            global_and_phase_x[self.global_mask_no_bkgd],
            RietveldPhases.global_x_no_bkgd_mask, apply_mask_to_input=False)
        self.phase_parameters.update_x(
            global_and_phase_x[self.phase_mask],
            mask[self.phase_mask])
        # RietveldPhases.global_x[RietveldPhases.global_x_no_bkgd_mask] \
        #     = self.global_and_phase_x[self.global_mask_no_bkgd]
        # self.phase_x = \
        #     self.global_and_phase_x[self.phase_mask]

    def phase_profile_x(self, x, mask):
        self.update_global_and_phase_x(x, mask)
        return self.phase_profile()

    def phase_profile_grad(self, global_and_phase_mask, epsilon=1e-6):
        num_params = np.sum(global_and_phase_mask)
        result = np.zeros((num_params, len(RietveldPhases.two_theta)))
        epsilons = epsilon*np.identity(num_params)

        # global_x_no_bkgd = RietveldPhases.global_x \
        #     [RietveldPhases.global_x_no_bkgd_mask]
        # global_and_phase_x = np.hstack((global_x_no_bkgd, self.phase_x))

        prev_state = np.copy(self.phase_profile_state)

        for i, eps in enumerate(epsilons):
            # print eps
            # print 'Before:', self.phase_parameters.x['values']
            self.increment_global_and_phase_x(eps, global_and_phase_mask)
            # print 'After:', self.phase_parameters.x['values']
            result[i, :] = (
                # self.phase_profile_x(self.global_and_phase_x['values'][global_and_phase_mask]+eps,
                #   global_and_phase_mask)
                self.phase_profile()-prev_state
                # -self.phase_profile_x(self.global_and_phase_x['values'][global_and_phase_mask]-eps,
                #    global_and_phase_mask)
                )/epsilon
            self.increment_global_and_phase_x(-eps, global_and_phase_mask)

        # print np.sum(result, axis=1)
        return result

    def update_phase_info(self, phase_dict):
        phase_dict['cif_path'] = self.file_path_cif
        phase_dict['phase_name'] = self.phase_settings['chemical_name']
        lattice_parameters = []
        for param in unit_cell.unit_cell_parameter_gen(self.phase_settings,
            np.ones(6, dtype=bool)):
            d = {}
            d['name'] = param[0]
            d['value'] = param[1]
            d['l_limit'] = param[3]
            d['u_limit'] = param[4]
            d['used'] = bool(np.any(np.array(param[2], dtype=bool)))
            d['uround'] = [bool(x) for x in np.nditer(param[2])]
            d['round'] = 2
            lattice_parameters.append(d)
        phase_dict['lattice_parameters'] = lattice_parameters
        phase_dict['lattice_parameter_tolerances'] = \
            self.phase_settings['lattice_dev']

    def as_dict(self):
        d = self.phase_parameters.as_dict()
        self.update_phase_info(d)
        return d

if __name__ == '__main__':
    RietveldPhases.set_profile('./data/profiles/d5_05005.xye')
    phase = RietveldPhases('./data/cifs/9015662-rutile.cif')
    profile = phase.phase_profile() + RietveldPhases.background_polynomial()
    profile_grad = phase.phase_profile_grad(
        np.ones(len(phase.phase_x)+1, dtype=bool))
