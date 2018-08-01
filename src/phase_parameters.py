from __future__ import division, print_function, absolute_import
from collections import OrderedDict
import numpy as np

import src.refinement_parameters as rp
import src.profiles as profiles
import src.cctbx_dep.unit_cell as unit_cell
import src.cctbx_dep.preferred_orientation as po

DEFAULT_U = ('cagliotti_u', 0.00, [False], -0.1, 0.1)
DEFAULT_V = ('cagliotti_v', 0.00, [False], -0.1, 0.1)
DEFAULT_W = ('cagliotti_w', 0.001, [True], 0.000001, 1)
DEFAULT_SCALE = ('scale', 0.1, [True], 0, float('inf'))
DEFAULT_PREF_OR_PARAMS = [('pref_or_0', 1.00, [True], 0.000001, 2)]
DEFAULT_PREF_OR_HKL = [0, 0, 1]
DEFAULT_PREF_OR_METHOD = 'march_dollase'
DEFAULT_PREF_OR_ELL = 4
DEFAULT_ETA_ORDER = 2
DEFAULT_LATTICE_DEV = 0.01
DEFAULT_PROFILE = 'PV'

class PhaseParameters(rp.RefinementParameters):
    '''
    A class used to keep track of phase-specific parameters used in computing
    powder diffraction profiles.
    '''
    def __init__(self, phase_settings,
        scale=DEFAULT_SCALE,
        U=DEFAULT_U,
        V=DEFAULT_V,
        W=DEFAULT_W,
        eta_order=DEFAULT_ETA_ORDER,
        profile=DEFAULT_PROFILE,
        pref_or=(
            DEFAULT_PREF_OR_PARAMS,
            DEFAULT_PREF_OR_HKL,
            DEFAULT_PREF_OR_METHOD,
            DEFAULT_PREF_OR_ELL,
        ),
        param_dict=None):
        self.phase_settings = phase_settings
        if self.phase_settings['recompute_peak_positions']:
            self.lattice_parameters = unit_cell.assemble_lattice_parameters(
                self.phase_settings)
            if param_dict is not None:
                self.lattice_parameters = self.copy_lp_urounds_from_param_dict(
                    param_dict)
        if self.phase_settings['preferred_orientation'] \
                and param_dict is not None:
            self.phase_settings['pref_orient_hkl'] \
                = param_dict['pref_orient_hkl']
            self.phase_settings['pref_orient_method'] \
                = param_dict['pref_orient_method']
            self.phase_settings['pref_orient_ell'] \
                = param_dict['pref_orient_ell']
            self.pref_orient = po.get_pref_orient_params(
                phase_settings,
                rp.get_param_from_dict(param_dict['pref_orient']))
        super(PhaseParameters, self).__init__(param_dict=param_dict)

        if param_dict is None:
            self.scale = scale
            self.U = U
            self.V = V
            self.W = W
            self.eta_order = eta_order
            self.eta = self.set_eta_order(self.eta_order)
            if self.phase_settings['preferred_orientation']:
                self.pref_orient = pref_or[0]
                self.phase_settings['pref_orient_hkl'] = pref_or[1]
                self.phase_settings['pref_orient_method'] = pref_or[2]
                self.phase_settings['pref_orient_ell'] = pref_or[3]
        assert profile in profiles.PROFILES
        self.profile = profile
        # self.profile = profiles.Profile(DEFAULT_PROFILE)
        # self.profile_parameters = [x for x in profile.]

    def set_eta_order(self, order):
        self.eta_order = self.validate_order(order,
            max_polynom_order=self.phase_settings["max_polynom_order"])
        self.eta = [x for x in self.eta_param_gen()]
        return self.eta

    def eta_param_gen(self):
        n = 0
        while n < self.eta_order:
            limit = np.power(0.01, n)
            if n == 0:
                yield ('eta_'+str(n), 0.5, [True], 0, 1)
            else:
                yield ('eta_'+str(n), 0.0, [True], -limit, limit)
            n += 1

    def param_gen(self):
        d = OrderedDict()
        d['scale'] = self.scale
        d['cagliotti_u'] = self.U
        d['cagliotti_v'] = self.V
        d['cagliotti_w'] = self.W
        d['eta'] = self.eta
        if self.phase_settings['recompute_peak_positions']:
            d['lattice_parameters'] = self.lattice_parameters
        if self.phase_settings['preferred_orientation']:
            d['pref_orient'] = self.pref_orient
        # if self.phase_settings['x']
        return d.iteritems()

    def from_dict(self, d):
        d2 = d.copy()
        d2['lattice_parameters'] = rp.as_dict(self.lattice_parameters)
        super(PhaseParameters, self).from_dict(d2)

    def copy_lp_urounds_from_param_dict(self, param_dict):
        from itertools import compress
        uc_mask = self.phase_settings['uc_mask']
        new_lps = []
        for i, lp in enumerate(self.lattice_parameters):
            lp = list(lp)
            lp[2] = list(compress(
                param_dict["lattice_parameters"], uc_mask))[i]['uround']
            new_lps.append(tuple(lp))
        return new_lps

    def as_dict(self):
        super_dict = super(PhaseParameters, self).as_dict()
        if isinstance(super_dict['pref_orient'], dict):
            super_dict['pref_orient'] = [super_dict['pref_orient']]
        for key in ('pref_orient_hkl', 'pref_orient_method', 'pref_orient_ell'):
            super_dict[key] = self.phase_settings[key]
        return super_dict
