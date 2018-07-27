from __future__ import division, print_function, absolute_import
from collections import OrderedDict
import copy
import json
import numpy as np
from gutter.client.default import gutter

from src.refinement_parameters import RefinementParameters, as_dict
import src.profiles as profiles
import src.cctbx_dep.unit_cell as unit_cell

DEFAULT_U = ('cagliotti_u', 0.00, [False], -0.1, 0.1)
DEFAULT_V = ('cagliotti_v', 0.00, [False], -0.1, 0.1)
DEFAULT_W = ('cagliotti_w', 0.001, [True], 0.000001, 1)
DEFAULT_SCALE = ('scale', 0.1, [True], 0, float('inf'))
DEFAULT_PREF_OR_PARAMS = [('pref_or_0', 1.00, [True], 0.000001, 2)]
DEFAULT_PREF_OR_HKL = [0, 0, 1]
DEFAULT_PREF_OR_METHOD = 'march_dollase'
DEFAULT_ETA_ORDER = 2
DEFAULT_LATTICE_DEV = 0.01
DEFAULT_PROFILE = 'PV'

class PhaseParameters(RefinementParameters):
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
        ),
        param_dict=None):
        # RefinementParameters.__init__(self)
        self.phase_settings = phase_settings
        self.recompute_peak_positions = \
            self.phase_settings['recompute_peak_positions']
        self.preferred_orientation = \
            self.phase_settings['preferred_orientation']
        # if self.phase_settings["recompute_peak_positions"]:
        if gutter.active('rpp', self):
            self.lattice_parameters = unit_cell.assemble_lattice_parameters(
                self.phase_settings)
            if param_dict is not None:
                self.lattice_parameters = self.copy_lp_urounds_from_param_dict(
                    param_dict)
            # print("from assemble lp:", lps)
            # print("from self.lps:", self.lattice_parameters)
        # else:
        #     self.eta = self.set_eta_order(len(self.eta))
        super(PhaseParameters, self).__init__(param_dict=param_dict)
        if param_dict is None:
            self.scale = scale
            self.U = U
            self.V = V
            self.W = W
            self.eta_order = eta_order
            self.eta = self.set_eta_order(self.eta_order)
            if gutter.active('pref_or', self):
                self.pref_or = pref_or[0]
                self.phase_settings['pref_orient_hkl'] = pref_or[1]
                self.phase_settings['pref_orient_method'] = pref_or[2]
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
        if gutter.active('rpp', self):
            d['lattice_parameters'] = self.lattice_parameters
        if gutter.active('pref_or', self):
            d['pref_or'] = self.pref_or
        # if self.phase_settings['x']
        return d.iteritems()

    def from_dict(self, d):
        d2 = d.copy()
        d2['lattice_parameters'] = as_dict(self.lattice_parameters)
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
