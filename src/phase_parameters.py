from __future__ import division, print_function, absolute_import
from collections import OrderedDict
import numpy as np
import json

from src.refinement_parameters import RefinementParameters
import src.profiles as profiles
import src.cctbx_dep.unit_cell as unit_cell

DEFAULT_U = ('cagliotti_u', 0.00, [False], -0.1, 0.1)
DEFAULT_V = ('cagliotti_v', 0.00, [False], -0.1, 0.1)
DEFAULT_W = ('cagliotti_w', 0.007, [True], 0.000001, 1)
DEFAULT_SCALE = ('scale', 0.1, [True], 0, float('inf'))
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
        param_dict=None):
        # RefinementParameters.__init__(self)
        super(PhaseParameters, self).__init__(param_dict=param_dict)
        self.phase_settings = phase_settings
        if param_dict is None:
            self.scale = scale
            self.U = U
            self.V = V
            self.W = W
            self.eta_order = eta_order
            self.eta = self.set_eta_order(self.eta_order)
        # else:
        #     self.eta = self.set_eta_order(len(self.eta))
        if self.phase_settings["recompute_peak_positions"]:
            self.lattice_parameters = unit_cell.assemble_lattice_parameters(
                self.phase_settings)
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
            limit = np.power(0.001, n)
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
        if self.phase_settings["recompute_peak_positions"]:
            d['lattice_parameters'] = self.lattice_parameters
        return d.iteritems()

    def load_json(json_str):
        json_dict = json.loads(json_str)
        self.phase_settings["cif_path"] = json_dict["cif_path"]