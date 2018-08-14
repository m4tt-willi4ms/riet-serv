from __future__ import division, print_function, absolute_import
import os
import numpy as np
import json
from twisted.test import proto_helpers
from twisted.trial import unittest

from rietveld_server import RietveldServer
import rietveld_client as rc
from src.rietveld_phases import RietveldPhases as RP
import src.cctbx_dep.target_wavelengths as target_wavelengths

class RietveldServerTestCase(unittest.TestCase):
    def setUp(self):
        self.tr = proto_helpers.StringTransport()
        self.proto = RietveldServer()
        self.proto.makeConnection(self.tr)
        rc.load_json_samples()
        self.samples = rc.SAMPLES.copy()

        ref_model_no_input_profile = json.loads(self.samples['refinery_model'])
        ref_model_no_input_profile['input_data_path'] = ''
        ref_model_no_input_profile = json.dumps(ref_model_no_input_profile)
        self.samples['refinery_model_no_input'] = ref_model_no_input_profile

        ref_model_no_kalpha2 = json.loads(self.samples['refinery_model'])
        wavelengths = ref_model_no_kalpha2['wavelength_c']
        ref_model_no_kalpha2['wavelength_c'] = [wavelengths[0], 0.0]
        ref_model_no_kalpha2 = json.dumps(ref_model_no_kalpha2)
        self.samples['refinery_model_no_kalpha2'] = ref_model_no_kalpha2
        # self.proto._set_global_parameters(
        #     json.loads(self.samples["global_parameters"]))

    def _test(self, cmd, expected, result=None):
        cmd_parts = cmd.split(self.proto.split_character)
        self.tr.clear()
        d = getattr(self.proto, 'call_' + cmd_parts[0])(*cmd_parts[1:])
        assert self.tr.value() == expected

    def test_initialize(self):
        self._test('initialize', 'Initializing\n')

    def test_exit(self):
        self._test('exit', 'Goodbye.\n')

    def test__set_global_parameters(self):
        self.proto._set_global_parameters(
            json.loads(self.samples["global_parameters"]))
        assert len(RP.global_parameters.bkgd) == 3

    def test_add_phase(self):
        # sample = json.loads(self.phase_parameters_sample)
        # for key in sample.keys():
        #     if isinstance(sample[key], dict):
        #         print(sample[key]['uround'])
        #     elif isinstance(sample[key], list):
        #         for item in sample[key]:
        #             print(item.get('uround'), None)
        self.proto.call_initialize()
        self.tr.clear()
        self.proto.call_add_phase(self.samples['phase_parameters'])
        phase_reply_dict = json.loads(self.tr.value().strip())
        assert len(self.proto.phase_list) == 1
        for lp, param in zip(
            phase_reply_dict['lattice_parameters'],
            self.proto.phase_list[0].phase_settings["unit_cell"].parameters()):
            assert np.isclose(lp['value'], param)

    def test_initialize_update_refinery_model_then_add_phase(self):
        self.proto.call_initialize()
        self.proto._set_global_parameters(
            json.loads(self.samples['global_parameters']))
        self.proto.call_update_refinery_model(self.samples['refinery_model'])
        self.proto.call_add_phase(self.samples['phase_parameters'])

    def test_update_refinery_model(self):
        self.proto.call_initialize()
        self.proto.call_load_profile(
            self.samples['refinery_model'],
            self.samples['global_parameters'])
        self.proto.call_add_phase(self.samples['phase_parameters'])
        self.tr.clear()
        self.proto.call_update_refinery_model(self.samples['refinery_model'])
        reply = self.tr.value().strip()
        assert reply == "True" + self.proto.split_character
        assert self.proto.refinery_model['refinement_method'] == 'trf'
        assert 'max_polynomial_degree' in self.proto.refinery_model
        assert self.proto.phase_list[0].phase_settings['max_polynom_order'] == \
            self.proto.refinery_model['max_polynomial_degree'] + 1
        exp_wavelengths = target_wavelengths.set_wavelength(
            self.proto.refinery_model['wavelength'],
            self.proto.refinery_model['wavelength_model']
        )
        for x, y, exp in zip(
                self.proto.refinery_model['wavelength_c'],
                self.proto.phase_list[0].phase_settings['wavelengths'],
                exp_wavelengths):
            assert np.isclose(x, exp, rtol=1e-3)
            assert np.isclose(y, exp, rtol=1e-3)


    def test_update_refinery_model_no_phase_no_profile(self):
        self.proto.call_initialize()
        self.tr.clear()
        self.proto.call_update_refinery_model(
            self.samples['refinery_model_no_input'])
        reply = self.tr.value().strip()
        #Should reply 'False' to indicate no plot_data available
        assert reply == "False" + self.proto.split_character
        assert len(RP.two_theta) == 1000
        assert RP.I is None
        ref_model = json.loads(self.samples['refinery_model_no_input'])
        assert RP.phase_settings['max_polynom_order'] == \
            ref_model['max_polynomial_degree'] + 1

    def test_remove_phase(self):
        self.proto.call_initialize()
        self.proto.call_add_phase(self.samples['phase_parameters'])
        self.proto.call_update_refinery_model(
            self.samples['refinery_model_no_input'])
        self.proto.call_get_plot_data()

    def test_update_refinery_model_change_roi_then_get_plot_data(self):
        self.proto.call_initialize()
        self.proto.call_load_profile(
            self.samples['refinery_model'],
            self.samples['global_parameters'])
        rm_changed_roi = json.loads(self.samples['refinery_model'])
        rm_changed_roi['two_theta_roi_window'] = [10, 20]
        assert len(RP.two_theta) == 7250
        self.proto.call_update_refinery_model(json.dumps(rm_changed_roi))
        assert len(RP.two_theta) == 501
        self.tr.clear()
        self.proto.call_get_plot_data()
        plot_data = json.loads(self.tr.value().split(
            self.proto.split_character)[0])
        assert len(plot_data['intensities']) == 501

    def test_load_profile(self):
        self.proto.call_initialize()
        self.proto.call_add_phase(self.samples['phase_parameters'])
        assert len(self.proto.phase_list[0].two_theta) == 1000
        self.proto.call_load_profile(
            self.samples['refinery_model'],
            self.samples['global_parameters'])
        assert len(self.proto.phase_list[0].two_theta) != 1000

    def _check_single_phase_single_profile(self):
        assert len(self.proto.phase_list) == 1
        assert len(self.proto.phase_list[0].two_theta) != 1000
        assert len(RP.global_parameters.bkgd) == 3
        assert np.isclose(
            self.proto.phase_list[0].phase_data["crystal_density"],
            3.982641
            )

    def test_add_phase_then_load_profile(self):
        self.proto.call_initialize()
        self.proto.call_add_phase(self.samples['phase_parameters'])
        self.proto.call_load_profile(self.samples['refinery_model'],
                             self.samples['global_parameters'])
        self._check_single_phase_single_profile()

    def _check_single_phase_no_profile(self):
        assert len(self.proto.phase_list) == 1
        assert len(self.proto.phase_list[0].two_theta) == 1000
        assert len(RP.global_parameters.bkgd) == 3
        assert np.isclose(
            self.proto.phase_list[0].phase_data["crystal_density"],
            3.982641
            )

    def test_reset_with_phase_and_profile_loaded(self):
        self.proto.call_initialize()
        self.tr.clear()
        self.proto.call_load_profile(self.samples['refinery_model'],
            self.samples['global_parameters'])
        # print(self.tr.value())
        self.proto.call_initialize()
        self.proto.call_add_phase(self.samples['phase_parameters'])
        # self.proto.call_start_refine(self.samples['refinery_model']
        #     self.samples['rietveld_client'])
        self._check_single_phase_no_profile()

    def test_update_refinery_model_no_kalpha2(self):
        self.proto.call_initialize()
        self.proto.call_update_refinery_model(
            self.samples['refinery_model_no_kalpha2'])
        self.proto.call_add_phase(self.samples['phase_parameters'])

    def test_start_refine(self):
        self.proto.call_initialize()
        self.proto.call_update_refinery_model(self.samples['refinery_model'])
        self.proto.call_add_phase(self.samples['phase_parameters'])
        self.proto.call_start_refine(
            self.samples['refinery_model'], self.samples['rietveld_state_pp'])