from __future__ import division, print_function, absolute_import
import os
import numpy as np
import json
from twisted.test import proto_helpers
from twisted.trial import unittest

from rietveld_server import RietveldServer
import rietveld_client as rc
from src.rietveld_phases import RietveldPhases as RP

class RietveldServerTestCase(unittest.TestCase):
    def setUp(self):
        self.tr = proto_helpers.StringTransport()
        self.proto = RietveldServer()
        self.proto.makeConnection(self.tr)
        rc.load_json_samples()
        self.samples = rc.SAMPLES
        self.proto._set_global_parameters(json.loads(self.samples["global_parameters"]))

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
        self.proto._set_global_parameters(json.loads(self.samples["global_parameters"]))
        assert len(RP.global_parameters.bkgd) == 3

    def test_add_phase(self):
        # sample = json.loads(self.phase_parameters_sample)
        # for key in sample.keys():
        #     if isinstance(sample[key], dict):
        #         print(sample[key]['uround'])
        #     elif isinstance(sample[key], list):
        #         for item in sample[key]:
        #             print(item.get('uround'), None)
        self.proto.call_add_phase(self.samples['phase_parameters'])
        assert len(self.proto.phase_list) == 1




# class ClientCalculationTestCase(unittest.TestCase):

#     def setUp(self):
#         self.tr = proto_helpers.StringTransport()
#         self.proto = RietveldServer()
#         self.proto.makeConnection(self.tr)


# #     def test_help(self):
# #         return self._test('help',
# # '''Valid commands: exit, help, img, loadJSON, reset, roi, update_R, writeJSON\r\n''')

#     def test_exit(self):
#         return self._test('exit', 'Goodbye.\n')

#     def test_initialize(self):
#         return self._test('initialize', 'Initializing\n')

#     def test_add_phase(self):
#         self.tr.clear()
#         self.proto.call_add_phase(rc.SAMPLES['phase_parameters'])
#         input_phase = json.loads(rc.SAMPLES['phase_parameters'])
#         output_phase = json.loads(
#             self.tr.value().split(self.proto.split_character)[0])
#         print json.dumps(output_phase, indent=4)
#         check_val = ('cif_path',
#                      'lattice_parameter_tolerances',
#                      )
#         check_len = ('eta', 'lattice_parameter_tolerances')
#         for key in check_len:
#             assert len(output_phase[key]) == len(input_phase[key])
#         for key in check_val:
#             if isinstance(output_phase[key], list):
#                 for i, item in enumerate(output_phase[key]):
#                     assert np.isclose(item, input_phase[key][i])
#             if isinstance(output_phase[key], str):
#                 assert output_phase[key] == input_phase[key]
