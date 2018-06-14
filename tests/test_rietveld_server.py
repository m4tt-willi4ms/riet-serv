from __future__ import division, print_function, absolute_import
import os
import numpy as np
import json
from twisted.test import proto_helpers
import pytest

from rietveld_server import RietveldServer
import rietveld_client as rc

@pytest.fixture(scope="module")
def tr():
    return proto_helpers.StringTransport()

@pytest.fixture(scope="module")
def server(tr):
    server = RietveldServer()
    server.makeConnection(tr)
    return server

def _test(tr, server, cmd, expected, result=None):
    cmd_parts = cmd.split(server.split_character)
    tr.clear()
    d = getattr(server, 'call_' + cmd_parts[0])(*cmd_parts[1:])
    assert tr.value() == expected

def test_initialize(tr, server):
    _test(tr, server, 'initialize', 'Initializing\n')

def test_exit(tr, server):
    _test(tr, server, 'exit', 'Goodbye.\n')


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
