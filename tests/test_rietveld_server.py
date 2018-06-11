from twisted.trial import unittest
from twisted.test import proto_helpers

from rietveld_server import RietveldServer


class ClientCalculationTestCase(unittest.TestCase):
    def setUp(self):
        self.tr = proto_helpers.StringTransport()
        self.proto = RietveldServer()
        self.proto.makeConnection(self.tr)

    def _test(self, cmd, expected, result=None):
        cmd_parts = cmd.split()
        self.tr.clear()
        d = getattr(self.proto, 'call_' + cmd_parts[0])(*cmd_parts[1:])
        self.assertEqual(self.tr.value(), expected)
        return d

#     def test_help(self):
#         return self._test('help',
# '''Valid commands: exit, help, img, loadJSON, reset, roi, update_R, writeJSON\r\n''')

    def test_exit(self):
        return self._test('exit', 'Goodbye.\n')

    def test_reset(self):
        return self._test('reset', 'Resetting\n')