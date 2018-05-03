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

#     def test_img(self):
#         return self._test('img test.tif', '''New image: .\\test.tif
#    R: 85 mm
#    two-theta: 30 deg
# ''')

#     def test_calc_last(self):
#         self.proto.dataReceived('reset')
#         self.proto.dataReceived('img test.tif')
#         self.proto.dataReceived('roi  (320,370,140,170)')
#         self.proto.dataReceived('roi_label YBCO (0,0,1)')
#         self.tr.clear()
#         self.proto.dataReceived('calc_last')
#         self.assertEqual(self.tr.value(), 'Starting calculation...\n')

#     def test_no_imgs(self):
#         self.proto.call_reset()
#         self.tr.clear()
#         self.proto.dataReceived('roi (320,370,140,170)')
#         self.assertEqual(self.tr.value(), 'Error: no file names specified.\n')
#         self.tr.clear()

#     def test_roi(self):
#         return self._test('roi YBCO',
# '''Received:
#    Material: YBCO
#    hkl: (0, 0, 1)
#    ROI Coordinates: (320, 370, 140, 170)
#    Assigned to test.tif\r\n''')