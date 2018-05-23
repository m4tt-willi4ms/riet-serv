
# Copyright (c) Twisted Matrix Laboratories.
# See LICENSE for details.


"""
An example client. Run simpleserv.py first before running this.
"""
import time
import json

from twisted.internet import reactor
from twisted.internet.protocol import ClientFactory
from twisted.protocols.basic import LineReceiver


class RietveldClient(LineReceiver):
    """Once connected, send a message, then print the result."""
    ref_model = json.dumps({
        # 'input_data_path': '.\\data\\profiles\\Jade-Al2O3-Sim.xye',
        'input_data_path': '.\\data\\profiles\\d5_05005.xye',
        'two_theta_roi_window': [25.0, 180.0],
        'wavelength': 'Cr',
        'wavelength_model': 2, #TODO: check 1
        })
    global_params = json.dumps({
        'two_theta_offset': {
        'name': 'two_theta_0',
        'value': 0.0,
        'uround': [True],
        'l_limit': -0.1,
        'u_limit': -0.1,
        },
        'bkgd': [{
        'name': 'bkgd_0',
        'value': 0.0,
        'uround': [True],
        'l_limit': -float('inf'),
        'u_limit': float('inf'),
        }]
        })
    phase_params = json.dumps({"lattice_parameter_tolerances":[0.0,0.0,0.0,0.0,0.0,0.0],"cif_path":"C:\\Users\\prati\\Desktop\\Rietveld_CCTBX\\data\\cifs\\1000032.cif","phase_name":"","scale":{"name":"scale","value":0.0,"round":1,"used":true,"uround":[false,false,false,false,false,false,false,false,false,false,false,false],"l_limit":0.0,"u_limit":3.40282347E+38},"cagliotti_u":{"name":"cagliotti_u","value":0.0,"round":4,"used":false,"uround":[false,false,false,false,false,false,false,false,false,false,false,false],"l_limit":-0.1,"u_limit":0.1},"cagliotti_v":{"name":"cagliotti_v","value":0.0,"round":4,"used":false,"uround":[false,false,false,false,false,false,false,false,false,false,false,false],"l_limit":-0.1,"u_limit":0.1},"cagliotti_w":{"name":"cagliotti_w","value":0.003,"round":3,"used":true,"uround":[false,false,false,false,false,false,false,false,false,false,false,false],"l_limit":1E-05,"u_limit":1.0},"eta":[{"name":"eta_0","value":0.5,"round":2,"used":true,"uround":[false,false,false,false,false,false,false,false,false,false,false,false],"l_limit":0.0,"u_limit":1.0},{"name":"eta_1","value":0.5,"round":4,"used":false,"uround":[false,false,false,false,false,false,false,false,false,false,false,false],"l_limit":0.0,"u_limit":1.0},{"name":"eta_2","value":0.5,"round":4,"used":false,"uround":[false,false,false,false,false,false,false,false,false,false,false,false],"l_limit":0.0,"u_limit":1.0}],"lattice_parameters":[{"name":"lattice_parameters_0","value":0.0,"round":2,"used":false,"uround":[false,false,false,false,false,false,false,false,false,false,false,false],"l_limit":-3.40282347E+38,"u_limit":3.40282347E+38},{"name":"lattice_parameters_1","value":0.0,"round":2,"used":false,"uround":[false,false,false,false,false,false,false,false,false,false,false,false],"l_limit":-3.40282347E+38,"u_limit":3.40282347E+38},{"name":"lattice_parameters_2","value":0.0,"round":2,"used":false,"uround":[false,false,false,false,false,false,false,false,false,false,false,false],"l_limit":-3.40282347E+38,"u_limit":3.40282347E+38},{"name":"lattice_parameters_3","value":0.0,"round":2,"used":false,"uround":[false,false,false,false,false,false,false,false,false,false,false,false],"l_limit":-3.40282347E+38,"u_limit":3.40282347E+38},{"name":"lattice_parameters_4","value":0.0,"round":2,"used":false,"uround":[false,false,false,false,false,false,false,false,false,false,false,false],"l_limit":-3.40282347E+38,"u_limit":3.40282347E+38},{"name":"lattice_parameters_5","value":0.0,"round":2,"used":false,"uround":[false,false,false,false,false,false,false,false,false,false,false,false],"l_limit":-3.40282347E+38,"u_limit":3.40282347E+38}]})
    messages = [
        'help',
        # 'help reset',
        # 'load_profile .\data\profiles\d5_05005.xye',
        'load_profile;' + ref_model + ';' + global_params,
        'add_phase;' + phase_params,
        # 'add_phase;' + phase_info
        # 'get_phase_info',
        # 'get_phase_profile',
        # 'writeJSON',
        # 'loadJSON',
        # 'reset',
        # 'exit',
        ]
    data_buffer = ""

    def sendMessage(self, msg):
        self.sendLine(msg)

    def connectionMade(self):
        pass
        # self.sendLine('help')
        # self.setRawMode()
        # print self.MAX_LENGTH
        # self.MAX_LENGTH = 4194304

    # def lineReceived(self, data):


    def dataReceived(self, data):
        print "Data received. Length:", len(data)
        self.data_buffer += data
        delay = 0.
        if data[-1] in [";", "\n"]:
            if len(self.messages) > 0:
                msg = self.messages.pop(0)
                self.clearLineBuffer()
                reactor.callLater(delay, self.sendMessage, msg)
            else:
                reactor.callLater(delay, self.transport.loseConnection)
            print self.data_buffer
            # print "Message sent."
            self.data_buffer = ""

    # def connectionLost(self, reason):
    #     print "connection lost"

class RietveldClientFactory(ClientFactory):
    protocol = RietveldClient

    def clientConnectionFailed(self, connector, reason):
        print "Connection failed - goodbye!"
        reactor.stop()

    def clientConnectionLost(self, connector, reason):
        print "Connection lost - goodbye!"
        reactor.stop()

    # def startedConnecting(self, connector):
    #     pass


# this connects the protocol to a server running on port 8000
def main():
    f = RietveldClientFactory()
    reactor.connectTCP("localhost", 8007, f)
    reactor.run()

# this only runs if the module was *not* imported
if __name__ == '__main__':
    main()
