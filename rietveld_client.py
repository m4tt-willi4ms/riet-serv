
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
    # ref_model = json.dumps({
    #     # 'input_data_path': '.\\data\\profiles\\Jade-Al2O3-Sim.xye',
    #     'input_data_path': '.\\data\\profiles\\d5_05005.xye',
    #     'two_theta_roi_window': [25.0, 180.0],
    #     'wavelength': 'Cr',
    #     'wavelength_model': 0, #TODO: check 1
    #     'wavelengthc': 1.0
    #     })
    with open('./data/server_input/rietveld_model_sample.json') as file:
        ref_model = json.dumps(json.load(file))
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
    with open('./data/server_input/phase_parameters_sample.json') as file:
        phase_params = json.dumps(json.load(file))

    with open('./data/server_input/rietveld_state_sample.json') as file:
        rietveld_state_sample = json.dumps(json.load(file))

    messages = [
        'help',
        # 'help reset',
        # 'load_profile .\data\profiles\d5_05005.xye',
        'reset',
        'add_phase;' + phase_params,
        'add_phase;' + phase_params,
        'add_phase;' + phase_params,
        'add_phase;' + phase_params,
        'load_profile;' + ref_model + ';' + global_params,
        'is_complete',
        'rounds_completed',
        'get_rietveld_state;2',
        'start_refine;' + ref_model + ';' + rietveld_state_sample,
        'is_complete',
        'rounds_completed',
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
        if data[-1] in ["\n"]:
            delay = 0.5
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
