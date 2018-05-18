
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
    phase_info = json.dumps({'input_cif_path': '.\data\cifs\9015662-rutile.cif'})
    ref_model = json.dumps({
        'input_data_path': '.\data\profiles\cement_15_03_11_0028.xye',
        'two_theta_roi_window': [25.0, 180.0],
        })
    global_params = json.dumps({
        'two_theta_offset': ('two_theta_0', 0.0, [True], -0.1, 0.2),
        'bkgd': [{
        'name': 'bkgd_0',
        'value': 0.0,
        'u_round': [True],
        'l_limit': -float('inf'),
        'u_limit': float('inf'),
        }]
        })
    messages = [
        'help',
        # 'help reset',
        # 'load_profile .\data\profiles\d5_05005.xye',
        'load_profile;' + ref_model + "; " + global_params
        # 'add_phase;' + phase_info
        # 'get_phase_info',
        # 'get_phase_profile',
        # 'writeJSON',
        # 'loadJSON',
        # 'reset',
        # 'exit',
        ]

    def sendMessage(self, msg):
        self.sendLine(msg)

    def connectionMade(self):
        self.setRawMode()

    def rawDataReceived(self, data):
        print "Received:", data
        self.clearLineBuffer()
        delay = 0.2
        if len(self.messages) > 0:
            reactor.callLater(delay,self.sendMessage, self.messages.pop(0))
        else:
            reactor.callLater(delay,self.transport.loseConnection)

    def connectionLost(self, reason):
        print "connection lost"

class RietveldClientFactory(ClientFactory):
    protocol = RietveldClient

    def clientConnectionFailed(self, connector, reason):
        print "Connection failed - goodbye!"
        reactor.stop()

    def clientConnectionLost(self, connector, reason):
        print "Connection lost - goodbye!"
        reactor.stop()


# this connects the protocol to a server running on port 8000
def main():
    f = RietveldClientFactory()
    reactor.connectTCP("localhost", 8007, f)
    reactor.run()

# this only runs if the module was *not* imported
if __name__ == '__main__':
    main()
