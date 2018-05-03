
# Copyright (c) Twisted Matrix Laboratories.
# See LICENSE for details.


"""
An example client. Run simpleserv.py first before running this.
"""
import time

from twisted.internet import reactor
from twisted.internet.protocol import ClientFactory
from twisted.protocols.basic import LineReceiver


class RietveldClient(LineReceiver):
    """Once connected, send a message, then print the result."""
    messages = [
        'help',
        # 'help reset',
        'load_profile .\data\profiles\d5_05005.xye',
        'add_phase .\data\cifs\9015662-rutile.cif',
        'get_phase_info',
        'get_phase_profile',
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
