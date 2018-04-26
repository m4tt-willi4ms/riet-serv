
# Copyright (c) Twisted Matrix Laboratories.
# See LICENSE for details.


"""
An example client. Run simpleserv.py first before running this.
"""
import time

from twisted.internet import reactor
from twisted.internet.protocol import ClientFactory
from twisted.protocols.basic import LineReceiver


class EchoClient(LineReceiver):
    """Once connected, send a message, then print the result."""
    messages = [
        'help',
        'help reset',
        # 'writeJSON',
        'loadJSON',
        'reset',
        # 'exit',
        ]

    def sendMessage(self, msg):
        self.sendLine(msg)

    def connectionMade(self):
        self.setRawMode()
        reactor.callLater(0,self.sendMessage, self.messages.pop(0))

    def rawDataReceived(self, data):
    #     # "As soon as any data is received, write it back."
        print "Received:", data
        self.clearLineBuffer()
        if len(self.messages) > 0:
            reactor.callLater(0.2,self.sendMessage, self.messages.pop(0))
        else:
            self.transport.loseConnection()

    def connectionLost(self, reason):
        print "connection lost"

class EchoFactory(ClientFactory):
    protocol = EchoClient

    def clientConnectionFailed(self, connector, reason):
        print "Connection failed - goodbye!"
        reactor.stop()

    def clientConnectionLost(self, connector, reason):
        print "Connection lost - goodbye!"
        reactor.stop()


# this connects the protocol to a server running on port 8000
def main():
    f = EchoFactory()
    reactor.connectTCP("localhost", 8007, f)
    reactor.run()

# this only runs if the module was *not* imported
if __name__ == '__main__':
    main()
