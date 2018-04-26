
# Copyright (c) Twisted Matrix Laboratories.
# See LICENSE for details.


"""
An example client. Run simpleserv.py first before running this.
"""
import time

from twisted.internet import reactor
from twisted.internet.protocol import ClientFactory
from twisted.protocols.basic import LineReceiver


# a client protocol

class EchoClient(LineReceiver):
    """Once connected, send a message, then print the result."""
    messages = [
        'help',
        'help reset',
        # 'writeJSON',
        'loadJSON',
        'reset',
        'help roi',
        'img test001.tif',
        'update_R 100',
        'img test002.tif',
        'update_R 60',
        'update_two_theta 30',
        'help update_wavelength',
        'update_wavelength 6.5',
        'update_omega 13',
        'img_dir .\\tests',
        'img image_00171.tif',
        'roi (320,370,140,170)',
        'roi_label YBCO (0,0,1)',
        'roi_label REO (1,0,2)',
        'calc_last',
        'is_calculating',
        'is_calculating',
        'is_calculating',
        'is_calculating',
        'writeJSON',
        'exit',
        ]

    def sendMessage(self, msg):
        self.sendLine(msg)

    def connectionMade(self):
        # d = defer.Deferred()
            # d.addCallback(self.sendMessage)
            # d.addErrback(self.printError)
        # for msg in self.messages:
        self.setRawMode()
        reactor.callLater(0,self.sendMessage, self.messages.pop(0))

    def rawDataReceived(self, data):
    #     # "As soon as any data is received, write it back."
        print "Received:", data
        time.sleep(0.05)
        self.clearLineBuffer()
        if len(self.messages) > 0:
            self.sendMessage(self.messages.pop(0))
        else:
            self.transport.loseConnection()
        # pass
        # self.transport.loseConnection()

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
    reactor.connectTCP("localhost", 8006, f)
    reactor.run()

# this only runs if the module was *not* imported
if __name__ == '__main__':
    main()
