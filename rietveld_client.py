"""
A test client to demonstrate some simple use of the Rietveld Server
"""
import time
import json
import os
from multiprocessing import Process

from twisted.internet import reactor
from twisted.internet.protocol import ReconnectingClientFactory
from twisted.protocols.basic import LineReceiver

with open(os.path.join(os.path.dirname(__file__),
    'data/server_input/rietveld_model_sample.json')) as file:
    ref_model = json.dumps(json.load(file))
with open(os.path.join(os.path.dirname(__file__),
    'data/server_input/global_parameters_sample.json')) as file:
    global_params = json.dumps(json.load(file))
with open(os.path.join(os.path.dirname(__file__),
    'data/server_input/phase_parameters_sample.json')) as file:
    phase_params = json.dumps(json.load(file))
with open(os.path.join(os.path.dirname(__file__),
    'data/server_input/rietveld_state_sample.json')) as file:
    rietveld_state_sample = json.dumps(json.load(file))

messages = [
    'help',
    # 'help reset',
    # 'load_profile .\data\profiles\d5_05005.xye',
    'reset',
    'load_profile;' + ref_model + ';' + global_params,
    'add_phase;' + phase_params,
    'add_phase;' + phase_params,
    'is_complete',
    'rounds_completed',
    'get_rietveld_state;2',
    'start_refine;' + ref_model + ';' + rietveld_state_sample,
    'is_complete',
    'rounds_completed',
    'is_complete',
    'rounds_completed',
    'get_plot_data',
    # 'add_phase;' + phase_info
    # 'get_phase_info',
    # 'get_phase_profile',
    # 'writeJSON',
    # 'loadJSON',
    # 'reset',
    # 'exit',
    ]

class RietveldClient(LineReceiver):
    """Once connected, send a message, then print the result."""
    data_buffer = ""

    def connectionMade(self):
        delay = 0.2
        if len(messages) > 0:
            msg = messages.pop(0)
            reactor.callLater(delay, self.sendLine, msg)
        else:
            # reactor.callLater(0, self.transport.loseConnection)
            reactor.stop()

    def dataReceived(self, data):
        self.data_buffer += data
        if data[-1] in ["\n"]:
            if not self.data_buffer[:9] == "Connected":
                print "Data received. Length:", len(self.data_buffer)
                self.transport.loseConnection()
                print self.data_buffer
            self.data_buffer = ""

    # def connectionLost(self, reason):
    #     print "connection lost"

class RietveldClientFactory(ReconnectingClientFactory):
    # def startedConnecting(self, connector):
    #     pass

    def buildProtocol(self, addr):
        self.resetDelay()
        return RietveldClient()

    # def clientConnectionLost(self, connector, reason):
    #     print 'Disconnected.'
    #     ReconnectingClientFactory.clientConnectionLost(self, connector, reason)

    def clientConnectionFailed(self, connector, reason):
        print 'Connection failed. Reason:', reason
        ReconnectingClientFactory.clientConnectionFailed(self, connector,
                                                         reason)


# this connects the protocol to a server running on port 8000
# def send_message_in_process():

def main():
    f = RietveldClientFactory()
    reactor.connectTCP("localhost", 8007, f)
    reactor.run()
        # import sys
        # del sys.modules['twisted.internet.reactor']
        # from twisted.internet import reactor

# this only runs if the module was *not* imported
if __name__ == '__main__':
    main()
