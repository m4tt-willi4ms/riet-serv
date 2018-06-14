"""
A test client to demonstrate some simple use of the Rietveld Server
"""
from __future__ import division, print_function, absolute_import
import time
import json
import os
from multiprocessing import Process

from twisted.internet import reactor
from twisted.internet.protocol import ReconnectingClientFactory
from twisted.protocols.basic import LineReceiver

SAMPLES = {
    'rietveld_model': None,
    'global_parameters': None,
    'phase_parameters': None,
    'rietveld_state': None,
}

def get_sample_path(key):
    return os.path.join(os.path.dirname(__file__),
            'data/server_input/' + key + '_sample.json')

def load_json_sample(key):
    with open(get_sample_path(key), 'r') as file:
        SAMPLES[key] = json.dumps(json.load(file))

def write_json_sample(key):
    with open(get_sample_path(key), 'w') as file:
        file.write(json.dumps(json.loads(SAMPLES[key]), indent=4))

def load_json_samples():
    for key in SAMPLES.keys():
        load_json_sample(key)

def write_json_samples():
    for key in SAMPLES.keys():
        write_json_sample(key)

load_json_samples()

MESSAGES = [
    'help',
    # 'help reset',
    # 'load_profile .\data\profiles\d5_05005.xye',
    'initialize',
    'add_phase;' + SAMPLES['phase_parameters'],
    'add_phase;' + SAMPLES['phase_parameters'],
    'load_profile;' + SAMPLES['rietveld_model'] + ';' + SAMPLES['global_parameters'],
    'is_complete',
    'rounds_completed',
    'get_rietveld_state;2',
    'start_refine;' + SAMPLES['rietveld_model'] + ';' + SAMPLES['rietveld_state'],
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
        delay = 0.02
        if len(MESSAGES) > 0:
            msg = MESSAGES.pop(0)
            reactor.callLater(delay, self.sendLine, msg)
        else:
            # reactor.callLater(0, self.transport.loseConnection)
            reactor.stop()

    def dataReceived(self, data):
        self.data_buffer += data
        if data[-1] in ["\n"]:
            if not self.data_buffer[:9] == "Connected":
                print("Data received. Length:", len(self.data_buffer))
                self.transport.loseConnection()
                print(self.data_buffer)
            self.data_buffer = ""

    # def connectionLost(self, reason):
    #     print "connection lost"

class RietveldClientFactory(ReconnectingClientFactory):
    def buildProtocol(self, addr):
        self.resetDelay()
        return RietveldClient()

    # def clientConnectionLost(self, connector, reason):
    #     print 'Disconnected.'
    #     ReconnectingClientFactory.clientConnectionLost(self, connector, reason)

    def clientConnectionFailed(self, connector, reason):
        print('Connection failed. Reason:', reason)
        ReconnectingClientFactory.clientConnectionFailed(self, connector,
                                                         reason)

def main():
    f = RietveldClientFactory()
    reactor.connectTCP("localhost", 8007, f)
    reactor.run()

if __name__ == '__main__':
    main()
