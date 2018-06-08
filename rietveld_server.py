from twisted.internet import reactor, protocol, defer, task
from twisted.protocols import basic
from twisted.python import log

import json, simplejson
import sys, os

import src.rietveld_phases as rp
import src.rietveld_refinery as rr

class RietveldServer(basic.LineReceiver):
    delimiter = b'\n'
    split_character = ';'
    calc_flag = False
    err_flag = False
    phase_list = []
    # phase_dict_list = []
    refinery_model = None
    rietveld_refinery = None
    rietveld_history = []

    def connectionMade(self):
        # pass
        self.sendLine(b'Connected. Type <help> [command] for info.')
        self.MAX_LENGTH = 4194304
        # self.setLineMode()
        # print self.__dict__
        # self.sendLine(b"Web checker console. Type 'help' for help.")

    def lineReceived(self, data):
        if not data: return
        data = data.decode("ascii")

        # Parse the command
        commandParts = [x.strip() for x in data.split(self.split_character)]
        print commandParts
        command = commandParts[0]#.lower()
        args = commandParts[1:]

        try:
            method = getattr(self, 'call_' + command)
        except AttributeError as e:
            self.sendLine(b'Error: no such command.')
        else:
            try:
                print command, args
                method(*args)
            except Exception as e:
                # pass
                self.sendLine(b'Error: ' + str(e).encode("ascii"))

    def call_exit(self):
        """exit: shuts down the TCP server"""
        self.sendLine(b'Goodbye.')
        self.transport.loseConnection()
        # reactor.stop()

    # def get_target():

    #    self.sendLine()

    def call_set_wavelength(self, target, mode):
        try:
            rp.RietveldPhases.set_wavelength(target, int(mode))

        except:
            log.err()

    def _set_refinery_model(self, ref_model):
        self.refinery_model = ref_model
        profile_path = self.refinery_model['input_data_path']
        min_two_theta = self.refinery_model['two_theta_roi_window'][0]
        max_two_theta = self.refinery_model['two_theta_roi_window'][1]
        rp.RietveldPhases.set_profile(profile_path,
            min_two_theta=min_two_theta,
            max_two_theta=max_two_theta,
            lines_to_strip_at_TOF=3,
            )
        wavelength_string = self.refinery_model["wavelength"]
        wavelength_model = self.refinery_model["wavelength_model"]
        custom_wavelength = self.refinery_model["wavelengthc"]
        rp.RietveldPhases.set_wavelength(wavelength_string, wavelength_model,
            custom_wavelength=custom_wavelength)

    def _set_global_parameters(self, global_parameters):
        rp.RietveldPhases.global_parameters.from_dict(global_parameters)
        # profile_filename = os.path.split(profile_path)[1]
        # self.sendLine(b'Loaded the profile {0}'.format(profile_filename))
        if rp.RietveldPhases.I is not None:
            self._bkgd_refine()

    def _bkgd_refine(self):
        self.rietveld_refinery = rr.RietveldRefinery(self.phase_list,
            bkgd_refine=True)
        self.rietveld_refinery.minimize()

    def call_load_profile(self, refinery_model_string, global_parameters):
        try:
            assert isinstance(refinery_model_string, unicode)
            self._set_refinery_model(json.loads(refinery_model_string))
            phase_temp = list(self.phase_list)
            self.phase_list = []
            for phase in phase_temp:
                self._add_phase(phase.as_dict())

            assert isinstance(global_parameters, unicode)
            self._set_global_parameters(json.loads(global_parameters))

            self.rietveld_refinery = rr.RietveldRefinery(self.phase_list)
            profile = self.rietveld_refinery.total_profile()
            rietveld_plot = rp.RietveldPhases.get_rietveld_plot(profile,
                compute_differences=True)
            reply = ""
            reply += json.dumps(rietveld_plot, indent=4) + ";"
            global_parameters = rp.RietveldPhases.global_parameters.as_dict()
            reply += json.dumps(global_parameters, indent=4) + ";"
            print "Message Length (in bytes):", len(reply.encode('utf-8'))
            self.sendLine(reply)
        except:
            log.err()

    # def _fit_added_phase(self, phase_parameters_JSON):

    def _add_phase(self, phase_dict):
        cif_path = phase_dict["cif_path"]
        # except KeyError:
        #    cif_path = phase_dict["input_cif_path"]
        assert type(cif_path) == unicode
        self.phase_list.append(rp.RietveldPhases(cif_path,
            phase_parameter_dict=phase_dict))

    def call_add_phase(self, phase_parameters_JSON):
        """add_phase: appends a phase to the server's phase list, given a
PhaseParameters object in json-serialized form.
        """
        try:
            phase_dict = json.loads(phase_parameters_JSON)
            self._add_phase(phase_dict)
            if rp.RietveldPhases.I is not None:
                self._bkgd_refine()

            phase_dict = json.dumps(self.phase_list[-1].as_dict(), indent=4)
            reply = ""
            reply += phase_dict + ";"
            if self.rietveld_refinery is not None:
                plot_data = json.dumps(rp.RietveldPhases.get_plot_data(
                    self.rietveld_refinery.total_profile()), indent=4)
            else:
                plot_data = json.dumps(rp.RietveldPhases.get_plot_data(
                    self.phase_list[-1].phase_profile()), indent=4)
            reply += plot_data + ";"
            # print len(reply.encode('utf-8'))
            self.sendLine(reply)
            # print 'Here', self.MAX_LENGTH
        # d = defer.Deferred()
        # d.addCallback(self._fit_added_phase)
        # d.addErrback(self._error)
        # # d.addCallback(lambda x: self._calc_complete())
        # reactor.callInThread(d.callback, json_string)
        except:
            log.err()

    # def call_get_current_profile(self):
    #    if self.rietveld_refinery is not None:

    def _calc_complete(self):
        self.calc_flag = False
        state = {}
        state['rietveld_data'] = self.rietveld_refinery.get_plot_data()
        state['global_state'] = rp.RietveldPhases.global_parameters.as_dict()
        state['phase_state'] = [phase.as_dict() for phase in self.phase_list]
        self.rietveld_history.append(state)

    def _refine_error(self):
        self.err_flag = True
        log.err()

    def _refine(self):
        self.rietveld_refinery.minimize()

    def _run(self):
        self.calc_flag = True
        self.err_flag = False
        d = defer.Deferred()
        d.addCallback(lambda x: self._refine())
        d.addErrback(lambda x: self._refine_error())
        d.addCallback(lambda x: self._calc_complete())
        reactor.callInThread(d.callback, None)

    def call_start_refine(self, refinery_model, rietveld_state):
        """start_refine(refinery_model, rietveld_state): returns a true or false
verifying that a refinement session (possibly consisting of several rounds) has
begun. The json-serialized input refinery_model, rietveld_state objects are used
to specify the refinement options and starting point of the refinement,
respectively
        """
        try:
            self._set_refinery_model(json.loads(refinery_model))

            rs = json.loads(rietveld_state)
            self._set_global_parameters(rs['global_state'])
            self.phase_list = []
            for phase in rs['phase_state']:
                self._add_phase(phase)

            self._bkgd_refine()

            self.rietveld_refinery = rr.RietveldRefinery(self.phase_list)
            self.rietveld_refinery.set_mask(['two_theta_0', 'bkgd', 'scale'])
            self._run()

            self.sendLine(str(True) + ";")
        except:
            log.err()
            self.sendLine(str(False) + ";")


    def call_rounds_completed(self):
        """rounds_completed: returns the number of rounds completed to date by
the refinement engine
        """
        self.sendLine(str(len(self.rietveld_history)) + ";")

    def call_get_rietveld_state(self, round_number=-1):
        """get_rietveld_state [round_number]: returns the json-serialized
rietveld_state object corresponding to the end of the round specified. (If no
round is specified, calling this method returns the last entry found in
rietveld_history.)
        """
        round_number = int(round_number)
        max_round = len(self.rietveld_history)
        if max_round > 0 and round_number < max_round:
            self.sendLine(json.dumps(self.rietveld_history[round_number]) + ";")
        else:
            self.sendLine(str(False) + ";")

    def call_is_complete(self):
        """is_complete: returns either true or false, depending on whether or
not the rietveld_refinement session has completed
        """
        self.sendLine(str(not self.calc_flag) + ";")

    def call_can_ping(self):
        """can_ping: returns True (for diagnostic purposes)
        """
        self.sendLine(str(True) + ";")

    def call_get_phase_profile(self, index=u'-1'):
        """get_phase_profile [index]: returns a json-serialized list containing
the phase profile data. If no index is specified, information for the
most-recently loaded phase is returned"""
        index = int(index)
        profile = json.dumps(list(self.phase_list[index].phase_profile()),
            indent=4)
        self.sendLine(profile)

    def call_reset(self):
        """reset: returns the rietveld_server to its initial state"""
        self.calc_flag = False
        self.err_flag = False
        self.phase_list = []
        # phase_dict_list = []
        self.refinery_model = None
        self.rietveld_refinery = None
        self.rietveld_history = []
        rp.RietveldPhases.global_parameters.reset_x()
        self.sendLine(b'Resetting')

    def call_help(self, command=None):
        """help [command]: List commands, or show help on the given command"""
        if command:
            doc = getattr(self, 'call_' + command).__doc__
            self.sendLine(doc.encode("ascii"))
        else:
            commands = [cmd[5:].encode("ascii")
                            for cmd in dir(self)
                            if cmd.startswith('call_')]
            self.sendLine(b"Valid commands: \r\n  " + "\r\n  ".join(commands))

    def call_is_calculating(self):
        if self.calc_flag:
            self.sendLine(b'Analysis in progress...')
        else:
            if self.err_flag:
                self.sendLine(b'Error in running analysis')
            else:
                self.sendLine(b'Analysis complete')

    def connectionLost(self, reason):
        pass
        # reactor.stop()

def main():
    factory = protocol.ServerFactory()
    factory.protocol = RietveldServer
    log.startLogging(sys.stdout)
    reactor.listenTCP(8007, factory)
    reactor.run()

# this only runs if the module was *not* imported
if __name__ == '__main__':
    main()