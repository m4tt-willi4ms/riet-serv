from __future__ import division, print_function, absolute_import
import json
import sys, os
import time
import numpy as np
import copy

from twisted.internet import reactor, protocol, defer
from twisted.protocols import basic
from twisted.python import log

import src.rietveld_phases as rp
import src.rietveld_refinery as rr
from src.rietveld_plot import RietveldPlot

MAX_PROFILE_VALUE = 10000000

class RietveldServer(basic.LineReceiver):
    delimiter = b'\n'
    split_character = ';'
    """The character used to separate commands, arguments when sent to the
server"""
    calc_flag = False
    err_flag = False
    phase_list = []
    # phase_dict_list = []
    refinery_model = None
    rietveld_refinery = None
    rietveld_history = []
    plot_data = None
    # plot = RietveldPlot()
    show_plot = False
    count = 0
    two_thetas = None
    wavelengths = None

    def connectionMade(self):
        # pass
        # self.sendLine(b'Connected. Type <help> [command] for info.')
        self.MAX_LENGTH = 4194304
        # self.setLineMode()
        # print self.__dict__
        # self.sendLine(b"Web checker console. Type 'help' for help.")

    def lineReceived(self, data):
        if not data: return
        data = data.decode("ascii")

        # Parse the command
        commandParts = [x.strip() for x in data.split(self.split_character)]
        command = commandParts[0]#.lower()
        args = commandParts[1:]

        try:
            method = getattr(self, 'call_' + command)
        except AttributeError as e:
            self.sendLine(b'Error: no such command.')
        else:
            try:
                print(command, args)
                method(*args)
            except Exception as e:
                # pass
                self.sendLine(b'Error: ' + str(e).encode("ascii"))

    def call_exit(self):
        """exit: shuts down the TCP server"""
        self.sendLine(b'Goodbye.')
        self.transport.loseConnection()

    def _set_refinery_model(self, ref_model):
        RietveldServer.refinery_model = ref_model
        profile_path = RietveldServer.refinery_model['input_data_path']
        min_two_theta = RietveldServer.refinery_model['two_theta_roi_window'][0]
        max_two_theta = RietveldServer.refinery_model['two_theta_roi_window'][1]
        if profile_path:
            rp.RietveldPhases.set_profile(profile_path,
                min_two_theta=min_two_theta,
                max_two_theta=max_two_theta,
                lines_to_strip_at_tof=3,
                )
        RietveldServer.wavelengths = \
            RietveldServer.refinery_model["wavelength_c"]
        if np.isclose(RietveldServer.wavelengths[-1], 0.0):
            RietveldServer.wavelengths = [RietveldServer.wavelengths[0]]
        rp.RietveldPhases.set_wavelengths(RietveldServer.wavelengths)
        if RietveldServer.phase_list is not None:
            for i, phase in enumerate(RietveldServer.phase_list):
                phase = self._add_phase(phase.as_dict(), i)

    def _set_global_parameters(self, global_parameters):
        rp.RietveldPhases.set_global_parameters(global_parameters)
        # profile_filename = os.path.split(profile_path)[1]
        # self.sendLine(b'Loaded the profile {0}'.format(profile_filename))
        # if rp.RietveldPhases.I is not None \
        #         and len(RietveldServer.phase_list) > 0:
        #     self._bkgd_refine()

    def _bkgd_refine(self):
        if rp.RietveldPhases.I is not None:
            bkgd_refinery = rr.RietveldRefinery(
                RietveldServer.phase_list, bkgd_refine=True)
            bkgd_refinery.minimize()
            RietveldServer.plotdata = bkgd_refinery.get_plot_data()

    def call_load_profile(self, refinery_model_string, global_parameters):
        try:
            assert isinstance(refinery_model_string, basestring)
            self._set_refinery_model(json.loads(refinery_model_string))
            RietveldServer.two_thetas = None
            phase_temp = copy.deepcopy(RietveldServer.phase_list)
            RietveldServer.phase_list = []

            assert isinstance(global_parameters, basestring)
            self._set_global_parameters(json.loads(global_parameters))

            for phase in phase_temp:
                self._add_phase(phase.as_dict())

            if RietveldServer.show_plot:
                RietveldServer.plot.setplotdata()
            self._update_plot_data()
            rietveld_plot = rp.RietveldPhases.get_rietveld_plot(
                RietveldServer.plot_data)
            reply = ""
            reply += json.dumps(rietveld_plot) + self.split_character
            global_parameters = rp.RietveldPhases.global_parameters.as_dict()
            reply += json.dumps(global_parameters) + self.split_character
            print("Message Length (in bytes):", len(reply.encode('utf-8')))
            # with open(os.path.join(
            #         os.path.dirname(__file__),'reply.json'), 'w') as file:
            #     file.write(json.dumps(rietveld_plot))
            self.sendLine(reply)
        except:
            log.err()

    # def _fit_added_phase(self, phase_parameters_JSON):

    def _add_phase(self, phase_dict, index=None, freeze_scale=False):
        cif_path = phase_dict["cif_path"]
        # except KeyError:
        #    cif_path = phase_dict["input_cif_path"]
        # print('scale:', phase_dict['scale']['value'])
        assert isinstance(cif_path, basestring)
        phase = rp.RietveldPhases(cif_path, phase_parameter_dict=phase_dict,
            wavelengths=RietveldServer.wavelengths,
            freeze_scale=freeze_scale
            )
        if index is None:
            RietveldServer.phase_list.append(phase)
        else:
            RietveldServer.phase_list[index] = phase
        self._bkgd_refine()

    def call_add_phase(self, phase_parameters_JSON):
        """add_phase: appends a phase to the server's phase list, given a
PhaseParameters object in json-serialized form.
        """
        try:
            # t0 = time.time()
            phase_dict = json.loads(phase_parameters_JSON)
            self._add_phase(phase_dict)
            new_phase = RietveldServer.phase_list[-1].as_dict()
            self.sendLine(json.dumps(new_phase))
            # t1 = time.time()
            # print("time taken:", str(t1-t0))
        except Exception as e:
            log.err()

    def call_update_phase(self, index, phase_parameters_JSON):
        try:
            index = int(index)
            phase_dict = json.loads(phase_parameters_JSON)
            self._add_phase(phase_dict, index=index)
            self.sendLine("True")
        except Exception as e:
            log.err()
            self.sendLine("False")


    def _append_history_entry(self):
        state = {}
        state['rietveld_data'] = \
            RietveldServer.rietveld_refinery.get_plot_data()
        state['global_state'] = rp.RietveldPhases.global_parameters.as_dict()
        state['phase_state'] = [phase.as_dict() for phase in \
            RietveldServer.phase_list]
        state['num_params_refined'] = \
            RietveldServer.rietveld_refinery.num_params
        state['goodness_of_fit'] = \
            RietveldServer.rietveld_refinery.GoF
        state['r_wp'] = RietveldServer.rietveld_refinery.R_wp
        state['time_elapsed'] = RietveldServer.rietveld_refinery.time_elapsed
        state['return_message'] = \
            RietveldServer.rietveld_refinery.result['message']

        RietveldServer.rietveld_history.append(state)

    def _calc_complete(self):
        RietveldServer.calc_flag = False

    def _refine_error(self):
        RietveldServer.err_flag = True
        log.err()

    def _update_plot_data(self):
        if rp.RietveldPhases.I is not None:
            if RietveldServer.calc_flag:
                checked_data = copy.deepcopy(
                    RietveldServer.rietveld_refinery.total_profile_state)
                if np.max(np.abs(checked_data)) < MAX_PROFILE_VALUE:
                    RietveldServer.plot_data = checked_data
            elif len(RietveldServer.phase_list) > 0:
                RietveldServer.rietveld_refinery = rr.RietveldRefinery(
                    RietveldServer.phase_list)
                RietveldServer.plot_data = \
                    RietveldServer.rietveld_refinery.total_profile_state
            else:
                RietveldServer.plot_data = np.zeros(
                    len(rp.RietveldPhases.I), dtype=float)
            # if not RietveldServer.calc_flag:
            #     self._bkgd_refine()
        else:
            RietveldServer.plot_data = np.sum(
                [phase.phase_profile() for phase in
                    RietveldServer.phase_list], axis=0)
            RietveldServer.two_thetas = rp.RietveldPhases.two_theta
        # with open('tmp.txt','w+') as f:
        #     f.write("profile length:" + str(len(profile)))
        # RietveldServer.plot_data = profile

    def _refine(self):
        RietveldServer.rietveld_refinery.minimize_all_rounds(
            callback_functions=[self._update_plot_data],
            end_of_round_callbacks=[self._append_history_entry])

    def _run(self):
        RietveldServer.calc_flag = True
        RietveldServer.err_flag = False
        d = defer.Deferred()
        d.addCallback(lambda x: self._refine())
        d.addErrback(lambda x: self._refine_error())
        d.addCallback(lambda x: self._calc_complete())
        reactor.callInThread(d.callback, None)

    def call_update_refinery_model(self, refinery_model):
        """update_refinery_model(refinery_model): returns a true or false, to
verify that the refinery model has been successfully updated.
        """
        try:
            self._set_refinery_model(json.loads(refinery_model))
            print(rp.RietveldPhases.I)
            if not RietveldServer.phase_list and rp.RietveldPhases.I is None:
                self.sendLine(str(False) + self.split_character)
            else:
                self.sendLine(str(True) + self.split_character)
        except Exception as e:
            log.err()
            self.sendLine("Error: " + e.message)


    def call_start_refine(self, refinery_model, rietveld_state):
        """start_refine(refinery_model, rietveld_state): returns a true or false
verifying that a refinement session (possibly consisting of several rounds) has
begun. The json-serialized input refinery_model, rietveld_state objects are used
to specify the refinement options and starting point of the refinement,
respectively
        """
        try:
            RietveldServer.rietveld_history = []
            self._set_refinery_model(json.loads(refinery_model))

            rs = json.loads(rietveld_state)
            RietveldServer.phase_list = []
            self._set_global_parameters(rs['global_state'])
            for phase in rs['phase_state']:
                self._add_phase(phase, freeze_scale=True)

            factr = RietveldServer.refinery_model.get('convergence_factor', 5)
            factr = 10**(-factr)
            maxiter = RietveldServer.refinery_model.get(
                'number_of_iterations', 150)
            RietveldServer.rietveld_refinery = rr.RietveldRefinery(
                RietveldServer.phase_list, factr=factr, maxiter=maxiter )
            # RietveldServer.rietveld_refinery.set_mask(
            #     ['two_theta_', 'bkgd', 'scale'])
            self._run()

            self.sendLine(str(True) + self.split_character)
        except Exception as e:
            log.err()
            self.sendLine(str(False) + self.split_character)


    def call_rounds_completed(self):
        """rounds_completed: returns the number of rounds completed to date by
the refinement engine
        """
        self.sendLine(str(len(RietveldServer.rietveld_history)) + self.split_character)

    def call_get_rietveld_state(self, round_number=-1):
        """get_rietveld_state [round_number]: returns the json-serialized
rietveld_state object corresponding to the end of the round specified. (If no
round is specified, calling this method returns the last entry found in
rietveld_history.)
        """
        round_number = int(round_number)
        max_round = len(RietveldServer.rietveld_history)
        if max_round > 0 and round_number < max_round:
            self.sendLine(json.dumps(
                RietveldServer.rietveld_history[round_number]) + self.split_character)
        else:
            self.sendLine(str(False) + self.split_character)

    def call_is_complete(self):
        """is_complete: returns either true or false, depending on whether or
not the rietveld_refinement session has completed
        """
        print("return:", not RietveldServer.calc_flag)
        self.sendLine(str(not RietveldServer.calc_flag) + self.split_character)

    def call_can_ping(self):
        """can_ping: returns True (for diagnostic purposes)
        """
        self.sendLine(str(True) + self.split_character)

    def call_get_plot_data(self):
        """get_plot_data: returns a JSON-serialized plot_data object
corresponding to the present state of the RietveldRefinery on the server
        """
        self._update_plot_data()

        try:
            reply = json.dumps(
                rr.RietveldPhases.get_plot_data(RietveldServer.plot_data,
                    two_thetas=RietveldServer.two_thetas)) + self.split_character
            # DEBUG
            if np.max(np.abs(RietveldServer.plot_data)) > MAX_PROFILE_VALUE:
                with open('reply{0}.json'.format(RietveldServer.count), 'w') as f:
                    f.write(reply)
                RietveldServer.count += 1
            self.sendLine(reply)
        except:
            log.err()

    def call_remove_phase(self, index=u'-1'):
        """remove_phase [index]: removes the phase specified by the index. If no
index is specified, the most-recently loaded phase is removed
        """
        try:
            index = int(index)
            assert index < len(self.phase_list)
            self.phase_list.pop(index)
            self.sendLine(str(True) + self.split_character)
        except Exception as e:
            log.err()
            self.sendLine(str(False) + self.split_character)


    def call_get_phase_profile(self, index=u'-1'):
        """get_phase_profile [index]: returns a json-serialized list containing
the phase profile data. If no index is specified, information for the
most-recently loaded phase is returned"""
        index = int(index)
        profile = json.dumps(
            list(RietveldServer.phase_list[index].phase_profile()),
            indent=4)
        self.sendLine(profile)

    def call_initialize(self):
        """reset: returns the rietveld_server to its initial state"""
        RietveldServer.calc_flag = False
        RietveldServer.err_flag = False
        RietveldServer.phase_list = []
        # phase_dict_list = []
        RietveldServer.refinery_model = None
        RietveldServer.rietveld_refinery = None
        RietveldServer.rietveld_history = []
        rp.RietveldPhases.two_theta = None
        rp.RietveldPhases.I = None
        rp.RietveldPhases.global_parameters.reset_x()
        self.sendLine(b'Initializing')

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
        if RietveldServer.calc_flag:
            self.sendLine(b'Analysis in progress...')
        else:
            if RietveldServer.err_flag:
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
