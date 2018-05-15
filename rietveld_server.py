from twisted.internet import reactor, protocol, defer, task
from twisted.protocols import basic
from twisted.python import log

import json
import sys, os

import src.rietveld_phases as rp
import src.rietveld_refinery as rr

class RietveldServer(basic.LineReceiver):
   delimiter = b'\n'
   split_character = ';'
   calc_flag = False
   err_flag = False
   phase_list = []
   phase_dict_list = []
   refinery_model = None
   rietveld_refinery = None

   def connectionMade(self):
      # pass
      self.sendLine(b'Connected. Type <help> [command] for info.')
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

   def call_load_profile(self, refinery_model_string, global_parameters):
      try:
         assert isinstance(refinery_model_string, unicode)
         self.refinery_model = json.loads(refinery_model_string)
         profile_path = self.refinery_model['input_data_path']
         min_two_theta = self.refinery_model['two_theta_roi_window'][0]
         max_two_theta = self.refinery_model['two_theta_roi_window'][1]
         rp.RietveldPhases.set_profile(profile_path,
            min_two_theta=min_two_theta,
            max_two_theta=max_two_theta,
            )
         # profile_filename = os.path.split(profile_path)[1]
         # self.sendLine(b'Loaded the profile {0}'.format(profile_filename))
         self.rietveld_refinery = rr.RietveldRefinery(self.phase_list,
            bkgd_refine=True)
         self.rietveld_refinery.minimize()
         profile = self.rietveld_refinery.total_profile()# reply = {}
         # reply['two_thetas'] = []
         # reply['intensities'] = list(self.rietveld_refinery.total_profile())
         # reply['errors'] = []

         def get_plot_data(nparray, two_thetas=[], errors=[]):
            d = {}
            d['two_thetas'] = two_thetas
            d['errors'] = errors
            d['intensities'] = list(nparray)
            return d
         reply = {}
         reply['input_data'] = get_plot_data(rp.RietveldPhases.I, 
            two_thetas=list(rp.RietveldPhases.two_theta), 
            errors=list(rp.RietveldPhases.sigma))
         reply['profile_data'] = get_plot_data(profile)
         reply['differences'] = get_plot_data(rp.RietveldPhases.I - profile)
         # reply['two_thetas'] = []
         # reply['intensities'] = list(self.rietveld_refinery.total_profile())
         # reply['errors'] = []
         self.sendLine(json.dumps(reply, indent=4))
      except:
         log.err()

   def call_add_phase(self, json_string):
      """add_phase: appends a phase to the server's phase list, given a
PhaseParameters object in json-serialized form.
      """
      try:
         phase_dict = json.loads(json_string)
         try:
            cif_path = phase_dict["cif_path"]
         except KeyError:
            cif_path = phase_dict["input_cif_path"]
         assert type(cif_path) == unicode
         self.phase_list.append(rp.RietveldPhases(cif_path))
         self.phase_dict_list.append(phase_dict)
         # print json.dumps(self.phase_list[-1].get_phase_info())
         # self.sendLine(json.dumps(self.phase_list[-1].get_phase_info()))
         self.sendLine(b'Added {0} to the server\'s phase list'.format(
            self.phase_list[-1].phase_settings["chemical_name"]))
      except:
         log.err()

   def call_get_phase_info(self, index=u'-1'):
      """get_phase_info [index]: returns a json-serialized dictionary containing
relevant phase information. If no index is specified, information for the
most-recently loaded phase is returned"""
      index = int(index)
      self.phase_list[index].update_phase_info(
         self.phase_dict_list[index])
      # info = json.dumps(self.phase_list[index].get_phase_info(), indent=4)
      info = json.dumps(self.phase_dict_list[index], indent=4)
      self.sendLine(info)

   def call_get_phase_profile(self, index=u'-1'):
      """get_phase_profile [index]: returns a json-serialized list containing
the phase profile data. If no index is specified, information for the
most-recently loaded phase is returned"""
      index = int(index)
      profile = json.dumps(list(self.phase_list[index].phase_profile()),
         indent=4)
      self.sendLine(profile)

   def call_reset(self):
      """reset: drops the loaded list of phases"""
      self.phase_list = []
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

   def _calc_complete(self):
      self.calc_flag = False

   def _refine_error(self, f):
      self.err_flag = True
      log.err()

   def _refine(self, f):
      pass

   def call_calc(self):
      self.calc_flag = True
      self.err_flag = False
      self.sendLine(b'Starting calculation...')
      d = defer.Deferred()
      d.addCallback(self._refine)
      d.addErrback(self._refine_error)
      d.addCallback(lambda x: self._calc_complete())
      reactor.callInThread(d.callback, -1)

   def call_writeJSON(self):
      """Dumps the current server state to state.json in the output directory."""
      pass
      # p = json.dumps(self.img_list, indent=4)
      # with open(self.output_dir+'state.json', 'w') as outfile:
      #    outfile.write(p)
      # self.sendLine(b"Output written to {}state.json".format(self.output_dir))

   def call_loadJSON(self):
      """Loads the state of the image list, as stored in state.json
in the current output directory."""
      # self.img_list = json.load(open(self.output_dir+'state.json'))
      # # print self.img_list
      self.sendLine(b"Loaded from state.json")

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