from twisted.internet import reactor, protocol, defer, task
from twisted.protocols import basic
from twisted.python import log

import json
import sys, os

import src.rietveld_phases as rp
import src.rietveld_refinery as rr

class RietveldServer(basic.LineReceiver):
   delimiter = b'\n'
   calc_flag = False
   err_flag = False
   phase_list = []

   def connectionMade(self):
      self.sendLine(b'Connected. Type <help> [command] for info.')
      # self.setLineMode()
      # print self.__dict__
      # self.sendLine(b"Web checker console. Type 'help' for help.")

   def lineReceived(self, data):
      if not data: return
      data = data.decode("ascii")

      # Parse the command
      commandParts = data.split()
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

   def call_load_profile(self, profile_path):
      try:
         assert type(profile_path) == unicode
         rp.RietveldPhases.set_profile(profile_path)
         profile_filename = os.path.split(profile_path)[1]
         self.sendLine(b'Loaded the profile {0}'.format(profile_filename))
      except:
         log.err()

   def call_add_phase(self, cif_path):
      """add_phase: appends a phase to the server's phase list"""
      try:
         assert type(cif_path) == unicode
         self.phase_list.append(rp.RietveldPhases(cif_path))
         self.sendLine(b'Added {0} to the server\'s phase list'.format(
            self.phase_list[-1].phase_settings["chemical_name"]))
      except:
         log.err()

   def call_get_phase_info(self, index=u'-1'):
      """get_phase_info [index]: returns a json-serialized dictionary containing
relevant phase information. If no index is specified, information for the
most-recently loaded phase is returned"""
      index = int(index)
      info = json.dumps(self.phase_list[index].get_phase_info(), indent=4)
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