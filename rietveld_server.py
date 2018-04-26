from twisted.internet import reactor, protocol, defer, task
from twisted.protocols import basic
from twisted.python import log

from ast import literal_eval
import json
import sys, os

import process_roi as pr
from analyze_image import analyze_image

DEFAULT_R = 85
DEFAULT_DIR = ".\\"
DEFAULT_TWO_THETA = 30
DEFAULT_OMEGA = 15.0
DEFAULT_WAVELENGTH = 1.54

class ImageDataServer(basic.LineReceiver):
   delimiter = b'\n'
   """This is just about the simplest possible protocol"""
   img_directory = ".\\"
   img_list = []
   R = DEFAULT_R # sample-detector distance in mm
   output_dir = DEFAULT_DIR
   two_theta = DEFAULT_TWO_THETA
   wavelength = DEFAULT_WAVELENGTH
   omega = DEFAULT_OMEGA
   calc_flag = False
   err_flag = True

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
            pass
            # self.sendLine(b'Error: ' + str(e).encode("ascii"))

   def call_exit(self):
      """exit: shuts down the TCP server"""
      self.sendLine(b'Goodbye.')
      self.transport.loseConnection()
      # reactor.stop()

   def call_reset(self):
      """reset: drops the loaded list of images, roi's and sets the image
parameters to their defaults"""
      self.img_list = []
      self.R = DEFAULT_R
      self.output_dir = DEFAULT_DIR
      self.sendLine(b'''Resetting image list.
   Default R is {} mm
   Default two-theta is {} deg
   Default omega is {} deg
   Default output directory is {}'''.format(
         DEFAULT_R, DEFAULT_TWO_THETA, DEFAULT_OMEGA, DEFAULT_DIR))

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

   def call_roi_label(self, mat, hkl):
      """roi_label <material> <hkl>: read in region-of-interest labels, in
order of decreasing intensity.
         Example: roi_label YBCO (0,0,1)"""
      if mat in pr.MATERIALS:
         pass
      else:
         error_msg = 'Material not recognized. Must be one of {0}' \
            .format(str(pr.MATERIALS))
         self.sendLine(error_msg)
         raise ValueError(error_msg)
      hkl = tuple(map(int, hkl.replace('(','').replace(')','').split(',')))
      assert len(hkl) == 3
      if len(self.img_list) != 0:
         if len(self.img_list[-1]['roi_list']) != 0:
            self.sendLine(
               b"""Received:\n   Material: {}\n   hkl: {}
   Assigned to {}, ROI {}""".format(
               mat,
               str(hkl),
               self.img_list[-1]['filename'],
               str(self.img_list[-1]['roi_list'][-1]['coords'])))
            self.img_list[-1]['roi_list'][-1]['mat'].append(mat)
            self.img_list[-1]['roi_list'][-1]['hkl'].append(hkl)
         else:
            self.sendLine(
               b"""Error: no region-of-interest specified for {}""".format(
                  sef.img_list[-1]['filename']))
      else:
         self.sendLine(
            b"""Error: no file names specified.""")

   def call_roi(self, coords):
      """roi <coords>: set region-of-interest coordinates.
         Example: roi (320,370,140,170)"""
      coords = tuple(map(int, coords.replace('(','').replace(')','').split(',')))
      assert len(coords) == 4
      if len(self.img_list) is not 0:
         self.sendLine(
            b"""Received:\n   ROI Coordinates: {}
   Assigned to {}""".format(str(coords), self.img_list[-1]['filename']))
         self.img_list[-1]['roi_list'].append(
            {
            'mat': [],
            'hkl': [],
            'coords': coords,
            })
      else:
         self.sendLine(
            b"""Error: no file names specified.""")

   def call_is_calculating(self):
      if self.calc_flag:
         self.sendLine(b'Image analysis in progress...')
      else:
         if self.err_flag:
            self.sendLine(b'Error in running image analysis')
         else:
            self.sendLine(b'Analysis complete')

   def _calc_complete(self):
      self.calc_flag = False

   def _analyze_image(self, index):
      self.img_list[index] = analyze_image(self.img_list[index])

   def _analyze_error(self, f):
      self.err_flag = True
      log.err()

   def call_calc_last(self):
      self.calc_flag = True
      self.err_flag = False
      self.sendLine(b'Starting calculation...')
      d = defer.Deferred()
      d.addCallback(self._analyze_image)
      d.addErrback(self._analyze_error)
      d.addCallback(lambda x: self._calc_complete())
      reactor.callInThread(d.callback, -1)

   def call_img(self, filename):
      """img <filename>: add a new file to the server's records
      """
      filename = str(filename)
      img_dict = {}
      img_dict['filename'] = filename
      img_dict['R'] = self.R
      img_dict['two_theta'] = self.two_theta
      img_dict['omega'] = self.omega
      img_dict['wavelength'] = self.wavelength
      img_dict['img_dir'] = self.img_directory
      img_dict['roi_list'] = []
      self.img_list.append(img_dict)
      self.sendLine(
         b"""\
New image: {}
   R: {} mm
   two-theta: {} deg
   omega: {} deg
   wavelength: {} Angstroms
""".format(
            img_dict['img_dir'] + img_dict['filename'],
            img_dict['R'],
            img_dict['two_theta'],
            img_dict['omega'],
            img_dict['wavelength'],
            ))

   def call_img_dir(self, path):
      """img_dir <path>: New images will be discovered in the path specified.
      """
      self.img_directory = str(os.path.abspath(path)) + "\\"
      self.sendLine(b"New image directory: {}".format(str(self.img_directory)))

   def call_update_R(self, R):
      """New images will be assigned the sample-detector distance (in mm) specified.
            Example: update_R 100"""
      self.R = float(R)
      self.sendLine(b"New sample-detector distance: {}".format(str(R)))

   def call_update_two_theta(self, two_theta):
      """New images will be assigned the two-theta angle (in degrees) specified.
            Example: update_two_theta 30"""
      self.two_theta = float(two_theta)
      self.sendLine(b"New two-theta value: {}".format(str(two_theta)))

   def call_update_wavelength(self, wavelength):
      """New images will be assigned the wavelength (in Angstroms) specified.
            Example: update_wavelength 1.54"""
      self.wavelength = float(wavelength)
      self.sendLine(b"New wavelength value: {}".format(str(wavelength)))

   def call_update_omega(self, omega):
      """New images will be assigned the omega angle (in degrees) specified.
            Example: update_omega 1.54"""
      self.omega = float(omega)
      self.sendLine(b"New omega value: {}".format(str(omega)))

   def call_output_dir(self, path):
      """output_dir <path>: server state and results will be stored in the
directory specified.
      """
      self.output_dir = str(path)
      self.sendLine(b"New output directory: {}".format(str(self.output_dir)))

   def call_writeJSON(self):
      """Dumps the current server state to state.json in the output directory."""
      p = json.dumps(self.img_list, indent=4)
      with open(self.output_dir+'state.json', 'w') as outfile:
         outfile.write(p)
      self.sendLine(b"Output written to {}state.json".format(self.output_dir))

   def call_loadJSON(self):
      """Loads the state of the image list, as stored in state.json
in the current output directory."""
      self.img_list = json.load(open(self.output_dir+'state.json'))
      # print self.img_list
      self.sendLine(b"Loaded {} images from state.json".format(
         len(self.img_list)))

   def connectionLost(self, reason):
      pass
      # reactor.stop()

def main():
   factory = protocol.ServerFactory()
   factory.protocol = ImageDataServer
   log.startLogging(sys.stdout)
   reactor.listenTCP(8006, factory)
   reactor.run()

# this only runs if the module was *not* imported
if __name__ == '__main__':
   main()