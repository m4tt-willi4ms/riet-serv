from __future__ import division
# from scitbx.array_family import flex
# import scitbx.math
# from libtbx.test_utils import approx_equal, eps_eq, Exception_expected
# import sys
import numpy as np
import time
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import \
   FigureCanvasTkAgg, NavigationToolbar2TkAgg
from  matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib import style, use
# use("TkAgg")
style.use("ggplot")

import matplotlib.animation as animation
import ttk
import tkFont
from multiprocessing import Process
from threading import Thread

import json,codecs

import sys, os
# sys.path.append(os.path.abspath("../.."))

from src.RietveldPhases import RietveldPhases
from src.RietveldRefinery import RietveldRefinery

from scitbx import lbfgsb
from cctbx.eltbx import wavelengths
from cctbx.sgtbx import lattice_symmetry

import Tkinter as tk
import tkFileDialog

LARGE_FONT = ("Verdana", 11)


fig = Figure(dpi=100)
# fig, ax = plt.subplots(nrows=3,ncols=1,sharex=True,dpi=100)
# xdata, ydata = [], []
# ln, = plt.plot([], [], 'ro', animated=True)
# fig.suptitle(plottitle)
canvas = None
ROI_mask = None
ROI_center = 32
ROI_delta_theta = 5
x_mask = np.zeros(0,dtype=bool)

subplot1 = fig.add_subplot(311) #plt.subplot(3,1,1)
      # plt.scatter(self.two_theta,self.y,label='Data',s=1, color='red')

      # self.current_profile, = self.subplot1.plot(self.two_theta,
      #    self.TotalProfile_state, label=r'$I_{\rm calc}$')
# subplot1.legend(bbox_to_anchor=(.8,.7))
# subplot1.ylabel(r"$I$")

subplot2 = fig.add_subplot(312) #plt.subplot(3,1,2)
# self.pltmask = np.abs(self.two_theta - two_theta_roi) \
#          < delta_theta
#       plt.scatter(self.two_theta[self.pltmask],self.y[self.pltmask],
#          label='Data',s=2, color='red')

# self.current_profile_masked, = subplot2.plot(self.two_theta[self.pltmask],
#    self.TotalProfile_state[self.pltmask], label=r'$I_{\rm calc}$')
      # plt.legend(bbox_to_anchor=(.8,.7))
# plt.ylabel(r"$I$")

subplot3 = fig.add_subplot(313) #plt.subplot(3,1,3)
      # self.residuals, = subplot3.plot(self.two_theta,
      #    self.Weighted_Squared_Errors(),'bo',ms=2)
# plt.ylabel(r"$\frac{1}{I} \, (I-I_{\rm calc})^2$")
# plt.xlabel(r'$2\,\theta$')

two_thetas = []
ys = []

Rt = []
RR = None

x_default = np.empty(0,dtype=RietveldPhases.custom_dtype)

max_refinement_rounds = 5

CU_wavelength = wavelengths.characteristic("CU").as_angstrom()

global_input_string = """\
Bkgd:          3
two_theta_0       0.      -0.5  0.5
"""

default_input_string = """\
U              0.00    0     0.1
V              -0.00   -0.1   0
W              0.01   0.0001     1
Amplitude         0.1 0      inf
eta:           3
unit_cell_a    0.01
unit_cell_b    0.01
unit_cell_c    0.01
unit_cell_alpha   0.005
unit_cell_beta    0.005
unit_cell_gamma   0.005
"""

minimizer_input_string = """\
factr       1e10
iprint      -1
maxiter     150
m           10
pgtol       1e-5
epsilon     1e-13
"""

def animate(i):
   profile_ys = json.load(
      codecs.open("current_profile.json", 'r', encoding='utf-8'))
   for x,y in zip(two_thetas,profile_ys):
      print x,y
   subplot1.clear()
   print "here"
   subplot1.axes.lines[0].set_ydata(RR.TotalProfile_state)
   # subplot2.axes.lines[0].set_ydata(RR.TotalProfile_state[ROI_mask])
   # subplot3.axes.lines[0].set_ydata(RR.Weighted_Squared_Errors_state)
   return subplot1.axes.lines,

# json.dump(self.TotalProfile().tolist(), 
#                codecs.open("current_profile.json", 'w', encoding='utf-8'), 
#                separators=(',', ':'), 
#                # sort_keys=True, 
#                indent=4)
      # self.fig.canvas.draw()

class RietveldGUI(tk.Tk):
   def __init__(self, *args, **kwargs):
      tk.Tk.__init__(self, *args, **kwargs)

      tk.Tk.iconbitmap(self, default=u"doc//ProtoLogo_R.ico")
      tk.Tk.wm_title(self, "Rietveld Refinement")
      tk.Tk.wm_geometry(self,'1200x600')

      self.container = tk.Frame(self)
      self.container.pack(side="top", fill="both", expand=True)

      self.container.grid_rowconfigure(0, weight=1)
      self.container.grid_columnconfigure(0, weight=1)

      self.menu = tk.Menu(self)
      self.config(menu=self.menu)
      
      self.FileMenu = tk.Menu(self.menu,tearoff=False)

      self.menu.add_cascade(label="File", menu=self.FileMenu)
      self.FileMenu.add_command(label="Load Profile...",command=self.getProfile)
      self.FileMenu.add_command(label="Load .cifs...", command=self.getCifs)
      self.FileMenu.add_command(label="Exit", command=self.exit)

      self.EditMenu = tk.Menu(self.menu,tearoff=False)
      self.menu.add_cascade(label="Edit", menu=self.EditMenu)
      self.EditMenu.add_command(label="Refinement Options")

      self.HelpMenu = tk.Menu(self.menu,tearoff=False)
      self.menu.add_cascade(label="Help", menu=self.HelpMenu)
      self.HelpMenu.add_command(label="Tutorial")
      self.HelpMenu.add_command(label="About")

      s = ttk.Style()
      s.theme_use('clam')
      # s.configure('ButtonText.TButton', font=('Verdana', 11))

      font_dict = {'family': 'Helvetica', 'size': 10}
      # k=v for (k,v) in font_dict.items()
      # dict_2 = {k: v for (k,v) in font_dict.items()}
      DEFAULT_FONT = tkFont.Font(**font_dict)
      tkFont.nametofont("TkDefaultFont").configure(**font_dict)
      tkFont.nametofont('TkHeadingFont').configure(**font_dict)
      tkFont.nametofont('TkMenuFont').configure(**font_dict)
      # default_font.configure(**font_dict)
      self.option_add("*Font",tkFont.nametofont("TkDefaultFont"))

      # s.theme_create( "yummy", parent="alt", settings={
      #   "TNotebook": {"configure": {"tabmargins": [2, 5, 2, 0] } },
      #   "TNotebook.Tab": {
      #       "configure": {"padding": [5, 1], "background": "grey" },
      #       "map":       {"background": [("selected", "blue")],
      #                     "expand": [("selected", [1, 1, 1, 0])] } } } )

      # s.theme_use("yummy")
      # s.theme_create( "MyStyle", parent="alt", settings={
      #   "TNotebook": {"configure": {"tabmargins": [2, 5, 2, 0] } },
      #   "TNotebook.Tab": {"configure": {"padding": [100, 100] },}})
      # s.theme_use("MyStyle")

      s.layout("Tab",
      [('Notebook.tab', {'sticky': 'nswe', 'children':
          [('Notebook.padding', {'side': 'top', 'sticky': 'nswe', 'children':
              #[('Notebook.focus', {'side': 'top', 'sticky': 'nswe', 'children':
                  [('Notebook.label', {'side': 'top', 'sticky': ''})],
              #})],
          })],
      })]
      )

      s.configure("Button",borderwidth=0)

      self.numPhases=0

      self.plotframe = PlotFrame(self.container,self,pady=10,padx=10)
      self.plotframe.grid(row=0,column=1,sticky='ew')

      # temp. to allow for auto-loading of profile
      self.getProfile()
      self.getCifs()

   def getCifs(self):
      # self.filePaths = tkFileDialog.askopenfilenames(
      #    initialdir = "./data/cifs")
      self.filePaths = [r".\data\cifs\Cement\1540705-Alite.cif"]
      # iconEntry.insert(0,fileName)
      global Rt
      for i,filePath in enumerate(self.filePaths):
         cif_file_name = os.path.split(filePath)[1]
         Rt.append(RietveldPhases(filePath, d_min,d_max, 
            input_string_or_file_name = default_input_string, I_max = y_max, 
            delta_theta=1.5,Intensity_Cutoff = 0.005))
         if self.numPhases == 0:
            self.paramframe.nb.add(RefinementParameterSet(self.paramframe.nb,
               self.paramframe,index=0,change_all=True),text="Phases: All")
         self.paramframe.nb.add(
            RefinementParameterSet(self.paramframe.nb,self.paramframe,index=i), 
            text=str(i+1))
         self.numPhases += 1
      # RietveldPhases.empty_x()
      # RietveldPhases.global_params_from_string(global_input_string,
      #    two_thetas,ys)
      global RR
      RR = RietveldRefinery(Rt,minimizer_input_string,
      store_intermediate_state=True, show_plots=False)

      updateplotprofile()


   def getProfile(self):
      # self.fileName = tkFileDialog.askopenfilename(
      #    initialdir = "./data/profiles")
      self.fileName = r".\\data\\profiles\\cement_15_03_11_0028.xye"
      self.winfo_toplevel().title("Rietveld Refinement (" + 
         os.path.split(self.fileName)[1]+")")

      global two_thetas, ys
      two_thetas = []
      ys = []
      with open(self.fileName) as file:
         for line in file.readlines()[1:]:
            two_thetatmp, ytmp, ztmp = line.split()
            # two_thetatmp, ytmp = line.split()
            if float(two_thetatmp) > 20:
               two_thetas.append(float(two_thetatmp))
               ys.append(float(ztmp)**2)

      global d_min, d_max, y_max
      d_min = CU_wavelength/2/np.sin(np.pi/360*two_thetas[-1])
      d_max = CU_wavelength/2/np.sin(np.pi/360*two_thetas[0])
      y_max = np.amax(ys)#/len(cifs)

      two_thetas = np.array(two_thetas)
      ys = np.array(ys)
      global ROI_mask
      ROI_mask = np.abs(two_thetas - ROI_center) < ROI_delta_theta

      RietveldPhases.global_params_from_string(global_input_string,
         two_thetas,ys)

      # global defaultphase
      # defaultphase = RietveldPhases('.//data//cifs//1000032.cif',
      #    d_min,d_max, input_string_or_file_name=default_input_string)

      self.paramframe = LoadFrame(self.container,self,padx=10)#,pady=10)
      self.paramframe.grid(row=0,column=0) 

      updateplotdata()

   def exit(self):
      self.destroy()

def updateplotdata():
   subplot1.clear()
   subplot1.scatter(two_thetas,ys,label='Data',s=1, color='red')
   subplot1.legend(bbox_to_anchor=(.8,.7))
   subplot1.set_ylabel(r"$I$")
   subplot2.clear()
   subplot2.scatter(two_thetas[ROI_mask],ys[ROI_mask],
      label='Data',s=1, color='red')
   subplot2.set_ylabel(r"$I$")
   canvas.show()

class VarLabelEntry(tk.Frame):
   def __init__(self,parent,text,x_label,index,*args, **kwargs):
      tk.Frame.__init__(self,parent)
      self.x_label = x_label
      self.index = index
      self.change_all = parent.change_all
      self.is_l_limits = self.x_label == "l_limits"
      self.is_u_limits = self.x_label == "u_limits"
      self.is_values = self.x_label == "values"
      self.passable = ["","-",".","-."]

      # self.value = tk.StringVar()
      # self.value.set(str(RietveldPhases.x[x_label][index]))
      self.label = tk.Label(self,text=text)

      # valid percent substitutions (from the Tk entry man page)
      # note: you only have to register the ones you need; this
      # example registers them all for illustrative purposes
      #
      # %d = Type of action (1=insert, 0=delete, -1 for others)
      # %i = index of char string to be inserted/deleted, or -1
      # %P = value of the entry if the edit is allowed
      # %s = value of entry prior to editing
      # %S = the text string being inserted or deleted, if any
      # %v = the type of validation that is currently set
      # %V = the type of validation that triggered the callback
      #      (key, focusin, focusout, forced)
      # %W = the tk name of the widget

      vcmd = (self.register(self.onValidate),
         '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
      self.entry = tk.Entry(self,width=8,validate="key",validatecommand=vcmd)
      self.entry.insert(0,str(RietveldPhases.x[x_label][index]))
      self.label.grid(row=0,column=0)
      self.entry.grid(row=0,column=1)

   def onValidate(self, d, i, P, s, S, v, V, W):
      try:
         if all(P != x for x in self.passable):
            val = float(P)
            if self.is_l_limits:
               if val <= RietveldPhases.x['u_limits'][self.index]:
                  RietveldPhases.x[self.x_label][self.index] = val
               else: 
                  raise ValueError
            elif self.is_u_limits:
               if val >= RietveldPhases.x['l_limits'][self.index]:
                  RietveldPhases.x[self.x_label][self.index] = val
               else: 
                  raise ValueError
            elif self.is_values:
               RietveldPhases.x[self.x_label][self.index] = val
         else:
            RietveldPhases.x[self.x_label][self.index] = 0
         global RR

      except(ValueError):
         return False
      if self.change_all:
         for Rp in Rt:
            RietveldPhases.x[self.x_label]
      if RR is not None:
         RR.Update_state()
         updateplotprofile()
      return True

class RoundsBoxes(tk.Frame):
   def __init__(self,parent,controller,default_round_start=1):
      tk.Frame.__init__(self, parent)
      self.RoundsLabel = tk.Label(self, text = "Rounds: ")
      self.RoundsLabel.grid(row=0,column=0,sticky='e')

      self.roundcheckbuttons = []
      self.roundstates = []

      for i in xrange(max_refinement_rounds):
         self.roundstates.append(tk.IntVar())
         if i >= default_round_start-1:
            self.roundstates[i].set(1)
         self.roundcheckbuttons.append(tk.Checkbutton(self,
            variable = self.roundstates[i], text=str(i+1)))
         self.roundcheckbuttons[i].grid(row=0,column=i+1)

class RefinementParameterControl(tk.Frame):
   def __init__(self, parent, controller, index, text="", 
      default_round_start=1,change_all=False,*args, **kwargs):
      tk.Frame.__init__(self, parent)
      self.parent = parent
      self.change_all = change_all
      self.index = index

      self.state = tk.IntVar()
      self.checkbutton = tk.Checkbutton(self, command=self.checkbutton_clicked,
         variable = self.state, text=text)
      self.checkbutton.grid(row=0,column=0,sticky='w')

      self.rounds = RoundsBoxes(self,parent,default_round_start)
      self.rounds.grid(row=0,column=1,sticky='w')

      # self.initial = VarLabelEntry(self,'Start at:', 'values', index)
      # self.initial.grid(row=0,column=1,sticky='e')
      
      # self.l_limit = VarLabelEntry(self,'Lower limit:', 'l_limits', index)
      # self.l_limit.grid(row=1,column=0,sticky='w')

      # self.u_limit = VarLabelEntry(self,'Upper limit:', 'u_limits', index)
      # self.u_limit.grid(row=1,column=1,sticky='e')

      # self.round_dropdownlist = Dropdown_Int_List(self, parent,
      #    text="Round:", min_int=1, max_int=max_refinement_rounds, 
      #    default_int=default_round)
      # self.round_dropdownlist.grid(row=0,column=2,rowspan=2,sticky='e')

      

      self.grid_columnconfigure(0,minsize=150)
      # self.grid_columnconfigure(1,minsize=120)
      # self.grid_columnconfigure(2,minsize=110)

      self.checkbutton_clicked()

   def checkbutton_clicked(self):
      if self.state.get() == 1:
         self.rounds.grid()
         # self.initial.grid()
         # self.l_limit.grid()
         # self.u_limit.grid()
         # self.round_dropdownlist.grid()
      if self.state.get() == 0:
         self.rounds.grid_remove()
         # self.initial.grid_remove()
         # self.l_limit.grid_remove()
         # self.u_limit.grid_remove()
         # self.round_dropdownlist.grid_remove()
      # if self.change_all:
      #    global Rt
      #    for phase in Rt:
      #       RietveldPhases.x['values'][index] = 0

class LatticeParameterControl(tk.Frame):
   def __init__(self,parent,controller,index):
      tk.Frame.__init__(self,parent)
      self.index = index
      self.crystal_system = Rt[index].structure.space_group().crystal_system()



class RefinementParameterPolynomControl(tk.Frame):
   def __init__(self, parent, controller, 
      text="", default_order=2, default_round_start=1,*args, **kwargs):
      tk.Frame.__init__(self, parent)
      self.state = tk.IntVar()
      self.parent = parent
      self.text = text

      self.checkbutton = tk.Checkbutton(self, command=self.checkbutton_clicked,
         text=text, variable = self.state)
      self.checkbutton.grid(row=0,column=0,sticky='w')

      self.order_dropdownlist = Dropdown_Int_List(self, parent,
         text="Order:", min_int=1, max_int=6, default_int=default_order)
      self.order_dropdownlist.grid(row=0,column=1)

      self.rounds = RoundsBoxes(self,parent,default_round_start)
      self.rounds.grid(row=0,column=2)

      # self.round_dropdownlist = Dropdown_Int_List(self, parent,
      #    text="Round:", min_int=1, max_int=max_refinement_rounds, 
      #    default_int=default_round)
      # self.round_dropdownlist.grid(row=0,column=2)

      self.grid_columnconfigure(0,minsize=90)
      # self.grid_columnconfigure(1,minsize=90)
      # self.grid_columnconfigure(2,minsize=110)

      self.checkbutton_clicked()

   def checkbutton_clicked(self):
      global RR
      if self.state.get() == 1:
         self.order_dropdownlist.grid()
         self.rounds.grid()

         if self.text == "Background":
            RR = RietveldRefinery(Rt,minimizer_input_string, \
               use_bkgd_mask=True,
               store_intermediate_state=True, show_plots=False)
            RR.minimize_Bkgd()
            updateplotprofile()

      if self.state.get() == 0:
         self.order_dropdownlist.grid_remove()
         self.rounds.grid_remove()

         if self.text == "Background":
            if len(subplot1.axes.lines) is not 0:
               for x in (subplot1,subplot2,subplot3):
                  for line in x.axes.lines:
                     line.remove()
            canvas.show()

def startthread(func):
   # print "here"
   t = Thread(target=func)
   t.start()
   return t

def updateplotprofile():
   if len(subplot1.axes.lines) is not 0:
      for x in (subplot1,subplot2,subplot3):
         for line in x.axes.lines:
            line.remove()
   subplot1.plot(two_thetas,RR.TotalProfile_state,
      label=r'$I_{\rm calc}$',color='blue')
   subplot1.legend(bbox_to_anchor=(.8,.7))
   subplot2.plot(two_thetas[ROI_mask],RR.TotalProfile_state[ROI_mask],
      color='blue')
   if RR is not None:
      subplot3.plot(two_thetas,RR.Weighted_Squared_Errors_state,
         label=r'$\frac{1}{I} \, (I-I_{\rm calc})^2$',color='green')
   subplot3.set_xlabel(r'$2\,\theta$')
   subplot3.set_ylabel(r"$\frac{(I-I_{\rm calc})^2}{I}$")
   print Rt[0].structure.space_group().crystal_system()
   print RietveldPhases.x['values'][Rt[0].unit_cell_indices[0]]
   Rt[0].update_unit_cell()
   print Rt[0].unit_cell
   print Rt[0].structure.space_group().average_unit_cell(Rt[0].unit_cell)
   # print lattice_symmetry.group(Rt[0].unit_cell)
   # print tmp.lattice_group_info()
   canvas.draw()

class Dropdown_Int_List(tk.Frame):
   def __init__(self, parent, controller, text="", default_int=2, min_int=0,
         max_int=5,*args, **kwargs):
      tk.Frame.__init__(self, parent)

      self.order_label = tk.Label(self,text=text)
      self.order_label.grid(row=0,column=0)

      self.order = tk.StringVar(self)
      self.orderMenu = ttk.Combobox(self, textvariable=self.order,
         state='readonly',width =2)
      order_values = ()
      default_index = 0
      for (index,i) in enumerate(xrange(min_int,max_int+1)):
         order_values += (str(i),)
         if i == default_int:
            default_index = index
      self.orderMenu['values'] = order_values
      self.orderMenu.current(default_index)
      self.orderMenu.bind("<<ComboboxSelected>>", self.order_selected)
      self.orderMenu.grid(row=0,column=1)

      self.grid_columnconfigure(0,minsize=50)
      self.grid_columnconfigure(1,minsize=50)

   def order_selected(self,tmp):
      print self.order.get()
      self.orderMenu.selection_clear()


class RefinementParameterSet(tk.Frame):
   def __init__(self, parent, controller, index=None, change_all=False,
       *args, **kwargs):
      tk.Frame.__init__(self, parent, padx=10,pady=10)

      if index is not None:
         phase = Rt[index]
      else:
         phase = Rt[0]

      RefinementParameterControl(self,parent,
         phase.Amplitude_index,text="Amplitude",change_all=change_all
         ,default_round_start=2).grid(row=0,column=0,sticky='w')

      RefinementParameterControl(self,parent,
         phase.W_index,text="Caglioti W",default_round_start=2) \
         .grid(row=1,column=0,sticky='w')

      RefinementParameterPolynomControl(self,parent,
         text="eta",default_order=1,default_round_start=3).grid(row=2,column=0,
         sticky='w')

      RefinementParameterControl(self,parent,
         phase.V_index,text="Caglioti V",default_round_start=2) \
         .grid(row=3,column=0,sticky='w')

      RefinementParameterControl(self,parent,
         phase.U_index,text="Caglioti U",default_round_start=2) \
         .grid(row=4,column=0,sticky='w')

      RefinementParameterControl(self,parent,
         phase.unit_cell_indices[3],text="Lattice Parameters",
         default_round_start=2).grid(row=5,column=0,sticky='w')

class LoadFrame(tk.Frame):
   def __init__(self, parent,controller,*args,**kwargs):
      # controller.minsize(width=400,height=400)
      tk.Frame.__init__(self,parent,*args,**kwargs)

      self.controller = controller
      self.globalnb = ttk.Notebook(self,height=100,width=450)
      self.globalFrame = tk.Frame(self.globalnb,padx=10,pady=10)
      #, width=100, height=100)
      # self.globalLabelFrame.p#grid(row=0,column=0)

      RefinementParameterPolynomControl(self.globalFrame,self,
         text="Background",default_order=2,default_round_start=1).grid(
         row=0,column=0,sticky='w')
      RefinementParameterControl(self.globalFrame,self,
         RietveldPhases.two_theta_0_index,text="Two theta offset") \
         .grid(row=1,column=0,sticky='w')
      
      self.globalnb.add(self.globalFrame,text="Global Parameters")
      self.globalnb.grid(row=0,column=0,columnspan=2,padx=10,pady=10)

      self.nb = ttk.Notebook(self,
         height=300,width=450)
      self.nb.enable_traversal()
      # 
      self.nb.grid(row=1,column=0,columnspan=2,padx=10,pady=10)
      # self.nb.grid(row=2,column=1,columnspan=2)

      self.RefineButton = ttk.Button(self, text='Refine', 
         command=self.refine, takefocus=False)
      self.RefineButton.grid(row=2,column=0, padx=10, pady=10)

      self.CancelButton = ttk.Button(self, text='Cancel', 
         command=self.cancel, takefocus=False)
      self.CancelButton.grid(row=2,column=1, padx=10, pady=10)

   def refine(self):
      print RietveldPhases.x['values'][RietveldPhases.two_theta_0_index]
      global RR
      t = startthread(RR.minimize_Amplitude_Bkgd_Offset)
      fig.suptitle("Progress: minimize_Amplitude_Bkgd_Offset")
      while t.isAlive():
         updateplotprofile()
         time.sleep(.5)
      fig.suptitle("Completed: minimize_Amplitude_Bkgd_Offset")
      updateplotprofile()
      t.join()
      RR.display_parameters(RR.minimize_Amplitude_Bkgd_Offset)
      RR.display_stats(RR.minimize_Amplitude_Bkgd_Offset)

   def cancel(self):
      self.controller.numPhases = 0
      global Rt, RR
      Rt = []
      RR = None
      RietveldPhases.empty_x()
      RietveldPhases.global_params_from_string(global_input_string,
         two_thetas,ys)
      for tab in self.nb.tabs():
         self.nb.hide(tab)
         self.nb.forget(tab)
      if len(subplot1.axes.lines) is not 0:
         for x in (subplot1,subplot2,subplot3):
            for line in x.axes.lines:
               line.remove()
         fig.suptitle("")
         subplot1.axes.legend_.remove()
         canvas.draw()

class PlotFrame(tk.Frame):
   def __init__(self, parent,controller,*args,**kwargs):
      tk.Frame.__init__(self,parent,*args,**kwargs)

      global canvas 
      canvas = FigureCanvasTkAgg(fig, master=self)
      canvas.get_tk_widget().pack(side=tk.TOP)
      # canvas.get_tk_widget().grid(row=0,column=0,sticky='n')

      toolbar = NavigationToolbar2TkAgg(canvas,self)
      toolbar.update()
      # canvas._tkcanvas.pack()
      # canvas._tkcanvas.grid()#row=0,column=0)#,sticky='n',pady=10)
      # canvas.show()


if __name__ == "__main__":
   root = RietveldGUI()
   root.mainloop()

input_strings = ["""\
U              0.00    0     0.1
V              -0.00   -0.1   0
W              0.01   0.0001     1
Amplitude         0.1 0      inf
eta:           3
unit_cell_a    0.01
unit_cell_b    0.01
unit_cell_c    0.01
unit_cell_alpha   0.005
unit_cell_beta    0.005
unit_cell_gamma   0.005
""",
"""\
U              0.00    0     0.1
V              -0.00   -0.1   0
W              0.01   0.0001     1
Amplitude         0.000001 0      inf
eta:           2
""",
"""\
U              0.00    0     0.1
V              -0.00   -0.1   0
W              0.01   0.0001     1
Amplitude         0.1 0      inf
eta:           2
""",
"""\
U              0.00    0     0.1
V              -0.00   -0.1   0
W              0.01   0.0001     1
Amplitude         0.1 0      inf
eta:           2
""",
"""\
U              0.00    0     0.1
V              -0.00   -0.1   0
W              0.01   0.0001     1
Amplitude         0.1 0      inf
eta:           2
""",
"""\
U              0.00    0     0.1
V              -0.00   -0.1   0
W              0.01   0     1
Amplitude         0.1 0      inf
eta:           2
""",
"""\
U              0.00    0     0.1
V              -0.00   -0.1   0
W              0.01   0.0001     1
Amplitude         0.1 0      inf
eta:           2
""",
"""\
U              0.00    0     0.1
V              -0.00   -0.1   0
W              0.01   0.0001     1
Amplitude         0.1 0      inf
eta:           2
"""]

global_input_string = """\
Bkgd:          3
two_theta_0       0.      -0.5  0.5
"""

bkgd_minimizer_input_string = """\
factr       1e10
iprint      -1
maxiter     150
m           10
pgtol       1e-5
epsilon     1e-13
"""

minimizer_input_string = """\
factr       1e2
iprint      -1
maxiter     150
m           10
pgtol       1e-5
epsilon     1e-13
"""

fine_minimizer_input_string = """\
factr       1e2
iprint      1
maxiter     150
m           15
pgtol       1e-5
epsilon     1e-13
"""

# tst_two_theta = []
# tst_y = []

# is_Sim_data = True #: Should be False unless simulated data 
#    #: (e.g. "Jade-AL2O3-Sim.xye") is used
# display_plots = True #: Only use to see sample plots
# # with open(r"17_05_23_0014_NIST SRM 1976b.xye") as file:
# # with open(r"16_01_07_0010_Aspirin_HighRez.xye") as file:
# # with open(r"16_03_09_0015_Silver Behenate.xye") as file:
# # os.path.dirname(__file__) + r
# with open(r"data//profiles//cement_15_03_11_0028.xye") as file:
#    for line in file.readlines()[1:]:
#       two_thetatmp, ytmp, ztmp = line.split()
#       # two_thetatmp, ytmp = line.split()
#       # if float(two_thetatmp) < 15:
#       tst_two_theta.append(float(two_thetatmp))
#       tst_y.append(float(ztmp)**2)
# tst_two_theta = np.array(tst_two_theta)
# # mask = np.ones(len(tst_two_theta),dtype=bool)
# mask = tst_two_theta > 20 
# # mask = np.logical_and(tst_two_theta >25,np.logical_or(tst_two_theta<33.75,
# #    tst_two_theta>34.3))
# # mask = np.logical_or(tst_two_theta<33.75,tst_two_theta>34.3)
# tst_two_theta = tst_two_theta[mask]
# tst_y = np.array(tst_y)[mask]

# def exercise_Rietveld_Refinery_Cement():
#    # RietveldPhase.fromstring(input_string) 
#    cifs = [
#       "1540705-Alite.cif", 
#       "9012789-Belite.cif", 
#       "1200009-Ferrite.cif", 
#       "1000039-AluminateCubic.cif", 
#       "9014308-AluminateOrtho.cif", 
#       "9007569-Arcanite.cif",
#       "1011094-FreeLime.cif", 
#       "1000053-Periclase.cif", 
#       ]
#    Rt = []

#    print "cifs: \n" 
#    for p in cifs:
#       print p
#    print "\nInput String: \n"
#    for i,p in enumerate(input_strings):
#       print "Phase " + str(i+1) + ": \n" + p
#    print "Global Input String: \n" + global_input_string
#    print "Minimizer Input String: \n" + minimizer_input_string
#    print "Fine Minimizer Input String: \n" + fine_minimizer_input_string


#    CU_wavelength = wavelengths.characteristic("CU").as_angstrom()
#    d_min = CU_wavelength/2/np.sin(np.pi/360*tst_two_theta[-1])
#    d_max = CU_wavelength/2/np.sin(np.pi/360*tst_two_theta[0])
#    tst_y_max = np.amax(tst_y)/len(cifs)

#    RietveldPhases.global_params_from_string(global_input_string,
#       tst_two_theta,tst_y)

#    for cif, input_string in zip(cifs,input_strings):
#    #    tt0 = time.time()
#       Rt.append(RietveldPhases(cif,input_string,d_min,d_max, \
#          I_max = tst_y_max, delta_theta=1.5,Intensity_Cutoff = 0.005))
   #    tt1 = time.time()
   #    print "TIME TO READ IN: " +str(tt1-tt0) + " seconds"

   # for i,Rp in enumerate(Rt):
   #    tmp = Rp.two_theta_peaks[np.abs(Rp.two_theta_peaks-34) <0.5]
   #    tmp2 = Rp.weighted_intensities[np.abs(Rp.two_theta_peaks-34) <0.5]
   #    print str(i) + ": " + str(tmp)
   #    print str(i) + ": " + str(tmp2)

   # numpeaks = 0
   # for i,Rp in enumerate(Rt):
   #    print Rp.two_theta_peaks.shape[0]
   #    numpeaks += Rp.two_theta_peaks.shape[0]
   # print numpeaks

   # First fit the background
   # RR = RietveldRefinery(Rt,bkgd_minimizer_input_string, \
   #    use_bkgd_mask=False,bkgd_delta_theta=0.05,
   #    store_intermediate_state=True, show_plots=True)
   # RR.display(RR.minimize_Bkgd)

   # #Now use the full dataset
   # RR = RietveldRefinery(Rt,minimizer_input_string,
   #    store_intermediate_state=True, show_plots=True)

   # # RR.display(RR.minimize_Bkgd)
   # # RR.display(RR.minimize_Bkgd_Offset)
   # # RR.display(RR.minimize_Amplitude)
   # # RR.display(RR.minimize_Amplitude)
   # RR.display(RR.minimize_Amplitude_Offset)
   # # RR.display(RR.minimize_Amplitude_Offset_unit_cell)
   # RR.display(RR.minimize_unit_cell)
   # # RR.display(RR.minimize_First_n_Phases)
   # # RR.display(RR.minimize_First_n_Phases,n=3)
   # # RR.display(RR.minimize_Amplitude_Offset_W)
   # RR.display(RR.minimize_Amplitude_Bkgd_Offset_W)
   # # RR.display(RR.minimize_Amplitude_Bkgd_Offset)
   # # RR.display(RR.minimize_only_Alite)
   # # RR.display(RR.minimize_All)
   # # RR.display(RR.minimize_All)
   # # RR.display(RR.minimize_All)
   # # RR.display(RR.minimize_All)
   # # RR.display(RR.minimize_All)

   # #For fine-tuning
   # RR2 = RietveldRefinery(RR.Phase_list,
   #    fine_minimizer_input_string,store_intermediate_state=True,show_plots=True)
   # RR2.display(RR2.minimize_All)
   # # RR2.display(RR2.minimize_All)
   # # RR2.display(RR2.minimize_All)
   # # RR2.display(RR2.minimize_All)


# def run():
#    exercise_Rietveld_Refinery_Cement()
#    print "OK"

# if (__name__ == "__main__"):
#    run()