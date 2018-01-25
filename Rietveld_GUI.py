from __future__ import division
# from scitbx.array_family import flex
# import scitbx.math
# from libtbx.test_utils import approx_equal, eps_eq, Exception_expected
# import sys
import numpy as np
import time
import pdb
import copy
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import \
   FigureCanvasTkAgg, NavigationToolbar2TkAgg
from  matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib import style, use
from matplotlib.widgets import SpanSelector
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
from src.RietveldRefinery import RietveldRefinery, RietveldPlot

from scitbx import lbfgsb
from cctbx.eltbx import wavelengths
from cctbx.sgtbx import lattice_symmetry

import Tkinter as tk
import tkFileDialog

LARGE_FONT = ("Verdana", 31)


# fig = Figure(dpi=100)
# fig, ax = plt.subplots(nrows=3,ncols=1,sharex=True,dpi=100)
# xdata, ydata = [], []
# ln, = plt.plot([], [], 'ro', animated=True)
# fig.suptitle(plottitle)
canvas = None
ROI_mask = None
ROI_center = 32
ROI_delta_theta = 3
x_mask = np.zeros(0,dtype=bool)

# def onselect(xmin, xmax):
#    x = RietveldPhases.two_theta
#    indmin, indmax = np.searchsorted(x,
#       (xmin, xmax))
#    indmax = min(len(x) - 1, indmax)

#    thisx = RietveldPhases.two_theta[indmin:indmax]
#    thisy = RietveldPhases.I[indmin:indmax]
#    # line2.set_data
#    # subplot2.axes.lines[0].set_data(thisx, thisy)
#    subplot2.set_xlim(thisx[0], thisx[-1])
#    subplot2.set_ylim(0, thisy.max())
#    fig.canvas.draw()

two_thetas = []
ys = []

Rt = []
x_list = []
Rt_list = []
mask_list = []
RR = None

# x_defaults = {}
# x_defaults['two_theta_0'] = np.array([('two_theta_0',0.0,-0.5,0.5)],
#    dtype=RietveldPhases.custom_dtype);
# x_defaults['U'] = np.array([('U',0.0,0,0.1)],
#    dtype=RietveldPhases.custom_dtype);
# x_defaults['V'] = np.array([('V',-0.00,-0.1,0)],
#    dtype=RietveldPhases.custom_dtype);

# x_default = np.empty(0,dtype=RietveldPhases.custom_dtype)

max_refinement_rounds = 5
num_phase_params = 6

CU_wavelength = wavelengths.characteristic("CU").as_angstrom()

class AutoScrollbar(tk.Scrollbar):
    # a scrollbar that hides itself if it's not needed.  only
    # works if you use the grid geometry manager.
    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            # grid_remove is currently missing from Tkinter!
            self.tk.call("grid", "remove", self)
        else:
            self.grid()
        tk.Scrollbar.set(self, lo, hi)
    def pack(self, **kw):
        raise TclError, "cannot use pack with this widget"
    def place(self, **kw):
        raise TclError, "cannot use place with this widget"

class RietveldGUI(tk.Tk):
   def __init__(self, *args, **kwargs):
      tk.Tk.__init__(self, *args, **kwargs)

      tk.Tk.iconbitmap(self, default=u"doc//ProtoLogo_R.ico")
      tk.Tk.wm_title(self, "Rietveld Refinement")
      tk.Tk.wm_geometry(self,'1260x600')

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


      self.history_frame = HistoryFrame(self.container)
      self.history_frame.grid(row=0,column=0)

      self.plot_frame = PlotFrame(self.container,self,padx=10,pady=10)
      self.plot_frame.grid(row=0,column=1,sticky='w')

      self.param_frame = ParamFrame(self.container,self,padx=10,
         nb_width=275
         )#,pady=10)
      self.param_frame.grid(row=0,column=2)
      # temp. to allow for auto-loading of profile
      # self.getProfile()
      # self.getCifs()
      # self.param_frame.refine()

   def getCifs(self,filePaths=None):
      if filePaths is None:
         # pass
         self.filePaths = tkFileDialog.askopenfilenames(
            initialdir = "./data/cifs")
      else:
         self.filePaths = filePaths
      # self.filePaths = [
      # r".\data\cifs\Cement\1540705-Alite.cif",
      # r".\data\cifs\Cement\1000039-AluminateCubic.cif",
      # r".\data\cifs\Cement\9014308-AluminateOrtho.cif",
      # r".\data\cifs\Cement\9004096-anhydrite.cif",
      # r".\data\cifs\Cement\9007569-Arcanite.cif",
      # r".\data\cifs\Cement\9005521-bassanite.cif",
      # r".\data\cifs\Cement\9012789-Belite.cif",
      # r".\data\cifs\Cement\9009667-calcite.cif",
      # r".\data\cifs\Cement\1200009-Ferrite.cif",
      # r".\data\cifs\Cement\1011094-FreeLime.cif",
      # r".\data\cifs\Cement\1000053-Periclase.cif",
      # r".\data\cifs\Cement\9000113-portlandite.cif",
      # ]
      # self.filePaths = [r".\data\cifs\9015662-rutile.cif"]


      global Rt, selections
      Rt = []
      selections = []
      for i,filePath in enumerate(self.filePaths):
         # cif_file_name = os.path.split(filePath)[1]
         Rt.append(RietveldPhases(filePath, #I_max=I_max/len(self.filePaths),
            delta_theta=0.8,intensity_cutoff=0.01))
         selections.append(np.zeros((num_phase_params,max_refinement_rounds+1),
            dtype=bool))
         self.param_frame.phase_names.append(Rt[-1].chemical_name)
         # if self.numPhases == 0:
         # self.param_frame.nb.add(
         #    RefinementParameterSet(self.param_frame.nb,self.param_frame,index=i),
         #    text=str(i+1)+" ")
         self.numPhases += 1

      self.param_frame.phase_combobox.grid_remove()
      self.param_frame.phase_combobox = ttk.Combobox(
         self.param_frame.phase_frame,
         textvariable=self.param_frame.selection,
         values=self.param_frame.phase_names,
         state='readonly',
         exportselection=0,
         width=min(len(max(self.param_frame.phase_names,key=len)),30),
         )
      self.param_frame.phase_combobox.grid(row=0,column=0,sticky='w')
      # self.param_frame.nb.add(RefinementParameterSet(self.param_frame.nb,
      #       self.param_frame),text="Phase Parameters")

      global RR,Rp
      RR = RietveldRefinery(Rt,Rp,
         store_intermediate_state=False, show_plots=False)

      self.param_frame.bkgd_control.state.set(1)
      self.param_frame.bkgd_control.checkbutton_clicked()
      # pdb.set_trace()

      Rp.updateplotprofile(RR.total_profile_state)


   def getProfile(self):
      self.fileName = tkFileDialog.askopenfilename(
         initialdir = "./data/profiles")
      # self.fileName = r".\\data\\profiles\\cement_15_03_11_0028.xye"
      # self.fileName = r".\\data\\profiles\\17_11_15_0004_CEMI425R_d6.xye"
      # self.fileName = r".\\data\\profiles\\d5_05005.xye"
      self.winfo_toplevel().title("Rietveld Refinement (" +
         os.path.split(self.fileName)[1]+")")

      RietveldPhases.set_profile(self.fileName, min_two_theta=25)
      global ROI_mask
      ROI_mask = np.abs(RietveldPhases.two_theta - ROI_center) < ROI_delta_theta


      global Rp
      Rp.setplotdata()

   def exit(self):
      self.destroy()

class HistoryFrame(tk.Frame):
   def __init__(self,parent,*args,**kwargs):
      tk.Frame.__init__(self,parent)

      self.results_string = tk.StringVar()
      self.results_title = tk.Label(self,
         textvariable=self.results_string,font=('Verdana',12))
      self.results_string.set("History")
      self.results_title.grid(row=0,column=0,sticky='n')

      self.results_box_scrollbar = AutoScrollbar(self)
      self.results_box_scrollbar.grid(row=1,column=1,sticky='ns',
         pady=10)

      self.results_box = tk.Listbox(self,
         activestyle='none',
         width=33,
         height=5,
         yscrollcommand=self.results_box_scrollbar.set
         )
      self.results_box.bind('<<ListboxSelect>>', self.onClick)
      self.results_box.grid(row=1,column=0,sticky='nsew',
         pady=10)


      self.results_box_scrollbar.config(command=self.results_box.yview)

      self.results_text_scrollbar = AutoScrollbar(self)
      self.results_text_scrollbar.grid(row=2,column=1,sticky='ns')

      self.results_text = tk.Text(self,
         height=15,
         width=33,
         yscrollcommand=self.results_text_scrollbar.set,
         state=tk.DISABLED,
         wrap=tk.WORD,
         )
      self.results_text.grid(row=2,column=0,sticky='ns')
      self.results_text.insert(tk.END,"")

      self.results_text_scrollbar.config(command=self.results_text.yview)

      self.grid_rowconfigure(2,minsize=300)

   def onClick(self,event):
      global RR,Rt,Rp,x_list,Rt_list
      # print self.results_box.curselection()[0]
      # pdb.set_trace()
      RR = RietveldRefinery(Rt_list[self.results_box.curselection()[0]],
         Rp,
         mask=mask_list[self.results_box.curselection()[0]],
         input_weights=RR.composition_by_weight)
      RR.revert_to_x(x_list[self.results_box.curselection()[0]])

      self.results_text.config(state=tk.NORMAL)
      self.results_text.delete(0.0,tk.END)
      self.results_text.insert(tk.END,RR.display_parameters())
      self.results_text.config(state=tk.DISABLED)
      # self.param_string.set(RR.display_parameters())

# class VarLabelEntry(tk.Frame):
#    def __init__(self,parent,text,x_label,index,*args, **kwargs):
#       tk.Frame.__init__(self,parent)
#       self.x_label = x_label
#       self.index = index
#       self.change_all = parent.change_all
#       self.is_l_limits = self.x_label == "l_limits"
#       self.is_u_limits = self.x_label == "u_limits"
#       self.is_values = self.x_label == "values"
#       self.passable = ["","-",".","-."]

#       # self.value = tk.StringVar()
#       # self.value.set(str(RietveldPhases.x[x_label][index]))
#       self.label = tk.Label(self,text=text)

#       # valid percent substitutions (from the Tk entry man page)
#       # note: you only have to register the ones you need; this
#       # example registers them all for illustrative purposes
#       #
#       # %d = Type of action (1=insert, 0=delete, -1 for others)
#       # %i = index of char string to be inserted/deleted, or -1
#       # %P = value of the entry if the edit is allowed
#       # %s = value of entry prior to editing
#       # %S = the text string being inserted or deleted, if any
#       # %v = the type of validation that is currently set
#       # %V = the type of validation that triggered the callback
#       #      (key, focusin, focusout, forced)
#       # %W = the tk name of the widget

#       vcmd = (self.register(self.onValidate),
#          '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
#       self.entry = tk.Entry(self,width=8,validate="key",validatecommand=vcmd)
#       self.entry.insert(0,str(RietveldPhases.x[x_label][index]))
#       self.label.grid(row=0,column=0)
#       self.entry.grid(row=0,column=1)

#    def onValidate(self, d, i, P, s, S, v, V, W):
#       global RR, Rp

#       try:
#          if all(P != x for x in self.passable):
#             val = float(P)
#             if self.is_l_limits:
#                if val <= RietveldPhases.x['u_limits'][self.index]:
#                   RietveldPhases.x[self.x_label][self.index] = val
#                else:
#                   raise ValueError
#             elif self.is_u_limits:
#                if val >= RietveldPhases.x['l_limits'][self.index]:
#                   RietveldPhases.x[self.x_label][self.index] = val
#                else:
#                   raise ValueError
#             elif self.is_values:
#                RietveldPhases.x[self.x_label][self.index] = val
#          else:
#             RietveldPhases.x[self.x_label][self.index] = 0

#       except(ValueError):
#          return False
#       if self.change_all:
#          for Rp in Rt:
#             RietveldPhases.x[self.x_label]
#       if RR is not None:
#          RR.update_state()
#       if Rp is not None:
#          Rp.updateplotprofile()
#       return True

class RoundsBoxes(tk.Frame):

   order_values = []
   rounds_options = []
   for i in xrange(max_refinement_rounds-2):
      order_values.append(str(i+1)+"-" + str(max_refinement_rounds))
      rounds_options.append((i+1,True))
      order_values.append(str(i+1)+"," + str(max_refinement_rounds))
      rounds_options.append((i+1,False))
   order_values.append(
         str(max_refinement_rounds-1) +"-" + str(max_refinement_rounds))
   rounds_options.append((max_refinement_rounds-1,True))
   order_values.append(str(max_refinement_rounds))
   rounds_options.append((max_refinement_rounds,True))

   def __init__(self,parent,controller,
      default_round_start=1,
      select_all_rounds=True):
      tk.Frame.__init__(self, parent)
      self.RoundsLabel = tk.Label(self, text = "Rounds: ")
      self.RoundsLabel.grid(row=0,column=0,sticky='e')

      self.roundcheckbuttons = []
      self.roundstates = []

      # for i in xrange(max_refinement_rounds):
      #    self.roundstates.append(tk.IntVar())
      #    if i >= default_round_start-1:
      #       self.roundstates[i].set(1)
      #    self.roundcheckbuttons.append(tk.Checkbutton(self,
      #       variable = self.roundstates[i], text=str(i+1)))
      #    self.roundcheckbuttons[i].grid(row=0,column=i+1)

      self.roundselection = tk.StringVar(self)
      self.orderMenu = ttk.Combobox(self, textvariable=self.roundselection,
         state='readonly',width =3)
      self.orderMenu['values'] = tuple(RoundsBoxes.order_values)

      assert type(default_round_start) == int
      if default_round_start > max_refinement_rounds:
         default_round_start = max_refinement_rounds
      if default_round_start <= 0:
         default_round_start = 1
      assert type(select_all_rounds) == bool
      if default_round_start >= max_refinement_rounds-1:
         select_all_rounds = True

      self.orderMenu.current(RoundsBoxes.rounds_options.index(
         (default_round_start,select_all_rounds)))
      self.orderMenu.bind("<<ComboboxSelected>>", self.order_selected)
      self.orderMenu.grid(row=0,column=1)

   def order_selected(self,event):
      pass

class RefinementParameterControl(tk.Frame):
   def __init__(self, parent, controller, text="",
      default_round_start=1,select_all_rounds=True,
      *args, **kwargs):
      tk.Frame.__init__(self, parent)
      self.parent = parent

      self.state = tk.IntVar()
      self.checkbutton = tk.Checkbutton(self, command=self.checkbutton_clicked,
         variable = self.state, text=text)
      self.checkbutton.grid(row=0,column=0,sticky='w')

      self.rounds = RoundsBoxes(self,parent,
         default_round_start=default_round_start,
         select_all_rounds=select_all_rounds)
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

class RefinementParameterRadioControl(tk.Frame):
   def __init__(self, parent, controller, text="",
      default_round_start=1,*args, **kwargs):
      tk.Frame.__init__(self, parent)
      self.parent = parent

      self.state = tk.IntVar()
      self.checkbutton = tk.Checkbutton(self, command=self.checkbutton_clicked,
         variable = self.state, text=text)
      self.checkbutton.grid(row=0,column=0,sticky='w')

      self.radioframe = tk.Frame(self)
      self.radiovar = tk.IntVar()
      tk.Radiobutton(self.radioframe,text="Angular",variable=self.radiovar,
         value=1,command=self.radiobutton_switched).pack()
      tk.Radiobutton(self.radioframe,text="Vertical",variable=self.radiovar,
         value=2,command=self.radiobutton_switched).pack()
      self.radioframe.grid(row=0,column=1)
      self.radiovar.set(1)

      self.rounds = RoundsBoxes(self,parent,default_round_start)
      self.rounds.grid(row=0,column=2,sticky='w')

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



      self.grid_columnconfigure(0,minsize=80)
      self.grid_rowconfigure(0,minsize=60)
      self.grid_columnconfigure(1,minsize=70)
      # self.grid_columnconfigure(2,minsize=110)

      self.checkbutton_clicked()

   def checkbutton_clicked(self):
      if self.state.get() == 1:
         self.radioframe.grid()
         self.rounds.grid()
         # self.initial.grid()
         # self.l_limit.grid()
         # self.u_limit.grid()
         # self.round_dropdownlist.grid()
      if self.state.get() == 0:
         self.radioframe.grid_remove()
         self.rounds.grid_remove()
         # self.initial.grid_remove()
         # self.l_limit.grid_remove()
         # self.u_limit.grid_remove()
         # self.round_dropdownlist.grid_remove()
      # if self.change_all:
      #    global Rt
      #    for phase in Rt:
      #       RietveldPhases.x['values'][index] = 0

   def radiobutton_switched(self):
      if self.radiovar.get() == 1:
         RietveldPhases.set_vertical_offset(False)
      elif self.radiovar.get() == 2:
         RietveldPhases.set_vertical_offset(True)
      RietveldPhases.two_theta_0['values'] = 0

# class LatticeParameterControl(tk.Frame):
#    def __init__(self,parent,controller,index):
#       tk.Frame.__init__(self,parent)
#       self.index = index
#       self.crystal_system = Rt[index].structure.space_group().crystal_system()

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
         text="Order:", min_int=0, max_int=6, default_int=default_order)
      self.order_dropdownlist.grid(row=0,column=1)

      self.rounds = RoundsBoxes(self,parent,default_round_start)
      self.rounds.grid(row=0,column=2)

      # self.round_dropdownlist = Dropdown_Int_List(self, parent,
      #    text="Round:", min_int=1, max_int=max_refinement_rounds,
      #    default_int=default_round)
      # self.round_dropdownlist.grid(row=0,column=2)

      self.grid_columnconfigure(0,minsize=60)
      self.grid_columnconfigure(1,minsize=90)
      # self.grid_columnconfigure(2,minsize=110)

      self.checkbutton_clicked()

   def checkbutton_clicked(self):
      global RR,Rp,canvas
      if self.state.get() == 1:
         self.order_dropdownlist.grid()
         self.rounds.grid()

         if self.text == "Bkgd.":
            try:
               RR = RietveldRefinery(Rt,Rp, \
                  bkgd_refine=True,
                  store_intermediate_state=False, show_plots=False)
               RR.minimize_bkgd()
            except AttributeError:
               pass
            # Rp.updateplotprofile(RietveldPhases.two_theta,
            #    RR.total_profile_state,wse=RR.weighted_squared_errors_state)
            # Rp.fig.canvas.draw()

      if self.state.get() == 0:
         self.order_dropdownlist.grid_remove()
         self.rounds.grid_remove()
         # if RR is not None:
         #    Rp.reset_plot_profile()
         # if self.text == "Background":
         #    if len(subplot1.axes.lines) is not 0:
         #       for x in (subplot1,subplot2,subplot3):
         #          for line in x.axes.lines:
         #             line.remove()
         #    canvas.show()

def startthread(func):
   # print "here"
   t = Thread(target=func)
   t.start()
   return t

class Dropdown_Int_List(tk.Frame):
   def __init__(self, parent, controller, text="", default_int=2, min_int=0,
         max_int=5,*args, **kwargs):
      tk.Frame.__init__(self, parent)
      self.parent = parent
      self.order_label = tk.Label(self,text=text+" ")
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

      # self.grid_columnconfigure(0,minsize=30)
      # self.grid_columnconfigure(1,minsize=50)

   def order_selected(self,tmp):
      # print self.order.get()
      if np.char.startswith(self.parent.parameter['labels'][0],"eta"):
         global Rt
         if self.parent.parent.index is not None:
            Rt[self.parent.parent.index].set_eta_order(int(self.order.get())+1)
         else:
            for i in xrange(len(Rt)):
               Rt[i].set_eta_order(int(self.order.get())+1)
      elif np.char.startswith(self.parent.parameter['labels'][0],"bkgd"):
         RietveldPhases.set_bkgd_order(int(self.order.get())+1)
      self.parent.checkbutton_clicked()
      self.orderMenu.selection_clear()

class LabelScale(tk.Frame):
   def __init__(self,parent,controller,
      text="",
      from_=0,
      to=200,
      tickinterval=50,
      length=100,
      resolution=10,
      initial=100,
      *args,**kwargs):
      tk.Frame.__init__(self,parent,*args,**kwargs)

      self.value = tk.IntVar()
      self.label = tk.Label(self,text=text)
      self.label.grid(row=0,column=0)

      self.vallabel=tk.Label(self,textvariable=self.value)
      self.vallabel.grid(row=0,column=1,sticky='w')

      self.scale = tk.Scale(self,
         from_=from_,
         to=to,
         # tickinterval=tickinterval,
         variable=self.value,
         orient=tk.HORIZONTAL,
         length=length,
         showvalue=0,
         resolution=resolution,
         )
      self.scale.grid(row=0,column=2)
      self.scale.set(initial)


      self.grid_columnconfigure(1,minsize=40)

class ParamFrame(tk.Frame):
   def __init__(self, parent,controller,
      nb_width=250,*args,**kwargs):
      # controller.minsize(width=400,height=400)
      tk.Frame.__init__(self,parent,*args,**kwargs)

      self.parent = parent
      self.controller = controller
      self.nb_width = nb_width
      self.global_nb = ttk.Notebook(self,height=100,width=self.nb_width)
      self.global_frame = tk.Frame(self.global_nb,padx=10,pady=10)
      self.numruns = 0
      #, width=100, height=100)
      # self.globalLabelFrame.p#grid(row=0,column=0)

      self.param_controls = []

      self.bkgd_control = \
         RefinementParameterPolynomControl(self.global_frame,self,
         text="Bkgd.",
         default_order=2,default_round_start=1)
      self.param_controls.append(self.bkgd_control)
      self.bkgd_control.grid(row=0,column=0,sticky='w')

      self.offset_control = \
         RefinementParameterRadioControl(self.global_frame,self,
         text=u"2\u03b8 Corr.")
      self.param_controls.append(self.offset_control)
      self.offset_control.grid(row=1,column=0,sticky='w')

      self.global_nb.add(self.global_frame,text="Global Parameters")
      self.global_nb.grid(row=0,column=0,columnspan=2,padx=10,pady=10)

      self.phase_nb = ttk.Notebook(self,height=220,width=self.nb_width)
      self.phase_frame = tk.Frame(self.phase_nb,padx=10,pady=10)

      self.phase_names = ['All']
      self.change_all = True

      self.selection = tk.StringVar()
      self.selection.set(self.phase_names[0])
      self.phase_combobox = ttk.Combobox(self.phase_frame,
         textvariable=self.selection,
         values=self.phase_names,
         state='readonly',
         exportselection=0,
         width=min(len(max(self.phase_names,key=len)),30),
         )
      self.phase_combobox.bind("<<ComboboxSelected>>",self.onPhaseSelected)
      self.phase_combobox.grid(row=0,column=0,sticky='w')

      self.scale_control = RefinementParameterControl(self.phase_frame,
         self.parent, text="Scale",change_all=self.change_all,
         default_round_start=1)
      self.param_controls.append(self.scale_control)
      self.scale_control.grid(row=1,column=0,sticky='w')

      self.W_control = RefinementParameterControl(self.phase_frame,self.parent,
         text="Caglioti W",default_round_start=2)
      self.param_controls.append(self.W_control)
      self.W_control.grid(row=2,column=0,sticky='w')

      self.eta_control = RefinementParameterPolynomControl(
         self.phase_frame,self.parent,
         text=u"\u03b7",
         default_order=1,
         default_round_start=2)
      self.param_controls.append(self.eta_control)
      self.eta_control.grid(row=3,column=0,sticky='w')

      self.V_control = RefinementParameterControl(self.phase_frame,self.parent,
         text="Caglioti V",
         default_round_start=3,
         select_all_rounds=False)
      self.param_controls.append(self.V_control)
      self.V_control.grid(row=4,column=0,sticky='w')

      self.U_control = RefinementParameterControl(self.phase_frame,self.parent,
         text="Caglioti U",
         default_round_start=4,
         select_all_rounds=False)
      self.param_controls.append(self.U_control)
      self.U_control.grid(row=5,column=0,sticky='w')

      self.lattice_control = \
         RefinementParameterControl(self.phase_frame,self.parent,
         text="Lattice Parameters",
         default_round_start=max_refinement_rounds,
         select_all_rounds=False)
      self.param_controls.append(self.lattice_control)
      self.lattice_control.grid(row=6,column=0,sticky='w')

      self.phase_nb.add(self.phase_frame,text="Phase Parameters")
      self.phase_nb.grid(row=1,column=0,columnspan=2,padx=10,pady=10)

      self.labelscale = LabelScale(self,parent,
         text="Max number of iterations: ",
         from_=0,
         to=300,
         initial=100,
         length=80)
      self.labelscale.grid(row=2,column=0,columnspan=2,padx=10,pady=10)

      self.RefineButton = ttk.Button(self, text='Refine',
         command=self.refine, takefocus=False)
      self.RefineButton.grid(row=3,column=0, padx=10, pady=10)

      self.CancelButton = ttk.Button(self, text='Reset',
         command=self.reset, takefocus=False)
      self.CancelButton.grid(row=3,column=1, padx=10, pady=10)

   def refine(self):
      global RR, Rt, Rp, x_list, Rt_list

      for rp in self.nb.children[self.nb.tabs()[0].split('.')[-1]] \
         .children.values():
            if isinstance(rp,RefinementParameterControl) \
               or isinstance(rp,RefinementParameterPolynomControl):
               if rp.parameter['labels'][0][0:3] == "uc_":
                  for i in xrange(len(Rt)):
                     if rp.state.get() == 1 and \
                        RR.composition_by_weight[i] >= RR.composition_cutoff:
                        Rt[i].recompute_peak_positions = True
                     else: Rt[i].recompute_peak_positions = False

      for i,tab in enumerate(self.nb.tabs()[1:]):
         for rp in self.nb.children[tab.split('.')[-1]].children.values():
            if isinstance(rp,RefinementParameterControl):
               if rp.parameter['labels'][0][0:3] == "uc_" and \
                  RR.composition_by_weight[i] >= RR.composition_cutoff:
                  if rp.state.get() == 1:
                     Rt[i].recompute_peak_positions = True

      maxiter = self.labelscale.scale.get()
      RR = RietveldRefinery(Rt,Rp,input_weights=RR.composition_by_weight,
         maxiter=maxiter)
      # RR.mask[0] = True

      for rp in self.nb.children[self.nb.tabs()[0].split('.')[-1]] \
         .children.values():
            if isinstance(rp,RefinementParameterControl) \
               or isinstance(rp,RefinementParameterPolynomControl):
               if rp.state.get() == 1:
                  for i in xrange(len(Rt)):
                     # print rp.parameter['labels'][0][0:3]
                     if RR.composition_by_weight[i] >= RR.composition_cutoff \
                        or np.any(
                           np.char.startswith(rp.parameter['labels'],"Amp")) \
                        or rp.parameter['labels'][0] == "W":
                        RR.mask = np.logical_or(RR.mask,
                           np.logical_and(RR.phase_masks[i],
                           np.char.startswith(RR.x['labels'],
                              rp.parameter['labels'][0][0:3])))

      for i,tab in enumerate(self.nb.tabs()[1:]):
         for rp in self.nb.children[tab.split('.')[-1]].children.values():
            if isinstance(rp,RefinementParameterControl) or \
               isinstance(rp,RefinementParameterPolynomControl):
               if rp.state.get() == 1:
                  # print rp.parameter['labels'][0][0:3]
                  if RR.composition_by_weight[i] >= RR.composition_cutoff \
                     or np.any(
                           np.char.startswith(rp.parameter['labels'],"Amp")) \
                     or rp.parameter['labels'][0] == "W":
                     RR.mask = np.logical_or(RR.mask,
                        np.logical_and(RR.phase_masks[i],
                        np.char.startswith(RR.x['labels'],
                           rp.parameter['labels'][0][0:3])))

      for rp in self.globalnb.children[self.globalnb.tabs()[0].split('.')[-1]] \
         .children.values():
         if rp.state.get() == 1:
            RR.mask = np.logical_or(RR.mask,
               np.logical_and(RR.global_mask,np.char.startswith(RR.x['labels'],
                  rp.parameter['labels'][0][0:3])))

      RR.minimize()

      x_list.append(copy.deepcopy(RR.x))
      mask_list.append(copy.deepcopy(RR.mask))
      Rt_list.append(copy.deepcopy(Rt))

      # RR.display_parameters(RR.minimize)#mplitude_Bkgd_Offset)
      RR.display_stats(RR.minimize)#mplitude_Bkgd_Offset)

      self.parent.master.results_box.insert(self.numruns,
         "Run " + str(self.numruns+1) + ": " + str(RR.num_params) +
         " parameters, GoF = " + str(round(RR.GoF,3)))
      self.parent.master.results_box.see(tk.END)
      # self.parent.master.param_string.set(RR.display_parameters())
      self.parent.master.results_text.config(state=tk.NORMAL)
      self.parent.master.results_text.delete(0.0,tk.END)
      self.parent.master.results_text.insert(tk.END,RR.display_parameters())
      self.parent.master.results_text.config(state=tk.DISABLED)
      self.numruns += 1

   def reset(self):
      self.controller.numPhases = 0
      global Rt, RR, Rp
      Rt = []
      if RR is not None:
         del RR
      RR = None
      # RietveldPhases.set_profile()
      for tab in self.nb.tabs():
         self.nb.hide(tab)
         self.nb.forget(tab)
      # Rp.reset_plot_profile()
      Rp.fig.suptitle("")
      self.controller.getCifs(self.controller.filePaths)

   def onPhaseSelected(self,event):
      if self.phase_combobox.current() == 0:
         self.change_all = True
      else:
         self.change_all = False
      print self.phase_combobox.current()

class PlotFrame(tk.Frame):
   def __init__(self, parent,controller,*args,**kwargs):
      tk.Frame.__init__(self,parent,*args,**kwargs)

      # global canvas, onselect, subplot1, span
      global Rp,canvas
      Rp = RietveldPlot()
      canvas = FigureCanvasTkAgg(Rp.fig, master=self)
      canvas.get_tk_widget().pack(side=tk.TOP,anchor='n')
      # canvas.get_tk_widget().grid(row=0,column=0,sticky='n')

      toolbar = NavigationToolbar2TkAgg(canvas,self)
      toolbar.update()
      # canvas._tkcanvas.pack()
      # canvas._tkcanvas.grid()#row=0,column=0)#,sticky='n',pady=10)
      # canvas.show()


if __name__ == "__main__":
   root = RietveldGUI()
   # w, h = root.winfo_screenwidth(), root.winfo_screenheight()
   # root.geometry("%dx%d+0+0" % (w, h))
   # root.state('zoomed')
   root.call('wm', 'attributes', '.', '-topmost', True)
   root.after_idle(root.call, 'wm', 'attributes', '.', '-topmost', False)
   root.mainloop()