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

Rt = []
x_list = []
Rt_list = []
mask_list = []
selections_list = []
final_params_list = []
RR = None

max_refinement_rounds = 5
num_displayed_params = 8

selections = np.zeros(
         (1,num_displayed_params,max_refinement_rounds),
         dtype=bool)

CU_wavelength = wavelengths.characteristic("CU").as_angstrom()

def set_refinement_masks():
   #assumes a refinement instance has been created
   masks = []
   for i in xrange(max_refinement_rounds):
      mask = copy.deepcopy(RR.mask)
      if selections[0,0,i]:
         mask = set_mask_by_label("bk",mask)
      if selections[0,1,i]:
         mask = set_mask_by_label("tw",mask)
      for j, phase_mask in enumerate(RR.phase_masks):
         if selections[j+1,2,i] or selections[0,2,i]:
            mask = set_mask_by_label("Sca",mask,phase_mask)
         if selections[j+1,3,i] or selections[0,3,i]:
            mask = set_mask_by_label("W",mask,phase_mask)
         if selections[j+1,4,i] or selections[0,4,i]:
            mask = set_mask_by_label("eta",mask,phase_mask)
         if selections[j+1,5,i] or selections[0,5,i]:
            mask = set_mask_by_label("V",mask,phase_mask)
         if selections[j+1,6,i] or selections[0,6,i]:
            mask = set_mask_by_label("U",mask,phase_mask)
         if selections[j+1,7,i] or selections[0,7,i]:
            mask = set_mask_by_label("uc_",mask,phase_mask)
      masks.append(mask)
   return masks

def set_refinement_mask():
   #assumes a refinement instance has been created
   global RR
   mask = copy.deepcopy(RR.mask)
   if np.any(selections[0,0,:]):
      mask = set_mask_by_label("bk",mask)
   if np.any(selections[0,1,:]):
      mask = set_mask_by_label("tw",mask)
   for j, phase_mask in enumerate(RR.phase_masks):
      if np.any(selections[j+1,2,:]) or np.any(selections[0,2,:]):
         mask = set_mask_by_label("Sca",mask,phase_mask)
      if np.any(selections[j+1,3,:]) or np.any(selections[0,3,:]):
         mask = set_mask_by_label("W",mask,phase_mask)
      if np.any(selections[j+1,4,:]) or np.any(selections[0,4,:]):
         mask = set_mask_by_label("eta",mask,phase_mask)
      if np.any(selections[j+1,5,:]) or np.any(selections[0,5,:]):
         mask = set_mask_by_label("V",mask,phase_mask)
      if np.any(selections[j+1,6,:]) or np.any(selections[0,6,:]):
         mask = set_mask_by_label("U",mask,phase_mask)
      if np.any(selections[j+1,7,:]) or np.any(selections[0,7,:]):
         mask = set_mask_by_label("uc_",mask,phase_mask)
   return mask

def set_mask_by_label(label,mask,phase_mask=None):
   if phase_mask is None:
      return np.logical_or(np.char.startswith(RR.x['labels'],label),mask)
   return np.logical_or(
      np.logical_and(np.char.startswith(RR.x['labels'],label),phase_mask),mask)

def set_checkbutton_states(selections,phase_index=None):
   if phase_index is None:
      phase_index = 0
   assert type(phase_index) == int

# class AutoScrollbar(tk.Scrollbar):
#     # a scrollbar that hides itself if it's not needed.  only
#     # works if you use the grid geometry manager.
#     def set(self, lo, hi):
#         if float(lo) <= 0.0 and float(hi) >= 1.0:
#             # grid_remove is currently missing from Tkinter!
#             self.tk.call("grid", "remove", self)
#         else:
#             self.grid()
#         tk.Scrollbar.set(self, lo, hi)
#     def pack(self, **kw):
#         raise TclError, "cannot use pack with this widget"
#     def place(self, **kw):
#         raise TclError, "cannot use place with this widget"

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

      self.num_phases=0


      self.history_frame = HistoryFrame(self.container)
      self.history_frame.grid(row=0,column=0)

      self.plot_frame = PlotFrame(self.container,self,padx=10,pady=10)
      self.plot_frame.grid(row=0,column=1,sticky='w')

      self.param_frame = ParamFrame(self.container,self,padx=10,nb_width=275)
      self.param_frame.grid(row=0,column=2)

      # temp. to allow for auto-loading of profile
      self.getProfile()
      self.getCifs()
      # self.param_frame.refine()

   def getCifs(self,filePaths=None):
      if filePaths is None:
         pass
         # self.filePaths = tkFileDialog.askopenfilenames(
         #    initialdir = "./data/cifs")
      else:
         self.filePaths = filePaths
      self.filePaths = [
      r".\data\cifs\Cement\1540705-Alite.cif",
      r".\data\cifs\Cement\1000039-AluminateCubic.cif",
      r".\data\cifs\Cement\9014308-AluminateOrtho.cif",
      # r".\data\cifs\Cement\9004096-anhydrite.cif",
      r".\data\cifs\Cement\9007569-Arcanite.cif",
      # r".\data\cifs\Cement\9005521-bassanite.cif",
      r".\data\cifs\Cement\9012789-Belite.cif",
      # r".\data\cifs\Cement\9009667-calcite.cif",
      r".\data\cifs\Cement\1200009-Ferrite.cif",
      r".\data\cifs\Cement\1011094-FreeLime.cif",
      r".\data\cifs\Cement\1000053-Periclase.cif",
      # r".\data\cifs\Cement\9000113-portlandite.cif",
      ]
      # self.filePaths = [
      # r".\data\cifs\1000032.cif",
      # r".\data\cifs\9015662-rutile.cif",
      # ]
      # self.filePaths = [r".\data\cifs\9015662-rutile.cif"]

      global Rt, selections
      Rt = []
      phase_names = []
      for filePath in self.filePaths:
         # cif_file_name = os.path.split(filePath)[1]
         phase = RietveldPhases(filePath, #I_max=I_max/len(self.filePaths),
            delta_theta=0.5,intensity_cutoff=0.01)
         Rt.append(phase)
         phase_names.append(phase.chemical_name)
         self.num_phases += 1

      self.param_frame.update_phase_combobox(phase_names)

      global RR,Rp
      RR = RietveldRefinery(Rt,Rp,
         store_intermediate_state=False, show_plots=False)

      self.history_frame.pie_figure.update_pie_chart(RR.composition_by_weight,
         phase_names)
      # Rp.updateplotprofile(RR.total_profile_state)

      selections = np.copy(np.broadcast_to(selections,
         (self.num_phases+1,)+selections.shape[1:]))

      self.param_frame.bkgd_control.state.set(1)
      self.param_frame.bkgd_control.checkbutton_clicked()
      # pdb.set_trace()



   def getProfile(self):
      # self.fileName = tkFileDialog.askopenfilename(
      #    initialdir = "./data/profiles")
      self.fileName = r".\\data\\profiles\\cement_15_03_11_0028.xye"
      # self.fileName = r".\\data\\profiles\\17_11_15_0004_CEMI425R_d6.xye"
      # self.fileName = r".\\data\\profiles\\Jade-Al2O3-Sim.xye"
      # self.fileName = r".\\data\\profiles\\d5_05005.xye"
      self.winfo_toplevel().title("Rietveld Refinement (" +
         os.path.split(self.fileName)[1]+")")

      # RietveldPhases.set_profile(self.fileName, min_two_theta=25)
      RietveldPhases.set_profile(self.fileName,
         min_two_theta=25,
         # number_of_columns=2,
         )

      global Rp
      Rp.setplotdata()

   def exit(self):
      self.destroy()

class HistoryFrame(tk.Frame):
   def __init__(self,parent,*args,**kwargs):
      tk.Frame.__init__(self,parent)
      self.parent = parent

      self.results_string = tk.StringVar()
      self.results_title = tk.Label(self,
         textvariable=self.results_string,font=('Verdana',12))
      self.results_string.set("History")
      self.results_title.grid(row=0,column=0,sticky='n')

      self.results_box_scrollbar = tk.Scrollbar(self)
      self.results_box_scrollbar.grid(row=1,column=1,sticky='ns',
         pady=10)

      self.results_box = tk.Listbox(self,
         activestyle='none',
         width=33,
         height=5,
         yscrollcommand=self.results_box_scrollbar.set
         )
      self.results_box.bind('<<ListboxSelect>>', self.onClick)
      self.results_box.bind('<Double-Button-1>', self.onDoubleClick)
      self.results_box.grid(row=1,column=0,sticky='nsew',
         pady=10)


      self.results_box_scrollbar.config(command=self.results_box.yview)

      # self.results_text_scrollbar = tk.Scrollbar(self)
      # self.results_text_scrollbar.grid(row=2,column=1,sticky='ns')

      # self.results_text = tk.Text(self,
      #    height=15,
      #    width=33,
      #    yscrollcommand=self.results_text_scrollbar.set,
      #    state=tk.DISABLED,
      #    wrap=tk.WORD,
      #    )
      # self.results_text.grid(row=2,column=0,sticky='ns')
      # self.results_text.insert(tk.END,"")

      # self.results_text_scrollbar.config(command=self.results_text.yview)

      self.pie_figure = PieFrame(self)
      self.pie_figure.grid(row=2,column=0,sticky='nesw')

      self.grid_rowconfigure(2,minsize=300)

   def onClick(self,event):
      global RR,Rt,Rp,x_list,Rt_list,selections_list
      # print self.results_box.curselection()[0]
      # pdb.set_trace()
      selected_index = self.results_box.curselection()[0]
      RR = RietveldRefinery(Rt_list[selected_index],
         Rp,
         # mask=mask_list[selected_index],
         # input_weights=RR.composition_by_weight,
         )
      RR.revert_to_x(x_list[selected_index])

      self.pie_figure.update_pie_chart(RR.composition_by_weight)

      # self.results_text.config(state=tk.NORMAL)
      # self.results_text.delete(0.0,tk.END)
      # self.results_text.insert(tk.END,RR.display_parameters())
      # self.results_text.config(state=tk.DISABLED)
      # self.param_string.set(RR.display_parameters())

   def onDoubleClick(self,event):
      global final_params_list
      text = final_params_list[self.results_box.curselection()[0]]
      self.popUp = PopUpParamBox(self,text)

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

   assert type(max_refinement_rounds) == int
   if max_refinement_rounds < 1:
      max_refinement_rounds = 1

   if max_refinement_rounds == 1:
      order_values.append(str(max_refinement_rounds))
      rounds_options.append((max_refinement_rounds-1,True))
   else:
      for i in xrange(max_refinement_rounds-2):
         order_values.append(str(i+1)+"-" + str(max_refinement_rounds))
         rounds_options.append((i,True))
         order_values.append(str(i+1)+"," + str(max_refinement_rounds))
         rounds_options.append((i,False))
      order_values.append(
            str(max_refinement_rounds-1) +"," + str(max_refinement_rounds))
      rounds_options.append((max_refinement_rounds-2,True))
      order_values.append(str(max_refinement_rounds))
      rounds_options.append((max_refinement_rounds-1,True))

   def __init__(self,parent,controller,
      default_round_start=1,
      select_all_rounds=True):
      self.parent = parent
      tk.Frame.__init__(self, parent)

      self.RoundsLabel = tk.Label(self, text = "Rounds: ")
      self.RoundsLabel.grid(row=0,column=0,sticky='e')

      # self.roundcheckbuttons = []
      self.round_states = []

      # for i in xrange(max_refinement_rounds):
      #    self.round_states.append(tk.IntVar())
      #    if i >= default_round_start-1:
      #       self.round_states[i].set(1)
      #    self.roundcheckbuttons.append(tk.Checkbutton(self,
      #       variable = self.round_states[i], text=str(i+1)))
      #    self.roundcheckbuttons[i].grid(row=0,column=i+1)

      self.round_selection = tk.StringVar(self)
      self.orderMenu = ttk.Combobox(self, textvariable=self.round_selection,
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
         (default_round_start-1,select_all_rounds)))
      self.orderMenu.bind("<<ComboboxSelected>>", self.order_selected)
      self.orderMenu.grid(row=0,column=1)

   def order_selected(self,event,rounds_mask=None):
      if rounds_mask is None:
         rounds_options = self.rounds_options[
            self.order_values.index(self.round_selection.get())]
         rounds_mask = np.zeros(max_refinement_rounds,dtype=bool)
         if rounds_options[1]:
            for i in xrange(max_refinement_rounds-rounds_options[0]):
               rounds_mask[rounds_options[0]+i] = True
         else:
            rounds_mask[rounds_options[0]] = True
            rounds_mask[-1] = True

      assert type(rounds_mask) == np.ndarray
      assert len(rounds_mask) == max_refinement_rounds

      combo_sel = int(self.parent.parent.master.master.phase_combobox.current())
      global selections
      if self.parent.index < 2: # i.e. is a global parameter control
         selections[:,self.parent.index,:] = rounds_mask
      else:
         selections[combo_sel,self.parent.index,:] = rounds_mask
      # print repr(selections) + '\n'
      # print repr(selections[0,:,:])
      # print np.all(np.equal(selections[0,:,:],selections[1,:,:]))

   def set_round_selection(self,rounds_mask):
      for i,item in enumerate(rounds_mask):
         if i < max_refinement_rounds-1:
            if item:
               if rounds_mask[i+1]:
                  self.orderMenu.current(RoundsBoxes.rounds_options.index(
                     (i,True)))
                  break
               else:
                  self.orderMenu.current(RoundsBoxes.rounds_options.index(
                     (i,False)))
                  break
         elif item:
            self.orderMenu.current(RoundsBoxes.rounds_options.index(
               (i,True)))

class RefinementParameterControl(tk.Frame):
   def __init__(self, parent, controller, index, text="",
      default_round_start=1,select_all_rounds=True,
      *args, **kwargs):
      tk.Frame.__init__(self, parent)
      self.parent = parent
      self.text = text
      self.index = index

      self.state = tk.IntVar()
      self.checkbutton = tk.Checkbutton(self, command=self.checkbutton_clicked,
         variable = self.state, text=text)
      self.checkbutton.grid(row=0,column=0,sticky='w')

      self.click_items = []

      self.rounds = RoundsBoxes(self,parent,
         default_round_start=default_round_start,
         select_all_rounds=select_all_rounds)
      self.click_items.append(self.rounds)

      self.is_phase_param = \
         any(list(map(lambda kv: isinstance(kv[1],ttk.Combobox),
            self.parent.children.iteritems())))

      # self.rounds.grid(row=0,column=2,sticky='w')

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
      checked = bool(self.state.get())

      global Rt
      if len(Rt) is not 0:
         if checked:
            self.rounds.order_selected(None)
         else:
            self.rounds.order_selected(None,
               rounds_mask=np.zeros(max_refinement_rounds,dtype=bool))

      # print list(map(lambda x: type(x.__class__),self.parent.children.values()))
      # print filter(lambda x: type(x) == ttk.Combobox,
      #    list(self.parent.children.values()))
      # is_phase_param = any(list(map(lambda kv: isinstance(kv[1],ttk.Combobox),
      #    self.parent.children.iteritems())))

      if checked:
         n = 1
         for item in self.click_items:
            item.grid(row=0,column=n)
            n += 1

         if self.text == "Bkgd.":
            try:
               RR = RietveldRefinery(Rt,Rp, \
                  bkgd_refine=True,
                  store_intermediate_state=False,
                  show_plots=False)
               RR.minimize_bkgd()
            except AttributeError:
               pass

         # print self.parent.phase_combobox.current()
         # self.rounds.grid()
         # self.initial.grid()
         # self.l_limit.grid()
         # self.u_limit.grid()
         # self.round_dropdownlist.grid()
      if not checked:
         for item in self.click_items:
            item.grid_remove()

         # print self.children
         # self.rounds.grid_remove()
         # self.initial.grid_remove()
         # self.l_limit.grid_remove()
         # self.u_limit.grid_remove()
         # self.round_dropdownlist.grid_remove()
      # if self.change_all:
      #    global Rt
      #    for phase in Rt:
      #       RietveldPhases.x['values'][index] = 0

class RadioRefinementParameterControl(RefinementParameterControl):
   def __init__(self, parent, controller, index, text="",
      default_round_start=1,*args, **kwargs):

      RefinementParameterControl.__init__(self,parent,controller, index,
         text=text, default_round_start=default_round_start,*args,**kwargs)

      self.radio_frame = tk.Frame(self)
      self.radiovar = tk.IntVar()
      tk.Radiobutton(self.radio_frame,text="Angular",variable=self.radiovar,
         value=1,command=self.radiobutton_switched).pack()
      tk.Radiobutton(self.radio_frame,text="Vertical",variable=self.radiovar,
         value=2,command=self.radiobutton_switched).pack()
      self.radiovar.set(1)
      self.click_items.insert(0,self.radio_frame)

      self.grid_columnconfigure(0,minsize=80)
      self.grid_rowconfigure(0,minsize=60)
      self.grid_columnconfigure(1,minsize=70)
      # self.grid_columnconfigure(2,minsize=110)

   # def checkbutton_clicked(self):
   #    if self.state.get() == 1:
   #       self.radio_frame.grid()
   #       self.rounds.grid()
   #    if self.state.get() == 0:
   #       self.radio_frame.grid_remove()
   #       self.rounds.grid_remove()

   def radiobutton_switched(self):
      if self.radiovar.get() == 1:
         RietveldPhases.set_vertical_offset(False)
      elif self.radiovar.get() == 2:
         RietveldPhases.set_vertical_offset(True)
      RietveldPhases.two_theta_0['values'] = 0

class PolynomRefinementParameterControl(RefinementParameterControl):
   def __init__(self, parent, controller, index,
      text="",
      default_round_start=1,
      select_all_rounds=True,
      default_order=2,
      *args, **kwargs):

      RefinementParameterControl.__init__(self, parent, controller, index,
         text=text,
         default_round_start=default_round_start,
         select_all_rounds=select_all_rounds, *args, **kwargs)

      self.order_dropdownlist = Dropdown_Int_List(self, parent,
         text="Order:", min_int=0, max_int=6, default_int=default_order)
      self.click_items.insert(0,self.order_dropdownlist)

      self.grid_columnconfigure(0,minsize=60)
      self.grid_columnconfigure(1,minsize=90)

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

      self.order_selected(None)

      # self.grid_columnconfigure(0,minsize=30)
      # self.grid_columnconfigure(1,minsize=50)

   def order_selected(self,tmp):
      # print self.order.get()
      order = int(self.order.get())+1
      if self.parent.index == 0:
         RietveldPhases.set_bkgd_order(order)
      elif self.parent.index == 4:
         phase_sel = self.parent.parent.master.master.phase_combobox.current()
         if phase_sel == 0:
            for i in xrange(len(Rt)):
               Rt[i].set_eta_order(order)
         else:
            Rt[phase_sel-1].set_eta_order(order)
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

      self.Scale = tk.Scale(self,
         from_=from_,
         to=to,
         # tickinterval=tickinterval,
         variable=self.value,
         orient=tk.HORIZONTAL,
         length=length,
         showvalue=0,
         resolution=resolution,
         )
      self.Scale.grid(row=0,column=2)
      self.Scale.set(initial)


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
         PolynomRefinementParameterControl(self.global_frame,self,0,
         text="Bkgd.",
         default_order=2,
         default_round_start=1,
         select_all_rounds=True,
         )
      self.param_controls.append(self.bkgd_control)
      self.bkgd_control.grid(row=0,column=0,sticky='w')

      self.offset_control = \
         RadioRefinementParameterControl(self.global_frame,self,1,
         text=u"2\u03b8 Corr.",
         default_round_start=1,
         select_all_rounds=True,
         )
      self.param_controls.append(self.offset_control)
      self.offset_control.grid(row=1,column=0,sticky='w')

      self.global_nb.add(self.global_frame,text="Global Parameters")
      self.global_nb.grid(row=0,column=0,columnspan=3,padx=10,pady=10)

      self.phase_nb = ttk.Notebook(self,height=220,width=self.nb_width)
      self.phase_frame = tk.Frame(self.phase_nb,padx=10,pady=10)

      self.phase_names = ['All']

      self.selection = tk.StringVar()
      self.selection.set(self.phase_names[0])
      self.phase_combobox = ttk.Combobox(self.phase_frame,
         textvariable=self.selection,
         values=self.phase_names,
         state='readonly',
         exportselection=0,
         width=min(len(max(self.phase_names,key=len)),30),
         )
      self.phase_combobox.bind("<<ComboboxSelected>>",self.on_phase_selected)
      self.phase_combobox.grid(row=0,column=0,sticky='w')

      self.Scale_control = RefinementParameterControl(self.phase_frame,
         self.parent,2,
         text="Scale",
         default_round_start=1,
         select_all_rounds=False,
         )
      self.param_controls.append(self.Scale_control)
      self.Scale_control.grid(row=1,column=0,sticky='w')

      self.W_control = RefinementParameterControl(self.phase_frame,self.parent,
         3,text="Caglioti W",
         default_round_start=2,
         select_all_rounds=False,
         )
      self.param_controls.append(self.W_control)
      self.W_control.grid(row=2,column=0,sticky='w')

      self.eta_control = PolynomRefinementParameterControl(
         self.phase_frame,self.parent, 4,
         text=u"\u03b7",
         default_order=2,
         default_round_start=2,
         select_all_round=False,
         )
      self.param_controls.append(self.eta_control)
      self.eta_control.grid(row=3,column=0,sticky='w')

      self.V_control = RefinementParameterControl(self.phase_frame,self.parent,
         5, text="Caglioti V",
         default_round_start=3,
         select_all_rounds=False,
         )
      self.param_controls.append(self.V_control)
      self.V_control.grid(row=4,column=0,sticky='w')

      self.U_control = RefinementParameterControl(self.phase_frame,self.parent,
         6, text="Caglioti U",
         default_round_start=4,
         select_all_rounds=False,
         )
      self.param_controls.append(self.U_control)
      self.U_control.grid(row=5,column=0,sticky='w')

      self.lattice_control = \
         RefinementParameterControl(self.phase_frame,self.parent, 7,
         text="Lattice Parameters",
         default_round_start=max_refinement_rounds,
         select_all_rounds=False,
         )
      self.param_controls.append(self.lattice_control)
      self.lattice_control.grid(row=6,column=0,sticky='w')

      self.phase_nb.add(self.phase_frame,text="Phase Parameters")
      self.phase_nb.grid(row=1,column=0,columnspan=3,padx=10,pady=10)

      self.iteration_Scale = LabelScale(self,parent,
         text="Max. number of iterations: ",
         from_=0,
         to=300,
         initial=100,
         length=80)
      self.iteration_Scale.grid(row=2,column=0,columnspan=3,padx=10,pady=10)

      self.RefineButton = ttk.Button(self, text='Refine',
         command=self.refine, takefocus=False, width=9)
      self.RefineButton.grid(row=3,column=0, padx=10, pady=10)

      self.RefineButton = ttk.Button(self, text='Refine All',
         command=self.refine_all, takefocus=False, width=9)
      self.RefineButton.grid(row=3,column=1, padx=10, pady=10)

      self.CancelButton = ttk.Button(self, text='Reset',
         command=self.reset, takefocus=False, width=9)
      self.CancelButton.grid(row=3,column=2, padx=10, pady=10)

   def refine(self,mask=None):
      global RR, Rt, Rp, x_list, Rt_list, selections_list, selections, \
         final_params_list

      maxiter = self.iteration_Scale.Scale.get()
      RR = RietveldRefinery(Rt,Rp, #input_weights=RR.composition_by_weight,
         maxiter=maxiter)

      if mask is None:
         mask = set_refinement_mask()

      RR.mask = mask
      assert len(RR.x) == len(RR.mask)
      # print RR.x[RR.mask]
      RR.minimize(callback_functions=(self.parent.master.update_idletasks,
         self.parent.master.update))
      time.sleep(0.1)
      # RR.display_parameters(RR.minimize)#mplitude_Bkgd_Offset)
      RR.display_stats(RR.minimize)#mplitude_Bkgd_Offset)

      final_params_list.append(RR.display_parameters())

      x_list.append(copy.deepcopy(RR.x))
      mask_list.append(copy.deepcopy(RR.mask))
      Rt_list.append(copy.deepcopy(Rt))
      selections_list.append(copy.deepcopy(selections))


      self.parent.master.history_frame.results_box.insert(self.numruns,
         "Run " + str(self.numruns+1) + ": " + str(RR.num_params) +
         " parameters, GoF = " + str(round(RR.GoF,3)))
      self.parent.master.history_frame.results_box.see(tk.END)

      self.parent.master.history_frame.pie_figure.update_pie_chart(
         RR.composition_by_weight)
      # self.parent.master.param_string.set(RR.display_parameters())
      # self.parent.master.history_frame.results_text.config(
      #    state=tk.NORMAL)
      # self.parent.master.history_frame.results_text.delete(0.0,tk.END)
      # self.parent.master.history_frame.results_text.insert(
      #    tk.END,RR.display_parameters())
      # self.parent.master.history_frame.results_text.config(
      #    state=tk.DISABLED)
      self.numruns += 1
      self.parent.master.update_idletasks()
      self.parent.master.update()

   def refine_all(self):
      masks = set_refinement_masks()
      for mask in masks:
         self.refine(mask)

   def reset(self):
      self.controller.num_phases = 0
      global Rt, RR, Rp, selections
      global x_list, mask_list, Rt_list, selections_list, final_params_list
      del Rt
      Rt = []
      if RR is not None:
         del RR
      RR = None

      del x_list
      x_list = []
      del mask_list
      mask_list = []
      del Rt_list
      Rt_list = []
      selections = copy.deepcopy(selections_list[-1])
      del selections_list
      selections_list = []
      del final_params_list
      final_params_list = []

      del self.parent.master.param_frame.phase_names
      self.parent.master.param_frame.update_phase_combobox([])

      # self.parent.master.history_frame.results_text.config(
      #    state=tk.NORMAL)
      # self.parent.master.history_frame.results_text.delete(0.0,tk.END)
      # self.parent.master.history_frame.results_text.config(
      #    state=tk.DISABLED)

      self.parent.master.history_frame.results_box.delete(0,tk.END)

      self.numruns = 0

      Rp.fig.suptitle("")
      self.controller.getCifs(self.controller.filePaths)

      # for p_c in self.parent.master.param_frame.param_controls:
      #    # p_c.state.set(0)
      #    p_c.checkbutton_clicked()

      self.on_phase_selected(None)


   def on_phase_selected(self,event):
      global selections
      phase_sel = int(self.phase_combobox.current())
      for i,control in enumerate(self.param_controls):
         rounds = selections[phase_sel,i,:]
         control.state.set(int(np.any(selections[phase_sel,i,:])))
         control.checkbutton_clicked()
         control.rounds.set_round_selection(rounds)
      self.phase_combobox.selection_clear()

   def update_phase_combobox(self,phase_names):
      self.phase_names = ['All'] + phase_names
      self.phase_combobox.grid_remove()
      self.phase_combobox = ttk.Combobox(
         self.phase_frame,
         textvariable=self.selection,
         values=self.phase_names,
         state='readonly',
         exportselection=0,
         width=min(len(max(self.phase_names,key=len))+1,30),
         )
      self.phase_combobox.bind(
         "<<ComboboxSelected>>",self.on_phase_selected)
      self.phase_combobox.grid(row=0,column=0,sticky='w')

class PlotFrame(tk.Frame):
   def __init__(self, parent,controller,*args,**kwargs):
      tk.Frame.__init__(self,parent,*args,**kwargs)

      global Rp
      Rp = RietveldPlot(width=5.5,height=4.5)
      self.canvas = FigureCanvasTkAgg(Rp.fig, master=self)
      self.canvas.get_tk_widget().pack(side=tk.TOP,anchor='n')
      # canvas.get_tk_widget().grid(row=0,column=0,sticky='n')

      toolbar = NavigationToolbar2TkAgg(self.canvas,self)
      toolbar.update()
      # toolbar.set_message("try this out")
      # canvas._tkcanvas.pack()
      # canvas._tkcanvas.grid()#row=0,column=0)#,sticky='n',pady=10)
      # canvas.show()

class PieFrame(tk.Frame):
   def __init__(self, parent,
      size=3, #the chart size, in inches
      pct_cutoff=5, #the percent cutoff when displaying compositions
      *args,**kwargs):
      tk.Frame.__init__(self,parent,*args,**kwargs)
      self.parent = parent
      self.pct_cutoff = pct_cutoff

      global RR
      self.labels = ['']
      self.compositions = [100]
      self.fig = Figure(figsize=(size,size),dpi=100)
      self.fig.suptitle('Composition (by Weight)')
      self.pie_chart = self.fig.add_subplot(111)

      self.canvas = FigureCanvasTkAgg(self.fig, master=self)
      self.canvas.get_tk_widget().pack(side=tk.TOP,anchor='n')

      self.phase_label = matplotlib.offsetbox.AnchoredText('Test',
         loc=8, frameon=True)
      self.phase_label.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
      self.phase_label.set_visible(False)

      self.update_pie_chart(self.compositions,labels=self.labels)
      self.canvas.mpl_connect('motion_notify_event',self.pick)
      self.canvas.mpl_connect('pick_event', self.hover)
      self.canvas.mpl_connect('axes_leave_event', self.leave)

      # self.pie_chart.annotate('test',xy=(0,0),xytext=(10,10))
      # self.phase_label
      # self.pie_chart.add_artist(self.phase_label)
      self.canvas.show()
      # canvas.get_tk_widget().grid(row=0,column=0,sticky='n')

      # toolbar = NavigationToolbar2TkAgg(self.canvas,self)
      # toolbar.update()
      # canvas._tkcanvas.pack()
      # canvas._tkcanvas.grid()#row=0,column=0)#,sticky='n',pady=10)
      # canvas.show()

   def pick(self,event):
      self.canvas.pick(event)

   def hover(self,event):
      if isinstance(event.artist,matplotlib.patches.Wedge):
         # print event.artist.get_label()
         self.phase_label.get_children()[0].set_text(event.artist.get_label())
         self.phase_label.set_visible(True)

         self.set_percent_visibilities(sel=self.wedges.index(event.artist))
         # else:
         #    self.percents[index].set_visible(False)
         self.canvas.draw()
         # for wedge in self.pie_chart_props:
         #    if event.artist.url == wedge[0].url:
         #       wedge[1].set_visible(True)

      #    print event.artist.properties()['label']
         # print event.artist.get_gid()
      # if event.inaxes == self.pie_chart:
      #    for prop in self.pie_chart_props:
      #       if prop[0].contains(event):
      #          print prop[1].get_text()

   def leave(self,event):
      self.phase_label.set_visible(False)
      self.canvas.draw()

   def update_pie_chart(self,compositions,labels=None):
      if labels is None:
         labels = self.labels
      else:
         self.labels = labels
      assert len(compositions) == len(labels)
      sorted_comps_and_labels = zip(*sorted(zip(compositions,labels),
                key=lambda x: x[0]))
      # self.compositions = sorted_comps_and_labels[0]
      # self.labels = sorted_comps_and_labels[1]

      self.pie_chart.clear()

      self.wedges, self.wedge_names, self.percents = \
         self.pie_chart.pie(sorted_comps_and_labels[0],
            labels=list(sorted_comps_and_labels[1]),
            shadow=True,
            startangle=90,
            wedgeprops={ 'picker': True }, #'clip_on': False,
            autopct='%.1f%%',
            )

      for name in self.wedge_names:
         name.set_visible(False)
      self.set_percent_visibilities()

      self.pie_chart.add_artist(self.phase_label)

      # self.fig.legend(patches,self.labels)
      # self.fig.tight_layout()
      self.canvas.draw_idle()

   def set_percent_visibilities(self,sel=None):
      if sel is not None:
         assert type(sel) == int
         sel_pct = self.percents[sel].get_text()[:-1]
      else: sel_pct = '0'
      for pct in self.percents:
         if float(pct.get_text()[:-1]) < self.pct_cutoff \
         and not pct.get_text()[:-1] == sel_pct:
            pct.set_visible(False)
         else: pct.set_visible(True)

   def picker(self,wedge,event):
      return wedge.contains(event), dict()

   # def autopct_cutoff(self,percents,cutoff=5):
   #    return ('%1.1f%%'% percents) if percents > cutoff else ''

class PopUpParamBox(tk.Toplevel):
   def __init__(self,parent,text):
      tk.Toplevel.__init__(self)
      self.resizable(width=False,height=False)
      self.parent = parent

      h = int(self.parent.master.winfo_height())
      w = int(self.parent.master.winfo_width())

      self.results_text_scrollbar = tk.Scrollbar(self)
      self.results_text_scrollbar.grid(row=0,column=1,sticky='ns')

      self.results_text = tk.Text(self,
         height=25,
         width=35,
         yscrollcommand=self.results_text_scrollbar.set,
         wrap=tk.WORD,
         )
      self.results_text_scrollbar.config(command=self.results_text.yview)
      self.results_text.grid(row=0,column=0,sticky='ns')

      self.results_text.insert(tk.END,text)
      self.results_text.config(state=tk.DISABLED)

      # self.message = tk.Message(self,text=text)
      # self.message.pack()

      self.ok = tk.Button(self, text="Got it", command=self.destroy)
      self.ok.grid(row=1,column=0,columnspan=2)

      self.geometry('+' + str(int(w/2)) + '+' + str(int(h/2)))

if __name__ == "__main__":
   root = RietveldGUI()
   # w, h = root.winfo_screenwidth(), root.winfo_screenheight()
   # root.geometry("%dx%d+0+0" % (w, h))
   # root.state('zoomed')
   root.call('wm', 'attributes', '.', '-topmost', True)
   root.after_idle(root.call, 'wm', 'attributes', '.', '-topmost', False)
   root.mainloop()
