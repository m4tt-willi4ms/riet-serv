# import matplotlib
# matplotlib.use("TkAgg")
# matplotlib.style.use("ggplot")
# from matplotlib.backends.backend_tkagg import \
#    FigureCanvasTkAgg, NavigationToolbar2TkAgg
# from  matplotlib.figure import Figure
# import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector
import numpy as np

from src.rietveld_phases import RietveldPhases

class RietveldPlot:
   """This class is used to group together objects associated with the
   matplotlib plots used in displaying the results of RietveldRefinery
   instances.
   """
   def __init__(self,
      width=4,height=3, #the figure width and height, in inches
      ):
      self.fig = plt.figure()#Figure(figsize=(width,height),dpi=100)
      self.subplot1 = self.fig.add_subplot(2,1,1)
      self.subplot2 = self.fig.add_subplot(2,1,2)

   def setplotdata(self,ROI_center=None,delta_theta=3):

      two_theta = RietveldPhases.two_theta
      I = RietveldPhases.I
      if ROI_center is None:
         ROI_center = two_theta[np.argmax(I)]
      ROI_mask = np.abs(two_theta-ROI_center) < delta_theta
      self.I_max = I.max()
      self.subplot1.clear()
      self.subplot1.scatter(two_theta,I,label='Data',s=1, color='blue')
      # self.subplot1.legend(bbox_to_anchor=(.8,.7))
      self.subplot1.set_ylabel(r"$I$",rotation=0)
      self.profile1, = self.subplot1.plot(two_theta,np.zeros(len(two_theta)),
         label=r'$I_{\rm calc}$',alpha=0.7,color='red')
      self.subplot1.axhline(y=-self.I_max/2,
         # xmin=two_theta[0],
         # xmax=two_theta[-1],
         color='black')
      self.diff1, = self.subplot1.plot(
         two_theta,-self.I_max/2*np.ones(len(two_theta)),
         label=r'$\Delta I$',color='green')
      self.subplot1.legend(bbox_to_anchor=(.8,.7))

      self.subplot2.clear()
      self.subplot2.scatter(two_theta,I,label='Data',s=1, color='blue')
      self.subplot2.set_ylabel(r"$I$",rotation=0)
      self.profile2, = self.subplot2.plot(two_theta,np.zeros(len(two_theta)),
         alpha=0.7,color='red')
      self.subplot2.axhline(y=-self.I_max/2,
         # xmin=two_theta[0],
         # xmax=two_theta[-1],
         color='black')
      self.diff2, = self.subplot2.plot(two_theta,
         -self.I_max/2*np.ones(len(two_theta)),
         label=r'$\Delta I$',color='green')

      self.subplot2.set_xlim(two_theta[ROI_mask][0],two_theta[ROI_mask][-1])
      self.subplot2.axes.set_xlabel(r'$2\,\theta$')

      plt.ion()
      self.fig.canvas.show()

   def updateplotprofile(self, profile, wse=None):
      self.profile1.set_ydata(profile)
      self.profile2.set_ydata(profile)
      if wse is not None:
         sigma = RietveldPhases.sigma
         self.diff1.set_ydata(-self.I_max/2+wse*sigma**2)
         self.diff2.set_ydata(-self.I_max/2+wse*sigma**2)

      self.fig.canvas.show()

      def onselect(xmin, xmax):
         x = RietveldPhases.two_theta
         indmin, indmax = np.searchsorted(x,
            (xmin, xmax))
         indmax = min(len(x) - 1, indmax)

         thisx = RietveldPhases.two_theta[indmin:indmax]
         # thisy = RietveldPhases.I[indmin:indmax]
         # line2.set_data
         # subplot2.axes.lines[0].set_data(thisx, thisy)
         self.subplot2.set_xlim(thisx[0], thisx[-1])
         self.subplot2.axes.set_ylim(
            top=1.07*max(np.max(profile[indmin:indmax]),
               np.max(RietveldPhases.I[indmin:indmax])))
         self.fig.canvas.draw_idle()

      self.span = SpanSelector(self.subplot1, onselect, 'horizontal',
                    rectprops=dict(alpha=0.5, facecolor='green'))

      self.subplot1.axes.relim()
      self.subplot1.axes.autoscale_view()
      self.subplot2.axes.relim()
      self.subplot2.axes.autoscale_view()
      # self.subplot3.axes.relim()
      # self.subplot3.axes.autoScale_view()

      # self.fig.canvas.draw_idle()
      self.fig.canvas.draw()

   def reset_plot_profile(self):
      if len(self.subplot1.axes.lines) is not 0:
         for x in (self.subplot1,self.subplot2):#,self.subplot3):
            for line in x.axes.lines:
               line.remove()
      self.fig.suptitle("")
      self.subplot1.axes.legend_.remove()
      self.fig.canvas.draw_idle()
