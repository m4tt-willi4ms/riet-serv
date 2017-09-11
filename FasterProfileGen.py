from __future__ import division
#
# Command to run this example:
#   iotbx.python iotbx_cif.py
#
# See also:
#   http://cctbx.sourceforge.net/iotbx_cif
#
# from scitbx.array_family import flex
import os, random, math
import iotbx.cif, cctbx.miller
import scitbx
from cctbx import xray
from cctbx import crystal
from cctbx.array_family import flex
from cctbx.eltbx import wavelengths
import time
import sys, subprocess
import numpy as np
import matplotlib.pyplot as plt
from scitbx import lbfgsb
import jsonpickle
from libtbx import easy_pickle

def two_thetabar_squared(two_theta,two_thetapeak,U,V,W):
   # tan_theta_mid = math.tan(math.pi/360*two_theta_mid)
   tan_thetapeak = np.tan(math.pi/360.0*two_thetapeak)
   # omegaUVW_squared = abs(U*(tan_thetapeak-tan_theta_mid)**2 \
   #    +V*(tan_thetapeak-tan_theta_mid)+W)
   #tantwo_thetapeak = math.tan(math.pi/360*two_thetapeak)
   # if hasattr(two_theta, "__len__"):
   #   # omegaUVW_squared = U*flex.tan(math.pi/360.0*(two_theta-two_thetapeak))**2+V*flex.tan(math.pi/360*(two_theta-two_thetapeak))+W
   #   tantwo_theta = flex.tan(math.pi/360.0*(two_theta))
   # else:
   #   # omegaUVW_squared = U*math.tan(math.pi/360.0*(two_theta-two_thetapeak))**2+V*math.tan(math.pi/360*(two_theta-two_thetapeak))+W
   #   tantwo_theta = math.tan(math.pi/360.0*(two_theta))
   omegaUVW_squared = abs(U*tan_thetapeak**2+V*tan_thetapeak+W)
   return (two_theta-two_thetapeak)**2/omegaUVW_squared

def PseudoVoigtProfile(eta,two_theta0,U,V,W,Amp,two_theta,two_theta_calc_peak,I_res_calc):
   tan_thetapeak = np.tan(math.pi/360.0*two_theta_calc_peak)
   omegaUVW_squared = abs(U*tan_thetapeak**2+V*tan_thetapeak+W)
   two_thetabar_squared = (two_theta-two_theta0-two_theta_calc_peak)**2/omegaUVW_squared
   # tt_bar_sq = two_thetabar_squared(two_theta-two_theta0,two_theta_calc_peak \
      # ,U,V,W)
   return I_res_calc*Amp*(eta/(1 \
         +two_thetabar_squared) +(1-eta)*np.exp(-np.log(2)*two_thetabar_squared))

def Profile_Calc(x,two_theta,Rel_Peak_Intensity,delta_theta):
   # print Rel_Peak_Intensity
   two_theta_peaks = Rel_Peak_Intensity[:,0]
   two_theta_peaks.shape = (len(Rel_Peak_Intensity),1)
   Intensities = Rel_Peak_Intensity[:,1]
   Intensities.shape = (len(Rel_Peak_Intensity),1)
   mask = np.abs(two_theta - two_theta_peaks)< delta_theta
   result = np.zeros(len(two_theta))
   for i in xrange(0,len(Rel_Peak_Intensity),1):
      result[mask[i]] += PseudoVoigtProfile(
          x[0], # eta
         x[1], # two_theta_0
         x[2], # U
         x[3], # V
         x[4], # W
         x[5], # Amplitude
         two_theta[mask[i]],two_theta_peaks[i],Intensities[i])
   return result+Background_Polynomial(two_theta,np.array([x[6],x[7],x[8]]))

def showplot(filename,two_theta,x,y,Rel_Peak_Intensity,delta_theta):
   plt.figure(figsize=(12, 8))
   plt.subplot(3,1,1)
   # for i in xrange(0,len(two_theta)):
   #    if delta_theta[i]:
   #       color = 'blue'
   #    else: color = 'green'
   #    plt.scatter(two_theta[i],y[i], color=color,s=1)
   plt.scatter(two_theta,y,label='Data',s=1, color='red')
   plt.title(r"Profile: $Al_2 O_3$ - "+ filename)
   # plt.axis([20,60,0,1500])

   plt.plot(two_theta,Profile_Calc(x,two_theta,Rel_Peak_Intensity,delta_theta), label=r'$I_{\rm calc}$')
   plt.legend(bbox_to_anchor=(.8,.7))
   plt.ylabel(r"$I$")

   plt.subplot(3,1,2)
   # for i in xrange(0,len(two_theta)):
   #    if delta_theta[i]:
   #       color = 'blue'
   #    else: color = 'green'
   #    plt.scatter(two_theta[i],y[i], color=color,s=1)
   pltmask = np.logical_and(two_theta >25, two_theta < 26.25)
   plt.scatter(two_theta[pltmask],y[pltmask],label='Data',s=1, color='red')
   # plt.title(r"Profile: $Al_2 O_3$")
   # plt.axis([20,60,0,1500])

   plt.plot(two_theta[pltmask],Profile_Calc(x,two_theta,Rel_Peak_Intensity,delta_theta)[pltmask], label=r'$I_{\rm calc}$')
   plt.legend(bbox_to_anchor=(.8,.7))
   plt.ylabel(r"$I$")

   plt.subplot(3,1,3)
   zf = flex.double(len(two_theta))
   y_calc = Profile_Calc(x,two_theta,Rel_Peak_Intensity,delta_theta)
   for i in xrange(len(two_theta)):
      zf[i] = 1/y[i]*(y[i]-y_calc[i])**2
   plt.scatter(two_theta,zf)
   plt.ylabel(r"$\frac{1}{I} \, (I-I_{\rm calc})^2$")
   plt.xlabel(r'$2\,\theta$')

   plt.show()#block=False)

def WSS(two_theta,x,y,Rel_Peak_Intensity,delta_theta):
   y_calc = Profile_Calc(x,two_theta,Rel_Peak_Intensity,delta_theta)
   return np.sum(1/y*(y-y_calc)**2)

def WSS_grad(two_theta,x,y,f,epsilon,Rel_Peak_Intensity,delta_theta,mask):
   grad = flex.double(len(x))
   # print str(x.as_numpy_array())
   for j in xrange(0,len(x),1):
      if mask[j]:
         x[j] += epsilon
         # print j, x[j]
         grad[j] = (WSS(two_theta,x,y,Rel_Peak_Intensity,delta_theta)-f)/epsilon
         x[j] -= epsilon
         # print x[j], grad[j]
      else: 
         grad[j] = 0.0
   return grad

def Background_Polynomial(two_theta,x_bkgd):
   powers = np.array(range(len(x_bkgd)))
   powers.shape = (len(x_bkgd),1)
   # print str(two_theta)
   # print str(np.dot(x_bkgd,np.power(two_theta,powers)))
   return np.dot(x_bkgd,np.power(two_theta,powers))

def Rel_Peak_Intensity(x,fn,lammbda="CUA1"):
   # print str(fn) + ': '
   with open(fn, 'r') as file:
      tmp_as_cif = file.read()
   tmp_structure = iotbx.cif.reader(
    input_string=tmp_as_cif).build_crystal_structures()[fn[0:-4]]
   #tmp_miller = iotbx.cif.reader(input_string=tmp_as_cif).build_miller_arrays()
   # Let's use scattering factors from the International Tables
   unit_cell = tmp_structure.unit_cell()
   # u_cart = tmp_structure.scatterers().extract_u_cart_plus_u_iso(unit_cell)
   # print u_cart[4][5]

   # tmp_structure.show_scatterers()
   # print tmp_structure.scatterers()[0].fp
   anomalous_flag = True
   d_min = 0.8

   f_miller_set = tmp_structure.build_miller_set(anomalous_flag, \
      d_min=d_min).sort()
   # f_miller_indices = f_miller_set.indices() 

   for scatterer in tmp_structure.scatterers():
      if (scatterer.label == "Al1"):
         scatterer.scattering_type = "Al"
      if (scatterer.label == "O1"):
         scatterer.scattering_type = "O"

   # wavelength = 1.54
   # wavelength = 1.540593 
   wavelength = wavelengths.characteristic(lammbda).as_angstrom()
   # wavelength2 = wavelengths.characteristic("CUA2").as_angstrom()
   s_wave = 'Wavelength: %f' % wavelength  + ' Angstroms'
   number_of_atoms = 0
   s_chem_formula = 'Chemical Formula: '
   # print tmp_structure.unit_cell_content()
   for item in tmp_structure.unit_cell_content().items():
      number_of_atoms += item[1]
      s_chem_formula += '%s%.0f ' % (item[0], item[1])
   s_num_atoms = 'Number of atoms in unit cell: %d' % number_of_atoms
   
   if ("--Verbose" in sys.argv[1:]):
      print s_wave
      print s_num_atoms
      print s_chem_formula   

   # print tmp_structure.scattering_type_registry

   # print tmp_structure.scattering_dictionary_as_string()
   # tmp_structure.show_scatterers()
   # tmp_structure.replace_scatterers(tmp_structure.scatterers(), \
      # site_symmetry_table='existing')
   # tmp_structure.show_scatterers()
   tmp_structure.scattering_type_registry(table="it1992"), # "it1992", "wk1995" "n_gaussian"\
      # types_without_a_scattering_contribution=["?"])
   # tmp_structure.set_inelastic_form_factors( \
   #    photon=wavelengths.characteristic(lammbda),table="sasaki")
   tmp_structure.set_b_iso(value=0.5)
   tmp_structure.use_u_iso = True
   # print tmp_structure.show_scatterer_flags_summary()
   # print tmp_structure.scatterers()[0].fp
   f_calc =  tmp_structure.structure_factors(d_min=d_min, \
      anomalous_flag=anomalous_flag).f_calc().sort()
   f_miller_indices = f_calc.indices()
   f_two_thetas_calc = unit_cell.two_theta(f_miller_indices,wavelength,deg=True)
   # f_calc.show_summary().show_array()
   # f_calc.show_array()

   f_calc2 = f_calc.data() #f_miller_set.structure_factors_from_scatterers(tmp_structure, algorithm='direct').f_calc().sort().data()

   f_calc_sq = f_calc.as_intensity_array().sort().data() #.apply_scaling(target_max=100).data()
   f_calc_mult = f_calc.multiplicities().sort().data()
   f_d_spacings = f_calc.d_spacings().sort().data()

   result_data = flex.double()
   # Ifactor = flex.double()
   Imax = 0.0
   assert len(f_calc_sq) == len(f_two_thetas_calc)
   if lammbda == "CUA2":
      K_alpha_2_factor = x[10]
   else: K_alpha_2_factor = 1
   cos_pref_angle_sq = x[9]
   
   for i in xrange(0,len(f_two_thetas_calc),1):
      two_theta = f_two_thetas_calc[i]
      factor = f_calc_mult[i]*K_alpha_2_factor*abs((1+cos_pref_angle_sq*math.cos(math.pi/180*two_theta)**2) \
         /(1+cos_pref_angle_sq)/math.sin(math.pi/360*two_theta)/math.sin(math.pi/180*two_theta))
      I = f_calc_sq[i]*factor
      # Ifactor.append(factor)
      if (I > Imax):
         Imax = I
      # result_data.append(I)
      result_data.append(I)
      # print f_miller_indices[i], I
   f_calc_sq_res = cctbx.miller.array(f_miller_set,data = result_data) #.sort() #.data() #.apply_scaling(target_max=100)

   # print f_calc_sq_res.is_real_array()
   # print list(f_calc_sq_res)
   f_calc_sq_res =  f_calc_sq_res.data() #.apply_scaling(target_max=100)
   # f_calc_sq_res = f_calc.get_rescaled_intensities(f_two_thetas_calc) #.apply_scaling(target_max=100)
   # print f_calc_sq_res

   if lammbda == "CUA1":
      f = open(fn[0:-4]+".txt",'w')
   else: 
      f = open(fn[0:-4]+".txt",'a')
   tmp_structure.show_summary(f=f)
   f.write(s_wave + '\n')
   f.write(s_chem_formula + '\n')
   f.write(s_num_atoms + '\n')
   f.write('{:13} {:8} {:11} {:12} {:9} {:9} {:7} {:4}\n'.format('hkl', \
      '2Theta','d-spacing', 'F', '|F|^2', 'I', 'I (%)', 'mult'))
   # print zip(f_miller_indices, f_two_thetas_calc, f_d_spacings, f_calc2, \
         # f_calc_sq, f_calc_sq_res, f_calc_mult)
   for miller_index,two_theta,d_spacing,fc,fcsq,fcsqres,m in \
      zip(f_miller_indices, f_two_thetas_calc, f_d_spacings, f_calc2, \
         f_calc_sq, f_calc_sq_res, f_calc_mult):
      if (fcsqres/Imax > 0.001):
         s = '{:13} {:7.3f} {:8.4f} {:14.2f} {:8.2f} {:10.2f} {:8.2f} {:3}' \
         .format(miller_index,two_theta,d_spacing,fc,fcsq,fcsqres,fcsqres/Imax*100,m)
         f.write(s +'\n')
         if ("--Verbose" in sys.argv[1:]):
            print s
   f.close()
   return zip(f_two_thetas_calc,result_data)

def driver1(x,labels,mask,use_fortran_library=False,use_bkgd_mask=False):
   two_theta = []
   y = []
   with open(r"17_05_23_0014_NIST SRM 1976b.xye") as file:
   # with open(r"Jade-Al2O3-Sim.xye") as file:
      for line in file.readlines()[1:]:
         two_thetatmp, ytmp, ztmp = line.split()
         # two_thetatmp, ytmp = line.split()
         # if float(two_thetatmp) < 32.0:
         two_theta.append(float(two_thetatmp))
         y.append(float(ytmp))
 
   two_theta = np.array(two_theta)
   y = np.array(y)

   t2 = time.time()
   fn = "1000032.cif"
   # fn = "9007634.cif"
   Relative_Peak_Intensity = Rel_Peak_Intensity(x,fn) \
                              + Rel_Peak_Intensity(x,fn,"CUA2")
   Relative_Peak_Intensity = np.array(Relative_Peak_Intensity)

   t3 = time.time()
   print "CIF Card Read-time: " + str(round(t3-t2,4)) + " seconds\n"

   if use_bkgd_mask == True:
      bkgd_mask = two_theta < 0
      for peak in Relative_Peak_Intensity[:,0]:
         bkgd_mask = np.logical_or(bkgd_mask,abs(two_theta-peak)<0.5)
      bkgd_mask = np.invert(bkgd_mask)
      two_theta = two_theta[bkgd_mask]
      y = y[bkgd_mask]
   

   # print zip(*Relative_Peak_Intensity)
   # print Relative_Peak_Intensity
   # max_overlapping_profiles = 5
   # delta_theta = flex.int(len(two_theta*max_overlapping_profiles))
   # delta_theta.reshape(flex.grid(len(two_theta),max_overlapping_profiles))

   delta_theta = 0.5
   # set_data_flags(two_theta,zip(*Relative_Peak_Intensity)[0],delta_theta,delta_theta)

   # delta_theta.to_list()

   # for i in xrange(0,len(two_theta),1):
   #    print 'Two-theta = {}: '.format(two_theta[i]),
   #    for j in xrange(0,delta_theta.all()[1],1):
   #       print delta_theta[i,j],
   #    print ''
   # print list(delta_theta)
   # print delta_theta[]
   # print delta_theta.as_numpy_array().tostring()

   n = len(x)
   nbd = flex.int(n)
   # x = flex.double(n)
   l = flex.double(n)
   u = flex.double(n)
   g = flex.double(n)
   if ("--Verbose" in sys.argv[1:]):
      iprint = 1#1000
   else:
      iprint = -1
   for i in xrange(0,n,1):
      nbd[i] = 0
   # for i in xrange(0,n,1):
   #   nbd[i] = 2
   #   l[i] = 1.0e0
   #   u[i] = 1.0e2
   # for i in xrange(1,n,2):
   #   nbd[i] = 2
   #   l[i] = -1.0e2
   #   u[i] = 1.0e2
   nbd[0] = 2
   l[0] = 0
   u[0] = 1
   # nbd[1] = 2
   # l[1] = -5
   # u[1] = 5
   # nbd[4] = 2
   # l[4] = -0.0005
   # u[4] = 0.0005
   nbd[2] = 2
   l[2] = -x[4]
   u[2] = x[4]
   nbd[3] = 2
   l[3] = -x[4]
   u[3] = x[4]
   # nbd[3] = 2
   # l[3] = 0
   # u[3] = 2
   nbd[9] = 2
   l[9] = 0
   u[9] = 1
   nbd[10] = 2
   l[10] = 0.4
   u[10] = 0.6
   print "Limits:"
   print "l"+str(list(l))
   print 'u'+str(list(u))
   # for i in xrange(n):
   #   x[i] = 3.0e0
   minimizer = lbfgsb.minimizer(
   n=n,
   m=5,
   l=l,
   u=u,
   nbd=nbd,
   factr=1.0e+6,
   pgtol=1.0e-7,
   iprint=iprint)

   f = WSS(two_theta,x.as_numpy_array(),y,Relative_Peak_Intensity,delta_theta)
   epsilon = 1e-9
   g = WSS_grad(two_theta,x,y,f,epsilon,Relative_Peak_Intensity,delta_theta,mask)

   
   #"x_bkgd_1", "x_bkgd_2"]
   print "\nBEFORE:"
   for i in xrange(0,len(labels),1):
      print labels[i] + ": " + str(x[i])

   showplot(fn,two_theta,x,y,Relative_Peak_Intensity,delta_theta)

   t0 = time.time()

   while True:
      if (minimizer.process(x, f, g, use_fortran_library)):
         Relative_Peak_Intensity = Rel_Peak_Intensity(x,fn) \
                              + Rel_Peak_Intensity(x,fn,"CUA2")
         Relative_Peak_Intensity = np.array(Relative_Peak_Intensity)
         f = WSS(two_theta,x,y,Relative_Peak_Intensity,delta_theta)
         g = WSS_grad(two_theta,x,y,f,epsilon,Relative_Peak_Intensity,delta_theta,mask)
      elif (minimizer.is_terminated()):
         break

   t1 = time.time()
   totaltime = t1-t0

   print "\nAFTER:"
   for i in xrange(0,len(x),1):
      print labels[i] + ": " + str(x[i])

   sumobs = 0
   for i in xrange(0,len(two_theta),1):
      sumobs = sumobs + y[i]

   print "\nR_wp: " + str(np.sqrt(f/sumobs))
   print "R_e: " + str(np.sqrt((len(two_theta)-len(x))/sumobs))
   print "Goodness-of-Fit = " + str(np.sqrt(f/(len(two_theta)-len(x))))

   print "\nTime elapsed: " + str(round(totaltime,4)) + " seconds\n"

   showplot(fn,two_theta,x,y,Relative_Peak_Intensity,delta_theta)

def Rietveld_Refine(x_initial,nbd_total,l_total,u_total,two_theta,y,Relative_Peak_Intensity,refine_flags):
   
   is_background_refinement = False
   if (refine_flags == [False,False,False,False,False,False,True,True,True]):
      is_background_refinement = True
   
   n = 0
   for flag in refine_flags:
      if (flag):
         n += 1

   x = flex.double(n)
   nbd = flex.int(n)
   l = flex.double(n)
   u = flex.double(n)
   g = flex.double(n)

   index = 0
   for i, flag in enumerate(refine_flags):
      if (flag):
         x[index] = x_initial[i]
         nbd[index] = nbd_total[i]
         l[index] = l_total[i]
         u[index] = u_total[i]
         index += 1

   if ("--Verbose" in sys.argv[1:]):
      iprint = 1000
   else:
      iprint = -1

   delta_theta = 0.5
   two_theta_peaks = zip(*Relative_Peak_Intensity)[0]
   print two_theta_peaks

   delta_theta = set_data_flags(two_theta_total,two_theta_peaks,delta_theta)

   print "l"+str(list(l))
   print 'u'+str(list(u))

   minimizer = lbfgsb.minimizer(
   n=n,
   m=5,
   l=l,
   u=u,
   nbd=nbd,
   factr=1.0e+8,
   pgtol=1.0e-5,
   iprint=iprint)

   f = WSS(two_theta,x,y,Relative_Peak_Intensity,delta_theta)

   epsilon = 1e-8
   g = WSS_grad(two_theta,x,y,f,epsilon,Relative_Peak_Intensity,delta_theta)

   labels = ["eta", "two_thetapeak", "U", "V", "W", "Amplitude", "x_bkgd_0", \
   "x_bkgd_1", "x_bkgd_2"]
   print "BEFORE:"
   for i in xrange(0,len(labels),1):
      print labels[i] + ": " + str(x[i])

   exit()

   showplot(two_theta,x,y,Relative_Peak_Intensity,delta_theta)

   t0 = time.time()


   while True:
      if (minimizer.process(x, f, g, use_fortran_library)):
         f = WSS(two_theta,x,y,Relative_Peak_Intensity,delta_theta)
         g = WSS_grad(two_theta,x,y,f,epsilon,Relative_Peak_Intensity,delta_theta)
      elif (minimizer.is_terminated()):
         break

   t1 = time.time()
   totaltime = t1-t0

   print "\nAFTER:"
   for i in xrange(0,len(x),1):
      print labels[i] + ": " + str(x[i])

   sumobs = 0
   for i in xrange(0,len(two_theta),1):
      sumobs = sumobs + y[i]

   print "\nR_wp: " + str(np.sqrt(f/sumobs))
   print "R_e: " + str(np.sqrt((len(two_theta)-len(x))/sumobs))
   print "Goodness-of-Fit = " + str(np.sqrt(f/(len(two_theta)-len(x))))

   print "\nTime elapsed: " + str(round(totaltime,4)) + " seconds\n"

   showplot(two_theta,x,y,Relative_Peak_Intensity,delta_theta)

def run():
   # driver1(use_fortran_library=("--fortran" in sys.argv[1:]))
   two_theta = []
   y = []

   with open(r"17_05_23_0014_NIST SRM 1976b.xye") as file:
      for line in file.readlines()[1:]:
         two_thetatmp, ytmp, ztmp = line.split()#, ztmp
         if float(two_thetatmp) < 80.0:
            two_theta.append(float(two_thetatmp))
            y.append(float(ytmp))

   fn = "1000032.cif"
   Relative_Peak_Intensity = Rel_Peak_Intensity(fn) \
                              + Rel_Peak_Intensity(fn,"CUA2")

   n=9
   x = flex.double(n) 

   x[0] = 0.0 #eta
   x[1] = 0.1 #2-theta_0
   x[2] = 0.0 #U
   x[3] = 0.0 #
   x[4] = 0.001
   x[5] = 0.001
   x[6] = 0 #Bkgd_Polynom_0
   x[7] = 0 #Bkgd_Polynom_1
   x[8] = 0 #Bkgd_Polynom_2

   nbd = flex.int(n)
   # x = flex.double(n)
   l = flex.double(n)
   u = flex.double(n)
   g = flex.double(n)

   for i in xrange(0,n,1):
      nbd[i] = 0
   # for i in xrange(0,n,1):
   #   nbd[i] = 2
   #   l[i] = 1.0e0
   #   u[i] = 1.0e2
   # for i in xrange(1,n,2):
   #   nbd[i] = 2
   #   l[i] = -1.0e2
   #   u[i] = 1.0e2
   nbd[0] = 2
   l[0] = 0
   u[0] = 1
   nbd[4] = 2
   l[4] = -0.005
   u[4] = 0.005

   n_runs = 2

   refine_flags = []
   refine_flags.append([False,False,False,False,False,False,True,True,True])
   refine_flags.append([True,True,True,True,True,True,False,False,False])

   for r_f in refine_flags:
      x = Rietveld_Refine(x,nbd,l,u,two_theta,y,Relative_Peak_Intensity, \
         r_f)

class Logger(object):
   def __init__(self):
     self.terminal = sys.stdout
     self.log = open("Al2O3_output.log", "w")

   def write(self, message):
     self.terminal.write(message)
     self.log.write(message)

if (__name__ == "__main__"):
   sys.stdout = Logger()

     # n: Number of refinement parameters
   n = 11
   x = flex.double(n)
   labels = []#, "x_bkgd_0", \
   x[0] = 0.5 # eta
   labels.append("eta")
   # mask.append(True)
   x[1] = 0.0 # two_theta_0
   labels.append("two_thetapeak")
   # mask.append(True)
   x[2] = 0.0 # U
   labels.append("U")
   # mask.append(True)
   x[3] = 0.0 # V
   labels.append("V")
   # mask.append(True)
   x[4] = 0.0006 # W
   labels.append("W")
   # mask.append(True)
   x[5] = 0.001 # Amplitude
   labels.append("Amplitude")
   # mask.append(True)
   x[6] = 0 #Bkgd0
   labels.append("Bkgd0")
   # mask.append(True)
   x[7] = 0 #Bkgd1
   labels.append("Bkgd1")
   # mask.append(True)
   x[8] = 0 #Bkgd2
   labels.append("Bkgd2")
   # mask.append(True)
   x[9] = 1 #Cos_preferred_angle_sq
   labels.append("Cos_preferred_angle_squared")

   x[10] = 0.48 #K_alpha2_factor
   labels.append("K_alpha2_factor")

   driver1(x,labels,[
      True,   #eta
      True,    #2-theta_0
      False,   #U
      False,   #V
      False,   #W
      True,    #Amplitude
      False,   #Background_Polynomial_0
      False,   #Background_Polynomial_1
      False,   #Background_Polynomial_2
      False,   #Cos_preferred_angle_squared
      False    #K_alpha2_factor
      ],use_fortran_library=("--fortran" in sys.argv[1:])
      )

   # driver1(x,labels,[
   #    True,   #eta
   #    False,    #2-theta_0
   #    False,   #U
   #    False,   #V
   #    True,   #W
   #    True,    #Amplitude
   #    False,   #Background_Polynomial_0
   #    False,   #Background_Polynomial_1
   #    False   #Background_Polynomial_2
   #    ],use_fortran_library=("--fortran" in sys.argv[1:])
   #    )

   driver1(x,labels,[
      False,   #eta
      False,    #2-theta_0
      False,   #U
      False,   #V
      False,   #W
      False,    #Amplitude
      True,   #Background_Polynomial_0
      True,   #Background_Polynomial_1
      True,   #Background_Polynomial_2
      False,   #Cos_preferred_angle_squared
      False    #K_alpha2_factor
      ],use_fortran_library=("--fortran" in sys.argv[1:])
      ,use_bkgd_mask = True)


   driver1(x,labels,[
      True,   #eta
      False,    #2-theta_0
      True,   #U
      True,   #V
      True,   #W
      True,    #Amplitude
      False,   #Background_Polynomial_0
      False,   #Background_Polynomial_1
      False,   #Background_Polynomial_2
      True,   #Cos_preferred_angle_squared
      True    #K_alpha2_factor
      ],use_fortran_library=("--fortran" in sys.argv[1:])
      )
   # run()