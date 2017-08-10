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
   tan_thetapeak = math.tan(math.pi/360.0*two_thetapeak)
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


def PseudoVoigtProfile(x,two_theta,two_theta_calc_peak,I_res_calc):
   eta=x[0]
   two_theta0 = x[1]
   U = x[2]
   V = x[3]
   W = x[4]
   Amp = x[5]
   tt_bar_sq = two_thetabar_squared(two_theta,two_theta_calc_peak \
      +two_theta0,U,V,W)
   return I_res_calc*Amp*(eta/(1 \
         +tt_bar_sq) +(1-eta)*np.exp(-np.log(2)*tt_bar_sq))

def set_data_flags(two_theta,two_theta_peaks,y_data_flags,delta_theta=0.5):
   max_overlapping_profiles = y_data_flags.all()[1]
   for i, tt in enumerate(two_theta):
      num_overlapping_profiles = 0
      for j,tt_peak in enumerate(two_theta_peaks):
         if abs(tt-tt_peak) < delta_theta and num_overlapping_profiles < max_overlapping_profiles:
            y_data_flags[i,num_overlapping_profiles] = j
            num_overlapping_profiles += 1

def Profile_Calc(x,two_theta,Rel_Peak_Intensity,y_data_flags):
   result = flex.double(len(two_theta))
   for i in xrange(0,len(two_theta),1):
      tmp_val = 0
      for j in xrange(0,y_data_flags.all()[1],1):
         if not y_data_flags[i,j] == 0:
            tmp_val += PseudoVoigtProfile(x,two_theta[i], \
                              Rel_Peak_Intensity[y_data_flags[i,j]][0],Rel_Peak_Intensity[y_data_flags[i,j]][1])
      result[i] = tmp_val
   return result


# def ycalc_plot(x,two_theta):
#    # eta=x[0]
#    # two_thetapeak = x[1]
#    # Uval = x[2]
#    # Vval = x[3]
#    # Wval = x[4]
#    # Amp = x[5]
#    # yf = Amp*(eta/(1+two_thetabar_squared_plot(two_theta,two_thetapeak,Uval,Vval,Wval))+(1-eta)*np.exp(-np.log(2)*two_thetabar_squared_plot(two_theta,two_thetapeak,Uval,Vval,Wval)))
#    return yf
#    # for 2t in two_theta:


def showplot(filename,two_theta,x,y,Rel_Peak_Intensity,y_data_flags):
   plt.figure(figsize=(12, 8))
   plt.subplot(3,1,1)
   # for i in xrange(0,len(two_theta)):
   #    if y_data_flags[i]:
   #       color = 'blue'
   #    else: color = 'green'
   #    plt.scatter(two_theta[i],y[i], color=color,s=1)
   plt.scatter(two_theta,y,label='Data',s=1, color='red')
   plt.title(r"Profile: $Al_2 O_3$ - "+ filename)
   # plt.axis([20,60,0,1500])

   plt.plot(two_theta,Profile_Calc(x,two_theta,Rel_Peak_Intensity,y_data_flags), label=r'$I_{\rm calc}$')
   plt.legend(bbox_to_anchor=(.8,.7))
   plt.ylabel(r"$I$")

   plt.subplot(3,1,2)
   # for i in xrange(0,len(two_theta)):
   #    if y_data_flags[i]:
   #       color = 'blue'
   #    else: color = 'green'
   #    plt.scatter(two_theta[i],y[i], color=color,s=1)
   plt.scatter(two_theta,y,label='Data',s=1, color='red')
   # plt.title(r"Profile: $Al_2 O_3$")
   # plt.axis([20,60,0,1500])

   plt.plot(two_theta,Profile_Calc(x,two_theta,Rel_Peak_Intensity,y_data_flags), label=r'$I_{\rm calc}$')
   plt.legend(bbox_to_anchor=(.8,.7))
   plt.ylabel(r"$I$")

   plt.subplot(3,1,3)
   zf = flex.double(len(two_theta))
   y_calc = Profile_Calc(x,two_theta,Rel_Peak_Intensity,y_data_flags)
   for i in xrange(len(two_theta)):
      zf[i] = 1/y[i]*(y[i]-y_calc[i])**2
   plt.scatter(two_theta,zf)
   plt.ylabel(r"$\frac{1}{I} \, (I-I_{\rm calc})^2$")
   plt.xlabel(r'$2\,\theta$')

   plt.show()#block=False)

def WSS(two_theta,x,y,Rel_Peak_Intensity,y_data_flags):
   f = 0
   # set_data_flags(two_theta,zip(*Rel_Peak_Intensity)[0],y_data_flags)
   # y_calc = Profile_Calc(x,two_theta,Rel_Peak_Intensity,y_data_flags)
   for i in xrange(0,len(two_theta),1):
      for j in xrange(0,y_data_flags.all()[1],1):
         y_calc = 0
         if not y_data_flags[i,j] == 0:
            y_calc += PseudoVoigtProfile(x,two_theta[i], \
               Rel_Peak_Intensity[y_data_flags[i,j]][0],Rel_Peak_Intensity[y_data_flags[i,j]][1])
      f+= 1/y[i]*(y[i]-y_calc)**2
   return f

def WSS_grad(two_theta,x,y,f,epsilon,Rel_Peak_Intensity,y_data_flags):

   grad = flex.double(len(x))
   for j in xrange(0,len(x),1):
      x[j] = x[j]+epsilon
      grad[j] = (WSS(two_theta,x,y,Rel_Peak_Intensity,y_data_flags)-f)/epsilon
      x[j] = x[j]-epsilon
      print grad[j]
   return grad

def Background_Polynomial(two_theta,x_bkgd):
   result = flex.double(len(two_theta))
   for i in xrange(0,len(two_theta)):
      for j in xrange(0, len(x_bkgd)):
         result[i] += x_bkgd[j]*pow(two_theta[i],j)
   return result

def Rel_Peak_Intensity(fn,lammbda="CUA1"):
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
   d_min = 1

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
   print s_wave
   number_of_atoms = 0
   s_chem_formula = 'Chemical Formula: '
   # print tmp_structure.unit_cell_content()
   for x in tmp_structure.unit_cell_content().items():
      number_of_atoms += x[1]
      s_chem_formula += '%s%.0f ' % (x[0], x[1])
   s_num_atoms = 'Number of atoms in unit cell: %d' % number_of_atoms
   print s_num_atoms
   print s_chem_formula

   # print tmp_structure.scattering_type_registry

   # print tmp_structure.scattering_dictionary_as_string()
   # tmp_structure.show_scatterers()
   # tmp_structure.replace_scatterers(tmp_structure.scatterers(), \
      # site_symmetry_table='existing')
   # tmp_structure.show_scatterers()
   tmp_structure.scattering_type_registry(table="wk1995"), # "it1992", "wk1995" "n_gaussian"\
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
      K_alpha_2_factor = 0.48
   else: K_alpha_2_factor = 1
   
   for i in xrange(0,len(f_two_thetas_calc),1):
      two_theta = f_two_thetas_calc[i]
      factor = f_calc_mult[i]*K_alpha_2_factor*abs((1+math.cos(math.pi/180*two_theta)**2) \
         /math.sin(math.pi/360*two_theta)/math.sin(math.pi/180*two_theta))
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
      if (fcsqres/Imax > 0.01):
         s = '{:13} {:7.3f} {:8.4f} {:14.2f} {:8.2f} {:10.2f} {:8.2f} {:3}' \
         .format(miller_index,two_theta,d_spacing,fc,fcsq,fcsqres,fcsqres/Imax*100,m)
         f.write(s +'\n')
         print s
   f.close()
   return zip(f_two_thetas_calc,result_data)

def driver1(use_fortran_library=False):
   two_theta = []
   y = []
   with open(r"17_05_23_0014_NIST SRM 1976b.xye") as file:
   # with open(r"Jade-Al2O3-Sim.xye") as file:
      for line in file.readlines()[1:]:
         two_thetatmp, ytmp, ztmp = line.split()
         # two_thetatmp, ytmp = line.split()
         if float(two_thetatmp) < 80.0:
            two_theta.append(float(two_thetatmp))
            y.append(float(ytmp))
   n = 6
   x = flex.double(n)
   x[0] = 0.0
   x[1] = 0.05
   x[2] = 0.0
   x[3] = 0.0
   x[4] = 0.0001
   x[5] = 0.001
   # x[6] = 0
   # x[7] = 0
   # x[8] = 0
   #Test7

   fn = "1000032.cif"
   # fn = "9007634.cif"
   Relative_Peak_Intensity = Rel_Peak_Intensity(fn) \
                              + Rel_Peak_Intensity(fn,"CUA2")

   # print zip(*Relative_Peak_Intensity)
   # print Relative_Peak_Intensity
   max_overlapping_profiles = 5
   y_data_flags = flex.int(len(two_theta*max_overlapping_profiles))
   y_data_flags.reshape(flex.grid(len(two_theta),max_overlapping_profiles))

   delta_theta = 0.5
   set_data_flags(two_theta,zip(*Relative_Peak_Intensity)[0],y_data_flags,delta_theta)

   # for i in xrange(0,len(two_theta),1):
   #    print 'Two-theta = {}: '.format(two_theta[i]),
   #    for j in xrange(0,y_data_flags.all()[1],1):
   #       print y_data_flags[i,j],
   #    print ''
   # print list(y_data_flags)
   # print y_data_flags[]
   # print y_data_flags.as_numpy_array().tostring()

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
   nbd[4] = 2
   l[4] = -0.0005
   u[4] = 0.0005
   # nbd[2] = 2
   # l[2] = -0.005
   # u[2] = 0.005
   # nbd[3] = 2
   # l[3] = -0.005
   # u[3] = 0.005
   # nbd[3] = 2
   # l[3] = 0
   # u[3] = 2
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
   factr=1.0e+7,
   pgtol=1.0e-5,
   iprint=iprint)

   f = WSS(two_theta,x,y,Relative_Peak_Intensity,y_data_flags)
   epsilon = 1e-8
   g = WSS_grad(two_theta,x,y,f,epsilon,Relative_Peak_Intensity,y_data_flags)

   labels = ["eta", "two_thetapeak", "U", "V", "W", "Amplitude"]#, "x_bkgd_0", \
   #"x_bkgd_1", "x_bkgd_2"]
   print "BEFORE:"
   for i in xrange(0,len(labels),1):
      print labels[i] + ": " + str(x[i])

   showplot(fn,two_theta,x,y,Relative_Peak_Intensity,y_data_flags)

   t0 = time.time()

   while True:
      if (minimizer.process(x, f, g, use_fortran_library)):
         f = WSS(two_theta,x,y,Relative_Peak_Intensity,y_data_flags)
         g = WSS_grad(two_theta,x,y,f,epsilon,Relative_Peak_Intensity,y_data_flags)
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

   showplot(fn,two_theta,x,y,Relative_Peak_Intensity,y_data_flags)



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

   y_data_flags = set_data_flags(two_theta_total,two_theta_peaks,delta_theta)



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

   f = WSS(two_theta,x,y,Relative_Peak_Intensity,y_data_flags)
   epsilon = 1e-8
   g = WSS_grad(two_theta,x,y,f,epsilon,Relative_Peak_Intensity,y_data_flags)

   labels = ["eta", "two_thetapeak", "U", "V", "W", "Amplitude", "x_bkgd_0", \
   "x_bkgd_1", "x_bkgd_2"]
   print "BEFORE:"
   for i in xrange(0,len(labels),1):
      print labels[i] + ": " + str(x[i])

   exit()

   showplot(two_theta,x,y,Relative_Peak_Intensity,y_data_flags)

   t0 = time.time()


   while True:
      if (minimizer.process(x, f, g, use_fortran_library)):
         f = WSS(two_theta,x,y,Relative_Peak_Intensity,y_data_flags)
         g = WSS_grad(two_theta,x,y,f,epsilon,Relative_Peak_Intensity,y_data_flags)
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

   showplot(two_theta,x,y,Relative_Peak_Intensity,y_data_flags)

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
      # print r_f
      x = Rietveld_Refine(x,nbd,l,u,two_theta,y,Relative_Peak_Intensity, \
         r_f)
# print "OK"

class Logger(object):
   def __init__(self):
     self.terminal = sys.stdout
     self.log = open("Al2O3_output.log", "w")

   def write(self, message):
     self.terminal.write(message)
     self.log.write(message)

if (__name__ == "__main__"):
   sys.stdout = Logger()
   driver1(use_fortran_library=("--fortran" in sys.argv[1:]))
   # run()