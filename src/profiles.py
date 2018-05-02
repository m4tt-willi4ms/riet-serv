import numpy as np

def pseudo_voigt(x_squared,eta):
   return (eta/(1+x_squared) \
         +(1-eta)*np.exp(-np.log(2)*x_squared))

def gaussian(x_squared,eta=None):
   return np.exp(-np.log(2)*x_squared)

def lorentz(x_squared,eta=None):
   return 1/(1+x_squared)

profiles = {
   'PV': pseudo_voigt,
   'Lorentz': lorentz,
   'Gaussian': gaussian,
}