import numpy as np

def peak_masks(two_theta, two_theta_0, two_theta_peaks, delta_theta):
   # print "called peak_masks()", inspect.stack()[1][3]
   return np.abs(two_theta - two_theta_0 - two_theta_peaks) < delta_theta

def get_masked_array(array, shape, mask):
   return np.broadcast_to(array, shape)[mask]
