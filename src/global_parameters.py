from refinement_parameter import refinement_parameter

DEFAULT_BKGD_ORDER = 3
DEFAULT_TWO_THETA_0 = refinement_parameter('two_theta_0', 0.0, True, -0.1, 0.1)

class GlobalParameters:
   self.bkgd_order = DEFAULT_BKGD_ORDER
   self.two_theta_0 = DEFAULT_TWO_THETA_0

   def __init__(self):
