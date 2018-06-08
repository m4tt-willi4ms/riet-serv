from collections import OrderedDict
import numpy as np

from refinement_parameters import RefinementParameters

DEFAULT_BKGD_ORDER = 3
DEFAULT_TWO_THETA_0 = ('two_theta_0', 0.0, [True], -0.1, 0.1)

class GlobalParameters(RefinementParameters):
    '''
    A class used to keep track of global parameters used in computing powder
    diffraction profiles.
    '''
    def __init__(self, phase_settings,
        two_theta_0=DEFAULT_TWO_THETA_0,
        bkgd_order=DEFAULT_BKGD_ORDER,
        param_dict=None):
        self.phase_settings = phase_settings
        self.bkgd_order = bkgd_order
        if param_dict is None:
            self.two_theta_0 = two_theta_0
            self.bkgd = [x for x in self.bkgd_param_gen()]
        super(GlobalParameters, self).__init__(param_dict=param_dict)

    def set_bkgd_order(self, order):
        r'''
        This method sets the order of the background polynomial, `bkgd`.

        (Strictly speaking, `order` is the number of parameters :math:`c_i` in
        the polynomial

        .. math:: P(2\theta) = \sum_{i=0}^{N} c_i \, (2\theta)^i

        not its degree.)

        Parameters
        ----------
        order : int
            The number of coefficients (:math:`N`) appearing in :math:`P(2\theta)`.

        Returns
        -------
        bkgd : np.array (custom dtype)
            A numpy array containing the coefficients :math:`c_i of the background
            polynomial.
        '''

        self.bkgd_order = self.validate_order(order,
            max_polynom_order=self.phase_settings["max_polynom_order"])
        self.bkgd = [x for x in self.bkgd_param_gen()]

    def bkgd_param_gen(self):
        n = 0
        while n < self.bkgd_order:
            # if cls.bkgd == None:
            yield ('bkgd_'+str(n), 0.0, [True], -float('inf'), float('inf'))
            # else:
            #    yield cls.bkgd[n]
            n += 1

    def param_gen(self):
        d = OrderedDict()
        d['two_theta_0'] = self.two_theta_0
        d['bkgd'] = self.bkgd
        return d.iteritems()

if __name__ == "__main__":
    phase_settings = {}
    phase_settings["max_polynom_order"] = 5
    t = GlobalParameters(phase_settings)
    t.assemble_x()