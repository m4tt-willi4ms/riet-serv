import numpy as np

from refinement_parameters import RefinementParameters

DEFAULT_BKGD_ORDER = 3
DEFAULT_TWO_THETA_0 = ('two_theta_0', 0.0, True, -0.1, 0.1)
DEFAULT_VERTICAL_OFFSET = False #:False = angular offset; True = Vertical Offset

class GlobalParameters(RefinementParameters):
    '''
    A class used to keep track of global parameters used in computing powder
    diffraction profiles.
    '''
    def __init__(self, validate_order_func=lambda x: x):
        self.validate_order_func = validate_order_func
        self.bkgd_order = DEFAULT_BKGD_ORDER
        self.two_theta_0 = DEFAULT_TWO_THETA_0
        self.bkgd = [x for x in self.bkgd_param_gen()]
        self.vertical_offset = DEFAULT_VERTICAL_OFFSET

    def set_bkgd_order(self, order):
        r'''
        This method sets the order of the background polynomial, `bkgd`.

        (Strictly speaking, `order` is the number of parameters :math:`c_i` in
        the polynomial

        .. math:: P(2\theta) = \sum_{i=0}^{N} c_i (2\theta)^i

        not its degree.)

        Parameters
        -----------
        order : int
            The number of coefficients (:math:`N`) appearing in :math:`P(2\theta)`.

        Returns
        -------
        bkgd : np.array (custom dtype)
            A numpy array containing the coefficients :math:`c_i of the background
            polynomial.
        '''

        self.bkgd_order = self.validate_order_func(order)
        self.bkgd = [x for x in self.bkgd_param_gen(order)]

    def bkgd_param_gen(self):
        n = 0
        while n < self.bkgd_order:
            # if cls.bkgd == None:
            yield ('bkgd_'+str(n), 0.0, True, -float('inf'), float('inf'))
            # else:
            #    yield cls.bkgd[n]
            n += 1

    def param_gen(self):
        yield self.two_theta_0
        yield self.bkgd

    def set_vertical_offset(self, value):
        assert type(value) == bool
        cls.vertical_offset = value

if __name__ == "__main__":
    t = GlobalParameters()
    t.assemble_x()