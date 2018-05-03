import numpy as np

DEFAULT_BKGD_ORDER = 3
DEFAULT_TWO_THETA_0 = ('two_theta_0', 0.0, True, -0.1, 0.1)
DEFAULT_VERTICAL_OFFSET = False #:False = angular offset; True = Vertical Offset

class GlobalParameters:
    '''
    A class used to keep track of global parameters used in computing powder
    profiles.
    '''
    def __init__(self, max_polynom_order=5):
        self.max_polynom_order = max_polynom_order
        self.bkgd_order = DEFAULT_BKGD_ORDER
        self.two_theta_0 = DEFAULT_TWO_THETA_0
        self.bkgd = [x for x i bkgd_param_gen()]
        self.vertical_offset = DEFAULT_VERTICAL_OFFSET

   def set_bkgd_order(self, order, validate_order_func):
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

        self.bkgd_order = validate_order_func(order)
        self.bkgd = [x for x in self.bkgd_param_gen(order)]

    def bkgd_param_gen(self):
        n = 0
        while n < self.bkgd_order:
            # if cls.bkgd == None:
            yield ('bkgd_'+str(n), 0.0, False, -float('inf'), float('inf'))
            # else:
            #    yield cls.bkgd[n]
            n += 1

    @classmethod
    def assemble_global_x(self):
        #TODO: generate dynamically
        cls.global_x = np.hstack((x for x in cls.global_param_gen()))
        cls.global_x_no_bkgd_mask = np.invert(
            np.char.startswith(cls.global_x['labels'], 'bkgd'))

        cls.two_theta_0 = cls.global_x[0]
        cls.bkgd = cls.global_x[1:1+cls.bkgd.shape[0]]

        @classmethod
    def global_param_gen(cls):
        yield cls.two_theta_0
        yield cls.bkgd