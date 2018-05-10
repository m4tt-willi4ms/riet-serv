import numpy as np
import copy

keys = ("labels", "values", "refine", "l_limits", "u_limits")

def validate_order(order, max_polynom_order=5):
    print max_polynom_order
    assert isinstance(order, int)
    if order > max_polynom_order:
        order = max_polynom_order
    elif order < 1:
        order = 1
    return order

def make_x_dict(parameters):
    x = {}
    def add_to_x(key, value):
        # for key, value in zip(keys, param):
        if key not in x:
            x[key] = [value]
        else:
            x[key].append(value)

    for param in parameters:
        if isinstance(param, tuple):
            for key, value in zip(keys, param):
                add_to_x(key, value)
        elif isinstance(param, list):
            for item in param:
                for key, value in zip(keys, item):
                    add_to_x(key, value)
        else: raise TypeError('Parameters must be either a parameter-tuple \
or a list of parameter-tuples.')

    for key in x.keys():
        x[key] = np.array(x[key])

    return x

class RefinementParameters:

    def __init__(self):
        self.x = {}
        self.validate_order = validate_order

    def assemble_x(self):
        '''
            Should not be called unless `self.param_gen()` is implemented.
        '''
        param_names = []
        parameters = []
        for name, param in self.param_gen():
            param_names.append(name)
            parameters.append(copy.deepcopy(param))

        self.x = make_x_dict(parameters)

        self.param_names = param_names
        self.parameters = parameters
        def param_gen():
            for name, param in zip(self.param_names, self.parameters):
                yield name, param

        self.param_gen = param_gen

        n = 0
        for name in self.param_names:
            l = np.sum(np.char.startswith(self.x['labels'], name))
            setattr(self, name, self.x['values'][n:n+l])
            n += l


    def update_x(self, x_new, mask=None):
        if mask is None:
            self.x['values'] = x_new
        else:
            self.x['values'][mask] = x_new[mask]

    def reset_x(self):
        self.x = {}
        for name, param in self.param_gen():
            setattr(self, name, param)

    def param_gen(self):
        """
            Must be implemented by children classes to return key, value pairs
            for each set of refinement parameters, where the key is the
            collective name for the parameters, and the value is either a tuple
            or list of tuples defining the individual parameter values.
        """
        raise NotImplementedError

    # def get_param_dict(self, param_tuple):
    #     d = {}
    #     for k, v in zip(keys, param_tuple):
    #         d[k[:-1]] = v
    #     return d
