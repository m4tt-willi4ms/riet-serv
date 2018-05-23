import numpy as np
import copy

keys = ("labels", "values", "refine", "l_limits", "u_limits")
param_keys = ("name", "value", "uround", "l_limit", "u_limit")
ignored_keys = ("lattice_parameters", "scale")

def validate_order(order, max_polynom_order=5):
    assert isinstance(order, int)
    if order > max_polynom_order:
        order = max_polynom_order
    elif order < 1:
        order = 1
    return order

def make_x_dict(parameters):
    x = {}

    def add_to_x(param):
        # for key, value in zip(keys, param):
        for key, value in zip(keys, param):
            if key not in x:
                x[key] = [value]
            else:
                x[key].append(value)

    for param in parameters:
        if isinstance(param, tuple):
            add_to_x(param)
        elif isinstance(param, list):
            for item in param:
                add_to_x(item)
        else: raise TypeError('Parameters must be either a parameter-tuple \
or a list of parameter-tuples.')

    for key in x.keys():
        x[key] = np.array([np.array(x_i) for x_i in x[key]])

    return x

class RefinementParameters(object):

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
        for name, param in zip(param_names, parameters):
            if isinstance(param, list):
                l = len(param)
            else: l = 1
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

    def set_param(self, param_name, val_tuple):
        assert len(val_tuple) == len(keys)
        setattr(self, param_name, val_tuple)

    def as_dict(self):
        result = {}
        if self.x:
            n=0
            def make_dict_from_x(index):
                lst = [self.x[k][index] for k in keys]
                lst[2] = [bool(x) for x in np.nditer(lst[2])]
                return dict(zip(param_keys, lst))
            for name, param in self.param_gen():
                val_array = getattr(self, name)
                l = len(val_array)
                if l > 1:
                    value = []
                    for i, val in enumerate(np.nditer(val_array)):
                        value.append(make_dict_from_x(n+i))
                elif l == 1:
                    value = make_dict_from_x(n)
                result[name] = value
                n += l
        else:
            for name, param in self.param_gen():
                result[name] = param
        return result


    def from_dict(self, d):
        def get_param_from_dict(d):
            return tuple([d[key] for key in param_keys])

        def get_param(val):
            result = None
            if isinstance(val, list):
                result = []
                for item in val:
                    result.append(get_param_from_dict(item))
            elif isinstance(val, dict):
                result = get_param_from_dict(val)
            return result

        for name, param in self.param_gen():
            if name in filter(lambda x: not x in ignored_keys, d.keys()):
                setattr(self, name, get_param(d[name]))
