import numpy as np

keys = ("labels", "values", "mask", "l_limits", "u_limits")

class RefinementParameters:
    x = None

    def assemble_x(self):
        x = {}
        parameters = []
        for param in self.param_gen():
            if type(param) == tuple:
                parameters.append(param)
            elif type(param) == list:
                for item in param:
                    parameters.append(item)

        for param in parameters:
            for key, value in zip(keys, param):
                if key not in x:
                    x[key] = [value]
                else:
                    x[key].append(value)

        for key in x.keys():
            x[key] = np.array(x[key])

        self.x = x

        return x

    def update_x(self, x_new, mask=None):
        if mask is None:
            self.x['values'] = x_new
        else:
            self.x['values'][mask] = x_new[mask]

    def param_gen(self):
        raise NotImplementedError
