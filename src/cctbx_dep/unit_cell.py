'''
This module groups together some custom methods used for masking, updating
the unit-cell parameters.
'''
from __future__ import division, print_function, absolute_import
import numpy as np

from cctbx import uctbx

def update_unit_cell(phase_settings, lattice_parameters):
    '''
        This function takes in a phase_settings dict containing the keywords
        `crystal_system` (and possibly `crystal_system_trigonal`), `uc_mask` ---
        as well as a numpy array of lattice parameters --- and updates the
        `unit_cell` phase_setting according to the values in `lattice_parameters`.

    '''
    # crystal_system = phase_settings["crystal_system"]
    # uc_mask = phase_settings["uc_mask"]

    inverse_filter = get_inverse_filter(phase_settings)
    uc_params = []
    for entry in inverse_filter:
        if entry < 6:
            uc_params.append(lattice_parameters[entry])
        else:
            uc_params.append(entry)

    assert len(uc_params) == 6
    unit_cell = uctbx.unit_cell(uc_params)
    phase_settings["unit_cell"] = unit_cell

def get_inverse_filter(phase_settings):
    '''
        Given a phase_settings dictionary with the `crystal_system` and
        `uc_mask` keys, an
        inverse filter is generated which will take the refinable parameters and
        map them to a full unit-cell list. For example, for a Hexagonal crystal
        system, get_inverse_filter() will return the list::

            [0, 0, 1, 90, 90, 120] .

        (It should be understood that any indices not in the range(0,6) will be
        interpreted as angles.)
    '''
    crystal_system = phase_settings["crystal_system"]
    uc_mask = phase_settings["uc_mask"]

    # n=0
    # inverse_filter = []
    # for i in xrange(6):
    #     if uc_mask[i]:
    #         inverse_filter.append(n)
    #         n += 1
    #     elif i < 3:
    #         inverse_filter.append(n)
    #     else:
    #         inverse_filter.append(90)


    if np.char.startswith(crystal_system, "M"):
        inverse_filter = [0, 1, 2]
        n = 3
        for i in xrange(3,6):
            if uc_mask[i]:
                inverse_filter.append(n)
                n += 1
            else:
                inverse_filter.append(90)
    elif np.char.startswith(crystal_system, "O"):
        inverse_filter = [0, 1, 2, 90, 90, 90]
    elif np.char.startswith(crystal_system, "Te"):
        inverse_filter = [0]
        for i in xrange(1,3):
            if uc_mask[i]:
                inverse_filter.append(0)
        else:
                inverse_filter.append(1)
        inverse_filter.extend([90, 90, 90])
    elif np.char.startswith(crystal_system, "Trig"):
        if phase_settings["crystal_system_trigonal"] == "R":
            inverse_filter = [0, 0, 0, 1, 1, 1]
        elif phase_settings["crystal_system_trigonal"] == "H":
            inverse_filter = [0, 0, 1, 90, 90, 120]
    elif np.char.startswith(crystal_system, "He"):
        inverse_filter = [0, 0, 1, 90, 90, 120]
    elif np.char.startswith(crystal_system, "C"):
        inverse_filter = [0, 0, 0, 90, 90, 90]
    else:
        inverse_filter = [0, 1, 2, 3, 4, 5]

    return inverse_filter

def set_unit_cell_mask(phase_settings):
    '''
        Given a phase_settings dictionary with the entries `unit_cell` and
        `crystal_system`, a unit-cell mask is generated to filter only those
        unit-cell parameters which are independently refinable.

    '''
    uc_params = phase_settings["unit_cell"].parameters()
    crystal_system = phase_settings["crystal_system"]
    if crystal_system == "Triclinic":
        uc_mask = [True, True, True, True, True, True]
    elif crystal_system == "Monoclinic":
        uc_mask = [True, True, True]
        for i in xrange(3, 6):
            if np.isclose(uc_params[i], 90):
                uc_mask.append(False)
            else: uc_mask.append(True)
    elif crystal_system == "Orthorhombic":
        uc_mask = [True, True, True, False, False, False]
    elif crystal_system == "Tetragonal":
        uc_mask = [True]
        if np.isclose(uc_params[1], uc_params[0]):
            uc_mask += [False, True]
        else: uc_mask += [True, False]
        uc_mask += [False, False, False]
    elif crystal_system == "Trigonal":
        if np.isclose(uc_params[3], uc_params[4]) and \
            np.isclose(uc_params[3], uc_params[5]):
            uc_mask = [True, False, False, True, False, False]
            phase_settings["crystal_system_trigonal"] = "R"
        else:
            uc_mask = [True, False, True, False, False, False]
            phase_settings["crystal_system_trigonal"] = "H"
    elif crystal_system == "Hexagonal":
        uc_mask = [True, False, True, False, False, False]
    elif crystal_system == "Cubic":
        uc_mask = [True, False, False, False, False, False]

    assert len(uc_mask) == 6
    phase_settings["uc_mask"] = uc_mask

    return uc_mask

def unit_cell_parameter_gen(phase_settings, uc_mask=None):
    """Returns all unit cell parameters specified by the mask
    (a list of six booleans).
    """
    uc_params = phase_settings["unit_cell"].parameters()
    if uc_mask is None:
        uc_mask = phase_settings["uc_mask"]
    lattice_dev = phase_settings["lattice_dev"]
    uc_labels = ["a", "b", "c", "alpha", "beta", "gamma"]
    for i in xrange(6):
        if uc_mask[i]:
            yield ('uc_'+uc_labels[i],
                     uc_params[i],
                     [True],
                     uc_params[i]*(1-lattice_dev[i]),
                     uc_params[i]*(1+lattice_dev[i])
                    )

def assemble_lattice_parameters(phase_settings):
    """Sets the independent self.lattice_parameters attributes,
    according to the crystal system.

    """
    set_unit_cell_mask(phase_settings)
    return [x for x in unit_cell_parameter_gen(phase_settings)]
