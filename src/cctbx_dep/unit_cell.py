'''
This module groups together some custom methods used for masking, updating
the unit-cell parameters.
'''
import numpy as np

from cctbx import uctbx

def update_unit_cell(phase_settings, lattice_parameters):
    '''
        This function takes in a phase_settings dict containing the keywords
        `crystal_system` (and possibly `crystal_system_trigonal`), `uc_mask` ---
        as well as a numpy array of lattice parameters --- and updates the
        `unit_cell` phase_setting according to the values in `lattice_parameters`.

    '''
    crystal_system = phase_settings["crystal_system"]
    uc_mask = phase_settings["uc_mask"]

    if np.char.startswith(crystal_system, "Tric"):
        unit_cell = uctbx.unit_cell(
            (float(x) for x in np.nditer(lattice_parameters)))
    elif np.char.startswith(crystal_system, "M"):
        a = lattice_parameters[0]
        b = lattice_parameters[1]
        c = lattice_parameters[2]
        if uc_mask[3]:
            alpha = lattice_parameters[3]
        else: alpha = 90
        if uc_mask[4]:
            beta = lattice_parameters[3]
        else: beta = 90
        if uc_mask[5]:
            gamma = lattice_parameters[3]
        else: gamma = 90
        unit_cell = uctbx.unit_cell((a, b, c, alpha, beta, gamma))
    elif np.char.startswith(crystal_system, "O"):
        a = lattice_parameters[0]
        b = lattice_parameters[1]
        c = lattice_parameters[2]
        unit_cell = uctbx.unit_cell((a, b, c, 90, 90, 90))
    elif np.char.startswith(crystal_system, "Te"):
        a = lattice_parameters[0]
        if uc_mask[1]:
            b = lattice_parameters[1]
        else: b = a
        if uc_mask[2]:
            c = lattice_parameters[1]
        else: c = a
        unit_cell = uctbx.unit_cell((a, b, c, 90, 90, 90))
    elif np.char.startswith(crystal_system, "Tr"):
        if phase_settings["crystal_system_trigonal"] == "R":
            a = lattice_parameters[0]
            alpha = lattice_parameters[1]
            unit_cell = uctbx.unit_cell((a, a, a, alpha, alpha, alpha))
        elif phase_settings["crystal_system_trigonal"] == "H":
            a = lattice_parameters[0]
            c = lattice_parameters[1]
            unit_cell = uctbx.unit_cell((a, a, c, 90, 90, 120))
    elif np.char.startswith(crystal_system, "He"):
        a = lattice_parameters[0]
        c = lattice_parameters[1]
        unit_cell = uctbx.unit_cell((a, a, c, 90, 90, 120))
    elif np.char.startswith(crystal_system, "C"):
        a = lattice_parameters[0]
        unit_cell = uctbx.unit_cell((a, a, a, 90, 90, 90))

    phase_settings["unit_cell"] = unit_cell

def get_unit_cell_mask(phase_settings):
    '''
        Given a phase_settings dictionary with the entries `unit_cell` and
        `crystal_system`, a unit-cell mask is generated to filter only those
        unit-cell parameters which are independently-refinable.

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
    """returns all unit cell parameters specified by the mask
    (a list of six booleans)
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
    uc_mask = get_unit_cell_mask(phase_settings)
    return [x for x in unit_cell_parameter_gen(phase_settings)]
