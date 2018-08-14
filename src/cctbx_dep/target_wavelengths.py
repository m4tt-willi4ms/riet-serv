"""This module makes use of the cctbx.eltbx wavelengths module to populate a
    dictionary of available characteristic wavelengths, and sets the appropriate
    phase_settings keys according the target specified.
"""
from __future__ import division, print_function, absolute_import
import numpy as np

from cctbx.eltbx import wavelengths

WAVELENGTHS_DICT = {
    "CoA1": 1.788965,
    "CoA2": 1.792850,
}
for w in wavelengths.characteristic_iterator():
    WAVELENGTHS_DICT[w.label()] = w.as_angstrom()

TARGETS = tuple(set([label[0:2] for label, val in WAVELENGTHS_DICT.items()]))

K_ALPHA_FACTORS = [1.0, 0.502]

def set_wavelength(
        target='Cu', wavelength_model=0, custom_wavelength=None):
    wavelengths = []
    if isinstance(target, basestring):
        assert target in TARGETS
        assert wavelength_model in (0, 1)
        if wavelength_model == 0:
            wavelengths = [
                WAVELENGTHS_DICT[target+"A1"],
                WAVELENGTHS_DICT[target+"A2"],
                ]
        elif wavelength_model == 1:
            wavelengths = [WAVELENGTHS_DICT[target+"A1"]]
    elif wavelength_model == 2 and custom_wavelength is not None:
        wavelengths = [float(custom_wavelength)]
    return wavelengths
