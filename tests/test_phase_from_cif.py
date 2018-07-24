import numpy as np
import pytest
import os

import src.cctbx_dep.phase_from_cif as pfc
import src.cctbx_dep.unit_cell as unit_cell
import src.cctbx_dep.target_wavelengths as tw
from src.rietveld_phases import RietveldPhases as Rp

uc = (4.7605, 4.7605, 12.9956, 90, 90, 120)

@pytest.fixture(scope="module")
def phase_settings():
    phase_settings = {
        "cif_path": os.path.join(os.path.split(__file__)[0],
        "../data/cifs/1000032.cif"),
        "min_two_theta": 10,
        "max_two_theta": 90,
        "intensity_cutoff": 0.005,
        }
    Rp.set_wavelengths(tw.set_wavelength('Cu', 0), phase_settings)
    return phase_settings


@pytest.fixture(scope="module")
def phase_data():
    return {}

def test_load_cif(phase_settings):
    #load_cif() is called in RietveldPhases' __init__
    pfc.load_cif(phase_settings)
    test_phase_settings = phase_settings
    # test_phase = Rp(test_file_path)
    # test_phase_settings = test_phase.phase_settings
    assert 'structure' in test_phase_settings
    labels = ['Al1','O1']
    element_symbols = ['Al','O']
    fps = [0.204601, 0.046428]
    fdps = [0.245574, 0.032228]
    occupancy = 1.0
    for x, l, s, fp, fdp in zip(test_phase_settings["structure"].scatterers(),
              labels, element_symbols, fps, fdps):
        assert x.label == l
        assert np.isclose(x.fp, fp)
        assert np.isclose(x.fdp, fdp)
        assert np.isclose(x.occupancy,occupancy)
        assert x.element_symbol() == s
        print(x.report_details(test_phase_settings["unit_cell"],''))

    assert 'unit_cell' in test_phase_settings
    print test_phase_settings["unit_cell"].parameters()
    for x,y in zip(test_phase_settings["unit_cell"].parameters(), uc):
        assert np.isclose(x,y)

    assert 'crystal_system' in test_phase_settings
    assert type(test_phase_settings["crystal_system"]) is str
    assert test_phase_settings["crystal_system"] == 'Trigonal'

def test_compute_relative_intensities(phase_settings, phase_data):
    phase_data.update(pfc.compute_relative_intensities(phase_settings))
    print(phase_data["weighted_intensities"][0:10])
    assert np.all(np.isclose(phase_data["weighted_intensities"][0:10],
    np.array(
        [[ 0.19601307],
         [ 0.61667176],
         [ 0.31832277],
         [ 1.09714791],
         [ 0.01885168],
         [ 0.8546938 ],
         [ 2.07552284],
         [ 0.05351527],
         [ 0.08408237],
         [ 0.18298366]])))


def test_set_two_theta_peaks(phase_settings, phase_data):
    pfc.load_cif(phase_settings)
    pfc.compute_relative_intensities(phase_settings)
    pfc.set_two_theta_peaks(phase_settings, phase_data)
    assert np.all(np.isclose(phase_data["two_theta_peaks"][0:10], np.array(
    [[ 25.56751716],
    [ 35.13879293],
    [ 37.76313851],
    [ 43.33853666],
    [ 46.16161879],
    [ 52.53263976],
    [ 57.47974186],
    [ 59.71809269],
    [ 61.10698338],
    [ 61.28247987]])))
