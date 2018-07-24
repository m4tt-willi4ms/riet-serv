from __future__ import division, print_function, absolute_import
import numpy as np
import time
from random import randrange
import profile, pstats
import pytest
import copy
import os

import src.rietveld_phases as rp
import src.phase_parameters as pp
from src.rietveld_phases import RietveldPhases as Rp
from src.rietveld_refinery import RietveldRefinery

import src.cctbx_dep.target_wavelengths

@pytest.fixture(scope="function")
def set_profile():
    return Rp.set_profile(
    r"./data/profiles/Jade-Al2O3-Sim.xye", number_of_columns=2)

@pytest.fixture(scope="module")
def test_phase():
    return Rp("./data/cifs/1000032.cif",
        intensity_cutoff=0.01,
        delta_theta=2.0)

def test_data_read_in(set_profile):
    assert len(Rp.two_theta) == 4250
    assert np.isclose(Rp.phase_settings["d_min"], 1.08934)
    assert np.isclose(Rp.phase_settings["d_max"], 17.63254)
    assert len(Rp.I) == 4250
    assert np.isclose(Rp.two_theta[-1],90)
    assert len(Rp.sigma) == 4250
    assert np.isclose(Rp.sigma[-1],np.sqrt(Rp.I[-1]))

def test_global_parameters_exist():
    Rp_dict = Rp.global_parameters.__dict__
    assert 'two_theta_offset' in Rp_dict
    assert 'bkgd' in Rp_dict

def test_parameters_exist(test_phase):
    tp_dict = test_phase.phase_parameters.__dict__
    print(tp_dict.keys())
    assert 'U' in tp_dict
    assert 'V' in tp_dict
    assert 'W' in tp_dict
    assert 'scale' in tp_dict
    assert 'eta' in tp_dict
    assert 'lattice_parameters' in tp_dict

def test_set_bkgd_order(test_phase):
    Rp.global_parameters.set_bkgd_order(3)
    Rp.assemble_global_x()
    assert len(Rp.bkgd) == 3
    assert len(test_phase.bkgd) == 3

def test_bkgd_param_gen():
    Rp.global_parameters.set_bkgd_order(3)
    gen = Rp.global_parameters.bkgd_param_gen()
    for n in xrange(0,3):
        assert next(gen)[0] == 'bkgd_' + str(n)
    assert pytest.raises(StopIteration,next,gen)

def test_background_polynomial():
    Rp.assemble_global_x()
    Rp.bkgd[0] = 1
    assert np.all(np.isclose(Rp.background_polynomial(),
        np.ones(len(Rp.two_theta),dtype=float)))

    Rp.bkgd = np.array([1,2,3],dtype=float)
    assert np.all(np.isclose(Rp.background_polynomial(),
        1.0+2.0*Rp.two_theta+3.0*Rp.two_theta_powers[2,:]))

    Rp.bkgd = np.array([0,0,0],dtype=float)
    assert np.all(np.isclose(Rp.background_polynomial(),
        np.zeros(len(Rp.two_theta),dtype=float)))

def test_U_default(test_phase):
    test_U = test_phase.phase_parameters.U
    for x, y in zip(test_U, pp.DEFAULT_U):
        assert x == y
    assert getattr(test_phase.phase_parameters, test_U[0]) == pp.DEFAULT_U[1]

def test_set_vertical_offset(test_phase):
    assert test_phase.phase_settings["vertical_offset"] == False
    assert 'cos_theta' not in test_phase.__dict__
    Rp.phase_settings["vertical_offset"] = True
    Rp.set_profile(r"./data/profiles/Jade-Al2O3-Sim.xye")
    assert 'cos_theta' in Rp.__dict__
    assert np.isclose(Rp.cos_theta[-1],1/np.sqrt(2))

def test_LP_intensity_scaling():
    assert len(Rp.two_theta) == len(Rp.LP_intensity_scaling())

    two_theta = Rp.two_theta
    assert np.all(np.isclose(Rp.LP_intensity_scaling(),
        (1+np.cos(np.pi/180*two_theta)**2) \
            /np.sin(np.pi/360*two_theta) \
            /np.sin(np.pi/180*two_theta)))

def test_assemble_global_x():
    #assemble_global_x() is called in RietveldPhases' __init__
    assert 'global_x' in Rp.__dict__
    assert 'global_x_no_bkgd_mask' in Rp.__dict__
    assert len(Rp.global_x) == len(Rp.bkgd) + 1 # + 1 for two_theta_0

def test_update_global_x():
    mask = np.zeros(len(Rp.global_x),dtype=bool)
    mask[0] = True
    new_x = copy.deepcopy(Rp.global_x)
    new_x[0] = 0.3
    Rp.global_parameters.update_x(new_x,mask)
    assert np.isclose(Rp.global_x[0],0.3)

def test_global_param_gen():
    gen = Rp.global_parameters.param_gen()
    exp_tt_0 = next(gen)
    assert exp_tt_0[0] == 'two_theta_offset'
    assert len(exp_tt_0[1]) == 5
    assert exp_tt_0[1][1] == 0.0
    exp_bkgd = next(gen)
    assert exp_bkgd[0] == 'bkgd'
    assert len(exp_bkgd[1]) == Rp.global_parameters.bkgd_order
    for val in exp_bkgd[1]:
        assert val[1] == 0.0
    assert pytest.raises(StopIteration, next, gen)

def test_phase_param_gen(test_phase):
    count = 0
    for x in test_phase.phase_parameters.param_gen():
        count += 1
    assert count == 6

def test_assemble_phase_x(test_phase):
    #assemble_phase_x() is called in RietveldPhases' __init__
    tp_dict = test_phase.__dict__
    assert 'phase_x' in tp_dict
    assert 'global_mask_no_bkgd' in tp_dict
    assert 'phase_mask' in tp_dict

    x = len(test_phase.global_mask_no_bkgd)
    assert x == len(test_phase.phase_x) + 1
    assert x == len(test_phase.phase_mask)

    assert str(type(test_phase.phase_parameters.cagliotti_u))[7:17] \
        == 'numpy.ndar'
    assert len(test_phase.phase_parameters.U) == 5

    assert len(test_phase.phase_parameters.eta) \
        == test_phase.phase_parameters.eta_order
    assert len(test_phase.phase_parameters.lattice_parameters) == 2

def test_eta_polynomial(test_phase):
    assert 'eta' in test_phase.phase_parameters.__dict__
    test_phase.phase_parameters.set_eta_order(2)
    test_phase.assemble_phase_x()
    test_eta = np.array([0.5,0.005])
    test_phase.phase_parameters.update_x(test_eta,
        np.char.startswith(test_phase.phase_parameters.x['labels'], 'eta'),
        apply_mask_to_input=False)
    tmp = test_phase.eta_polynomial()
    assert np.isclose(tmp[-1],test_eta[0]+Rp.two_theta[-1]*test_eta[1])
    test_phase.phase_parameters.update_x(np.array([0.5,0]),
        np.char.startswith(test_phase.phase_parameters.x['labels'], 'eta'),
        apply_mask_to_input=False)

def test_compute_relative_intensities(test_phase):
    tp_dict = test_phase.phase_data
    assert 'crystal_density' in tp_dict
    assert 'f_miller_set' in tp_dict
    assert 'd_spacings' in tp_dict
    assert 'relative_intensities' in tp_dict

    assert np.all(np.isclose(test_phase.phase_data["d_spacings"], np.array(
        [ 3.48114434,  2.55177292,  2.38025   ,  2.0860757 ,  1.74057217,
          1.60196736,  1.54715716,  1.51527741,  1.51135808,  1.40499663,
          1.37423798,  1.3364588 ,  1.27588646,  1.23944057,  1.23455034,
          1.19354167,  1.190125  ,  1.16038145,  1.16038145,  1.14760201,
          1.12613194,  1.12452033,  1.09932894])))
    print(repr(test_phase.phase_data["relative_intensities"]))
    assert np.all(np.isclose(test_phase.phase_data["relative_intensities"],
        np.array(
        [ 0.19601307,  0.61667176,  0.31832277,  1.09714791,  0.8546938 ,
          2.07552284,  0.05351527,  0.08408237,  0.18298366,  1.14648955,
          1.82420591,  0.03712813,  0.05924289,  0.6978057 ,  0.34759558,
          0.03013175,  0.29775878,  0.0213698 ,  0.02146123,  0.24001883,
          0.18601828,  0.1423521 ,  0.40275275])))

def test_set_profile():
    xyepath = os.path.join(os.path.dirname(__file__),
        '..\\data\\profiles\\cement_15_03_11_0028.xye')
    Rp.set_profile(xyepath, min_two_theta=50, max_two_theta=52)
    assert Rp.two_theta is not None
    print(repr(Rp.two_theta))
    assert np.all(np.isclose(Rp.two_theta, np.array(
        [ 50.0051,  50.0253,  50.0454,  50.0655,  50.0857,  50.1058,
          50.1259,  50.1461,  50.1662,  50.1864,  50.2065,  50.2266,
          50.2468,  50.2669,  50.2871,  50.3072,  50.3273,  50.3475,
          50.3676,  50.3877,  50.4079,  50.428 ,  50.4482,  50.4683,
          50.4884,  50.5086,  50.5287,  50.5489,  50.569 ,  50.5891,
          50.6093,  50.6294,  50.6495,  50.6697,  50.6898,  50.71  ,
          50.7301,  50.7502,  50.7704,  50.7905,  50.8107,  50.8308,
          50.8509,  50.8711,  50.8912,  50.9113,  50.9315,  50.9516,
          50.9718,  50.9919,  51.012 ,  51.0322,  51.0523,  51.0725,
          51.0926,  51.1127,  51.1329,  51.153 ,  51.1731,  51.1933,
          51.2134,  51.2336,  51.2537,  51.2738,  51.294 ,  51.3141,
          51.3343,  51.3544,  51.3745,  51.3947,  51.4148,  51.4349,
          51.4551,  51.4752,  51.4954,  51.5155,  51.5356,  51.5558,
          51.5759,  51.5961,  51.6162,  51.6363,  51.6565,  51.6766,
          51.6967,  51.7169,  51.737 ,  51.7572,  51.7773,  51.7974,
          51.8176,  51.8377,  51.8579,  51.878 ,  51.8981,  51.9183,
          51.9384,  51.9585,  51.9787,  51.9988])))