from __future__ import division, print_function, absolute_import
import numpy as np
import time
from random import randrange
import profile, pstats
import pytest
import copy

import src.rietveld_phases as rp
import src.phase_parameters as pp
from src.rietveld_phases import RietveldPhases as Rp
from src.rietveld_refinery import RietveldRefinery

import src.cctbx_dep.target_wavelengths

Rp.set_profile(
    r"./data/profiles/Jade-Al2O3-Sim.xye", number_of_columns=2)

def test_data_read_in():
    assert len(Rp.two_theta) == 4250
    assert np.isclose(Rp.phase_settings["d_min"], 1.08934)
    assert np.isclose(Rp.phase_settings["d_max"], 17.63254)
    assert len(Rp.I) == 4250
    assert np.isclose(Rp.two_theta[-1],90)
    assert len(Rp.sigma) == 4250
    assert np.isclose(Rp.sigma[-1],np.sqrt(Rp.I[-1]))

def test_global_parameters_exist():
    Rp_dict = Rp.__dict__
    assert 'two_theta_offset' in Rp_dict
    assert 'bkgd' in Rp_dict

test_phase = Rp("./data/cifs/1000032.cif")

def test_parameters_exist():
    tp_dict = test_phase.phase_parameters.__dict__
    print(tp_dict.keys())
    assert 'U' in tp_dict
    assert 'V' in tp_dict
    assert 'W' in tp_dict
    assert 'scale' in tp_dict
    assert 'eta' in tp_dict
    assert 'lattice_parameters' in tp_dict

def test_set_bkgd_order():
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

def test_U_default():
    test_U = test_phase.phase_parameters.U
    for x, y in zip(test_U, pp.DEFAULT_U):
        assert x == y
    assert getattr(test_phase.phase_parameters, test_U[0]) == pp.DEFAULT_U[1]

def test_set_vertical_offset():
    assert test_phase.phase_settings["vertical_offset"] == False
    assert 'cos_theta' not in test_phase.__dict__
    Rp.phase_settings["vertical_offset"] = True
    Rp.set_profile(r"./data/profiles/Jade-Al2O3-Sim.xye")
    assert 'cos_theta' in Rp.__dict__
    assert np.isclose(Rp.cos_theta[-1],-360/np.pi/np.sqrt(2))

def test_LP_intensity_scaling():
    assert len(Rp.two_theta) == len(Rp.LP_intensity_scaling())

    Rp.two_theta_0 = 0.1
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

def test_phase_param_gen():
    count = 0
    for x in test_phase.phase_parameters.param_gen():
        count += 1
    assert count == 6

def test_assemble_phase_x():
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

def test_eta_polynomial():
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

def test_compute_relative_intensities():
    tp_dict = test_phase.phase_data
    assert 'crystal_density' in tp_dict
    assert 'f_miller_set' in tp_dict
    assert 'd_spacings' in tp_dict
    assert 'relative_intensities' in tp_dict

    assert np.all(np.isclose(test_phase.phase_data["d_spacings"], np.array(
        [ 3.48114434,  2.55177292,  2.38025   ,  2.0860757 ,  1.74057217,
          1.60196736,  1.54715716,  1.51527741,  1.51135808,  1.40499663,
          1.37423798,  1.3364588 ,  1.27588646,  1.23944057,  1.23455034,
          1.19354167,  1.190125  ,  1.14760201,  1.12613194,  1.12452033,
          1.09932894])))
    assert np.all(np.isclose(test_phase.phase_data["relative_intensities"],
        np.array(
        [ 0.19906132,  0.61069362,  0.29855677,  1.18211397,  0.86296333,
          2.11656703,  0.05664854,  0.08304719,  0.18443051,  1.14926644,
          1.85168419,  0.03841942,  0.0587068 ,  0.69865033,  0.35418387,
          0.03049309,  0.29644002,  0.24097501,  0.18603062,  0.14235153,
          0.40296412])))

def exercise_RietveldPhases():
    cifs = ["1000032.cif","1507774.cif"]

    # Testing LP_intensity_scaling
    assert np.isclose(Rt[0].LP_Intensity_Scaling(20.0),31.7054214503)

    # Testing two_theta outputs
    assert np.all(np.isclose(Rt[0].two_theta_peaks[0:10], np.array(
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

    #Testing Weighted Intensity outputs
    assert np.all(np.isclose(Rt[0].weighted_intensities[0:10], np.array(
        [[ 0.19906132],
         [ 0.61069362],
         [ 0.29855677],
         [ 1.18211397],
         [ 0.01740687],
         [ 0.86296333],
         [ 2.11656703],
         [ 0.05664854],
         [ 0.08304719],
         [ 0.18443051]])))

    #Testing PseudoVoigtProfile
    if display_plots:
        for RV in Rt:
            #Select a random peak:
            if len(RV.two_theta_peaks) != 0:
                rnd_index = randrange(0,len(RV.two_theta_peaks),1)
                RV.showPVProfilePlot("Test",rnd_index, autohide=False)

    # #Testing Refinery
    # RR = RietveldRefinery(Rt,minimizer_input_string,
    #    store_intermediate_state=True) #,show_plots=False)

    # # RR.display(RR.minimize_Scale_Offset)
    # RR.display(RR.minimize_Scale_Offset_W)
    # RR.display(RR.minimize_All)



def run():
    # exercise_RietveldPhases()
    print("OK")

if (__name__ == "__main__"):
    # pr = profile.Profile()
    # pr.enable()
    # run()
    unittest.main(warnings='ignore')
    # profile.run('run(); print')
    # pr.disable()
    # s = StringIO.StringIO()
    # sortby = 'cumulative'
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # ps.print_stats()
    # print s.getvalue()