import numpy as np
import __init__
from cctbx.miller import sym_equiv_indices

def get_pref_orient_params(phase_settings, params):
    method = phase_settings['pref_orient_method']
    if method == 'march_dollase':
        return params
    elif method == 'sph_harm':
        return None

def get_sym_equiv_angles(phase_settings, phase_data):
    sg = phase_settings['structure'].space_group()
    hkl = phase_settings['pref_orient_hkl']
    uc_rec = phase_settings['unit_cell'].reciprocal()
    angles = []
    for index in phase_data['f_miller_indices']:
        sei = sym_equiv_indices(sg, map(int,list(index)))
        index_angles = []
        for sym_equiv_index in sei.indices():
            index_angles.append(uc_rec.angle(
                hkl, (0, 0, 0), sym_equiv_index.h()))
        angles.append(np.array(index_angles))
        phase_data['md_cos_factors'] = [np.cos(np.pi/180*sym_equiv_angles)**2
            for sym_equiv_angles in angles]
    return angles

def pref_orient_function(r, factors):
    # tmp = np.apply_along_axis(
    #     lambda theta, r: np.power(
    #         1/r + (r**2-1/r)*np.cos(np.pi/180*theta)**2, -1.5),
    #     0, sym_equiv_angles, r)
    tmp = np.power(1/r + (r**2-1/r)*factors, -1.5)
    return np.sum(tmp) / len(factors)

def update_pref_orient_factors(phase_settings, phase_data, pref_or):
    if phase_settings['pref_orient_method'] == 'march_dollase':
        factors = phase_data['md_cos_factors']
        result = []
        for factor_array in factors:
            result.append(pref_orient_function(pref_or[0], factor_array))
        result = np.tile(
            np.array(result), len(phase_settings['wavelengths']))[:, np.newaxis]
        phase_data['pref_orient_factors'] = result
        return result

def set_pref_orient_settings(phase_settings, phase_dict):
    phase_settings['pref_orient_hkl'] \
        = phase_dict.get('pref_orient_hkl', [0, 0, 1])
    phase_settings['pref_orient_method'] \
        = phase_dict.get('pref_orient_method', 'march_dollase')
    phase_settings['pref_orient_ell'] \
        = phase_dict.get('pref_orient_ell', 2)
