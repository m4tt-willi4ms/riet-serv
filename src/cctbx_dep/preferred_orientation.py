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
    return angles

def get_pref_orient_function(phase_settings, phase_data):
    method = phase_settings['pref_orient_method']
    indices = phase_data['f_miller_indices']
    if method == 'march_dollase':
        phase_data['sym_equiv_angles'] = get_sym_equiv_angles(
            phase_settings, phase_data)
        def md_function(r):
            pass
        return md_function

def _compute_md_coefficients(sg, uc_rec, h0, h):
    result = []
    for index in list(h):
        pass

def pref_orient_function(r, sym_equiv_angles):
    tmp = np.apply_along_axis(
        lambda y, x: np.power(1/x+(x**2-1/x)*np.cos(np.pi/180*y)**2, -1.5), 0,
        sym_equiv_angles, r)
    return np.sum(tmp) / len(sym_equiv_angles)

def update_pref_orient_factors(phase_data, r):
    angles = phase_data['sym_equiv_angles']
    result = []
    for sym_equiv_angles in angles:
        result.append(pref_orient_function(r, sym_equiv_angles))
    result = np.array(result)
    phase_data['pref_orient_factors'] = result
    return result

# h0 = (2, 0, 2)
# uc_rec = uc.reciprocal()
# print uc_rec.parameters()
# for x in sei.indices():
#    ang = uc_rec.angle(h0, (0, 0, 0), x.h())
#    print x.h(), ang
