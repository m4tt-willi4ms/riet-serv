import numpy as np

def get_pref_orient_params(phase_settings, params):
    method = phase_settings['pref_orient_method']
    if method == 'march_dollase':
        return params

def get_pref_orient_coefficients(phase_settings, phase_data, params):
    method = phase_settings['pref_orient_method']
    hkl = phase_settings['pref_orient_hkl']
    sg = phase_settings['structure'].space_group()
    
    if method == 'march_dollase':
