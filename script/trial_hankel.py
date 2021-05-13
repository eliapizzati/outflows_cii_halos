# -*- coding: utf-8 -*-
"""
Created on Thu May 13 15:51:20 2021

@author: anna
"""


from pyhank import qdht
import numpy as np
import scipy.special
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

import scipy.special as scipy_bessel
import time

import natconst as nc

from post_sol_modules import get_ionization_states, get_surface_density, get_intensity_raw, \
                      get_intensity_convolved, get_chi2

from model_classes import load_from_file


from load_data import obs_data_list, names, names_CII_halo, names_wo_CII_halo,\
                      names_other,  observational_data_fuji




#
                      
def _spline(x0, y0, x):
    f = interp1d(x0, y0, axis=0, fill_value='extrapolate', kind='cubic')
    return f(x)


def hankel_transform(radial_grid, f):
                      
    max_radius = np.max(radial_grid)
    n_points = radial_grid.size
    
    alpha = scipy_bessel.jn_zeros(0, n_points + 1)
    alpha = alpha[0:-1]
    alpha_n1 = alpha[-1]
    S = alpha_n1

    v_max = alpha_n1 / (2 * np.pi * max_radius)    

    r = alpha * max_radius / alpha_n1
    v = alpha / (2 * np.pi * max_radius)
    kr = 2 * np.pi * v      
    
    fr = _spline(radial_grid, f, r)  
    
    jp = scipy_bessel.jv(0, (alpha[:, np.newaxis] @ alpha[np.newaxis, :]) / S)
    jp1 = np.abs(scipy_bessel.jv( 1, alpha))
    T = 2 * jp / ((jp1[:, np.newaxis] @ jp1[np.newaxis, :]) * S)
    jr = jp1 / max_radius
    jv = jp1 / v_max

    fv = jv * np.matmul(T, (fr / jr))
                      
    return kr, fv
                      
                      
                      
#############
                      
                      
                      
data = obs_data_list[1]

params = dict([("DM_model", "NFW"),
   ("beta", 5.0), 
   ("SFR", data.params_obs["SFR"]),
   ("f_esc_ion", 0.), 
   ("f_esc_FUV", 0.), 
   ("v_c", data.params_obs["v_c"]),
   ("M_vir", data.params_obs["M_vir"]),
   ("redshift", data.params_obs["redshift"]),
   ("Zeta", 1.0),
   ("alfa", 1.0),
   ("R_in", 0.3)])    

profiles = load_from_file(params, class_type = "profiles")
 
ionization_state = get_ionization_states(profiles, params)

sigma_CII = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500, add_CMB_suppression=True)

intensity_raw = get_intensity_raw(sigma_CII, params, data.params_obs)

beam_interp = np.interp(intensity_raw.h/1e3/nc.pc, data.x_beam/1e3/nc.pc, data.beam, right=0.)

beam_interp[beam_interp<0.] = 0.

t_h = time.perf_counter()

kr_int, h_int = qdht(intensity_raw.h/1e3/nc.pc, intensity_raw.var)
kr_beam, h_beam = qdht(intensity_raw.h/1e3/nc.pc, beam_interp)

r_h, int_conv_h = qdht(kr_int, h_int*h_beam)

int_conv_h /= (2*np.pi)

time_h = time.perf_counter() - t_h

t_h2 = time.perf_counter()

kr_int2, h_int2 = hankel_transform(intensity_raw.h/1e3/nc.pc, intensity_raw.var)
kr_beam2, h_beam2 = hankel_transform(intensity_raw.h/1e3/nc.pc, beam_interp)

r_h2, int_conv_h2 = hankel_transform(kr_int2, h_int2*h_beam2)

int_conv_h2 /= (2*np.pi)

time_h2 = time.perf_counter() - t_h2

t_conv = time.perf_counter()

intensity_conv = get_intensity_convolved(intensity_raw, params, data.params_obs, data, add_central_contribution=False)

time_conv = time.perf_counter() - t_conv

print("time h", time_h)
print("time h2", time_h2)
print("time conv", time_conv)
fig, ax = plt.subplots()

ax.plot(intensity_conv.h/1e3/nc.pc, intensity_conv.var)

ax.plot(intensity_raw.h/1e3/nc.pc, beam_interp)
ax.plot(intensity_raw.h/1e3/nc.pc, intensity_raw.var)

ax.plot(r_h, int_conv_h*intensity_conv.var[0]/int_conv_h[0])
ax.set_yscale("log")
#ax.set_xlim(0.,15.)

#ax.set_ylim(-5,40.)
