# -*- coding: utf-8 -*-
"""
Created on Wed May 12 14:35:50 2021

@author: anna
"""

import itertools
import os
import time

import matplotlib.pyplot as plt

import my_dir
import natconst as nc
import plot_config as pltc
from mcmc import data, other_params, get_emission_fast, h, f_beam, grid
from post_sol_modules import get_ionization_states, get_surface_density, get_intensity_raw, get_intensity_convolved
from sol_modules import get_profiles

betas = [7.5]
SFRs = [50.]
v_cs = [200.]



fig_int_conv, ax_int_conv = pltc.plot_configurator(plot_type="int", xlim=15)

ax_int_conv.set_ylim((1e-3,1e2))

folder = "plot_fast_emission"
    

if not os.path.exists(os.path.join(my_dir.plot_dir, folder)):
    os.mkdir(os.path.join(my_dir.plot_dir, folder))

for  beta, SFR, v_c in itertools.product(betas, SFRs, v_cs):

    theta = [beta, SFR,v_c]
    print(theta)
     
    time_profile = time.perf_counter()
        
    intensity = get_emission_fast(theta, data, other_params, h, grid, f_beam)    
            
    time_profile = (time.perf_counter() - time_profile)
        
    print("total emission time (s)=", time_profile)
        
    ax_int_conv.plot(h, intensity)
    
    
    params = dict([("DM_model", "NFW"),
           ("beta", 5.5), 
           ("SFR", SFR),
           ("f_esc_ion", 0.), 
           ("f_esc_FUV", 0.), 
           ("v_c", v_c),
           ("redshift", data.params_obs["redshift"]),
           ("Zeta", 1.0),
           ("alfa", 1.0),
           ("R_in", 0.3)])   
    
    profiles = get_profiles(params, resol=500)
    #profiles.plot(savefig=True)
    ionization_state = get_ionization_states(profiles, params)
    #ionization_state.plot(savefig=True)
    sigma_CII = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500, add_CMB_suppression=True)
    intensity_raw = get_intensity_raw(sigma_CII, params, data.params_obs)
    #intensity_raw.plot(savefig=True)
    intensity_conv = get_intensity_convolved(intensity_raw, params, data.params_obs, data, add_central_contribution=False)

    
    ax_int_conv.plot(intensity_conv.h/(1000*nc.pc), intensity_conv.var)

    
alpine = ax_int_conv.errorbar(data.x/(1000*nc.pc), data.data, yerr=data.err, \
        markerfacecolor='maroon',markeredgecolor='maroon', marker='o',\
        linestyle='', ecolor = 'maroon')
                        

   
fig_int_conv.legend(loc="lower center", ncol=8, fontsize="small")
plt.savefig(os.path.join(my_dir.plot_dir, folder, "emission_mus.png"))
    
plt.show()
        
    
    
    
