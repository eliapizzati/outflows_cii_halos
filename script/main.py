"""
Created on Sun Mar  7 17:16:03 2021

@author: anna
"""



import os
import numpy as np
import mydir
import natconst as nc

from post_sol_modules import get_ionization_states, \
                          get_surface_density, get_intensity_raw, get_intensity_convolved

from model_classes import load_from_file
from load_data import obs_data_list
import plot_config as pltc

import time




"""
PARAMETER DESCRIPTION

 - PARAMS
 
     - SFR: star formation rate in Msun/yr
     
     - beta: beta parameter
     
     - f_esc: escape fraction of ionizing photons

     - v_c: circular velocity in km/s
     
     - Zeta: metallicity (solar unity)
     
     - alfa: alpha parameter
     
     - R_in: inner radius in kpc
 
 - PARAMS_OBS
 
     - line_frequency: frequency of the rest-frame CII line in Hz
         
     - redshift: redshift of the CII line
     
     - line_FWHM: FWHM of the CII line
         
     - sersic_effective_radius: effective radius in kpc
     
     - sersic_index: sersic index
     
     - exp_effective_radius: effective radius in kpc
 
=======================
 
WORKFLOW:
    
    profiles = get_profiles(params)
    
    ionization_state = get_ionization_states(profiles, params)
    
    
"""

load_sol_from_file = False


betas = np.asarray([0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4])
#betas = [3.0,3.2,3.4]


datas = obs_data_list[1:4]


for data in datas:

    params = dict([("beta", 1.0), 
               ("SFR", data.params_obs["SFR"]),
               ("f_esc", 0.), 
               ("v_c", data.params_obs["v_c"]),
               ("Zeta", 1.0),
               ("alfa", 1.0),
               ("R_in", 0.3)])    
    
    time_start = time.perf_counter()

    for beta_el in betas:
    
        params.update(beta = beta_el)
    
        if load_sol_from_file == False:
        
            from sol_modules import get_profiles
            profiles = get_profiles(params, resol=1000)
        
            profiles.to_file()
            profiles.plot(ax=pltc.axs_sol)
            print(profiles.check_nans())
        
        else:
        
            profiles = load_from_file(params, type_class = "sol_profiles")
         
            ionization_state = get_ionization_states(profiles, params)
        
            sigma_CII = get_surface_density(profiles, ionization_state, params)
        
            intensity_raw = get_intensity_raw(sigma_CII, params, params_obs)
        
            intensity_conv = get_intensity_convolved(intensity_raw, params, data.params_obs, data)
         
            ionization_state.to_file()
            ionization_state.plot(ax=pltc.axs_ion)
             
            sigma_CII.to_file()
            sigma_CII.plot(ax=pltc.ax_sigma)
            
            intensity_raw.to_file()
            intensity_raw.plot(ax=pltc.ax_int_raw)
            
            intensity_conv.to_file()
            intensity_conv.plot(ax=pltc.ax_int_conv)
            


    time_elapsed = (time.perf_counter() - time_start)

    print("total time elapsed (s)=", time_elapsed)


