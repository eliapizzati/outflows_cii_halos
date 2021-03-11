"""
Created on Sun Mar  7 17:16:03 2021

@author: anna
"""



import os
import numpy as np
import matplotlib.pyplot as plt
import mydir
import natconst as nc

from model_modules import get_ionization_states, \
                          get_surface_density, get_intensity_raw, get_intensity_convolved

from model_classes import load_from_file
from load_data import observational_data, params_obs

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
 
     - line_frequency: frequency of the CII line in Hz
         
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

params = dict([("class", "sol"),
               ("type", "n"),
               ("SFR", 20.),
               ("beta", 1.0), 
               ("f_esc", 0.), 
               ("v_c", 175.),
               ("Zeta", 1.0),
               ("alfa", 1.0),
               ("R_in", 0.3)])    

load_sol_from_file = False

    
betas = np.asarray([0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4])

time_start = time.perf_counter()

for beta_el in betas:
    
    params.update(beta = beta_el)
    
    
    if load_sol_from_file == True:
        profiles = load_from_file(params, type_class = "sol_profiles")
    else:
        from model_modules import get_profiles, 
        profiles = get_profiles(params)
        
    ionization_state = get_ionization_states(profiles, params)

    sigma_CII = get_surface_density(profiles, ionization_state, params)

    intensity_raw = get_intensity_raw(sigma_CII, params, params_obs)

    intensity_conv = get_intensity_convolved(intensity_raw, params, params_obs, observational_data)

    
    profiles.to_file()
    profiles.plot()
    
    print(profiles.check_nans())
    
    ionization_state.to_file()
    ionization_state.plot()
    
    
    sigma_CII.to_file()
    sigma_CII.plot()
    
    intensity_raw.to_file()
    intensity_raw.plot()
    
    intensity_conv.to_file()
    intensity_conv.plot(data=observational_data)
    
time_elapsed = (time.perf_counter() - time_start)

print(time_elapsed)


