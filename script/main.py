"""
Created on Sun Mar  7 17:16:03 2021

@author: anna
"""



import os
import numpy as np
import matplotlib.pyplot as plt
import mydir
import natconst as nc

from model_modules import get_profiles, get_ionization_states, get_surface_density, get_intensity_raw, get_intensity_convolved
#from data_modules import *
#from utils import * 
from data_classes import obs_data
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


time_start = time.clock()

profiles = get_profiles(params)
    
ionization_state = get_ionization_states(profiles, params)

sigma_CII = get_surface_density(profiles, ionization_state, params)

intensity_raw = get_intensity_raw(sigma_CII, params, params_obs)

intensity_conv = get_intensity_convolved(intensity_raw, params, params_obs)


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
intensity_conv.plot()

time_elapsed = (time.perf_counter() - time_start)

print(time_elapsed)


