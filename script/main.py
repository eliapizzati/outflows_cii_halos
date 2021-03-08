"""
Created on Sun Mar  7 17:16:03 2021

@author: anna
"""



import os
import numpy as np
import matplotlib.pyplot as plt
import mydir
import natconst as nc

from model_modules import profile, get_profiles, get_ionization_states, get_surface_density, get_intensity_raw, get_intensity_convolved
#from data_modules import *
#from utils import * 

import time



"""
PARAMETER DESCRIPTION

 - PARAMS
 
     - class: "sol", "ion", "flux", "cool"
     
     - type: if class = "sol", "v", "n", "T";
              if class = "ion", "x_H", "x_CII";
              else type = None
              
     - SFR: star formation rate in Msun/yr
     
     - beta: beta parameter
     
     - f_esc: escape fraction of ionizing photons

     - v_c: circular velocity in km/s
     
     - Zeta: metallicity (solar unity)
     
     - alfa: alpha parameter
     
     - R_in: inner radius in kpc
 
 - PARAMS_OBS
 
     - line_frequency: frequency of the CII line in Hz
         
     - redshift: redshift of the observed 
     
     - line_FWHM:
         
     - r_beam:
         
     - beam: 
 
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

params_obs = dict([("line_frequency", 1900*1e9),
                   ("redshift", 5.96),
                   ("line_FWHM", 296*nc.km),
                   ("r_beam", None), 
                   ("beam", None)])


time_start = time.clock()

profiles = get_profiles(params)
    
ionization_state = get_ionization_states(r, profiles, params)

#h, sigma_CII = get_surface_density(r, profiles, ionization_state, params)


profiles.to_file()
profiles.plot()

print(profiles.check_nans())

ionization_state.to_file()
ionization_state.plot()


time_elapsed = (time.clock() - time_start)

print(time_elapsed)


