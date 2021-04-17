# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 17:01:47 2021

@author: anna
"""

"""
PARAMETER DESCRIPTION

 - PARAMS
 
     - SFR: star formation rate in Msun/yr
     
     - beta: beta parameter
     
     - f_esc: escape fraction of ionizing photons

     - v_c: circular velocity in km/s
     
     - redshift: redshift 
     
     - Zeta: metallicity (solar unity)
     
     - alfa: alpha parameter
     
     - R_in: inner radius in kpc
 
 - PARAMS_OBS

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