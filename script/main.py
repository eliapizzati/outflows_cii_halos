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


betas = np.asarray([0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4])
#betas = [3.0,3.2,3.4]


datas = obs_data_list


for data in datas:
    
    print("#####################")
    print(data.params_obs["name"])
    print("#####################")


    params = dict([("beta", 1.0), 
               ("SFR", data.params_obs["SFR"]),
               ("f_esc", 0.), 
               ("v_c", data.params_obs["v_c"]),
               ("Zeta", 1.0),
               ("alfa", 1.0),
               ("R_in", 0.3)])    
    
    time_start = time.perf_counter()

    fig_sol, axs_sol = pltc.plot_configurator(plot_type="sol")
        
    fig_ion, axs_ion = pltc.plot_configurator(plot_type="ion")

    fig_int_raw, ax_int_raw = pltc.plot_configurator(plot_type="int")

    fig_int_conv, ax_int_conv = pltc.plot_configurator(plot_type="int")
    
    ax_int_raw.set_ylim((1e-3,1e2))
    ax_int_conv.set_ylim((1e-3,1e2))


    fig_sol.suptitle("{0:}, v_c = {1:.1f} km/h, SFR = {2:.1f}".format(data.params_obs["name"], data.params_obs["v_c"], data.params_obs["SFR"]))
    fig_ion.suptitle("{0:}, v_c = {1:.1f} km/h, SFR = {2:.1f}".format(data.params_obs["name"], data.params_obs["v_c"], data.params_obs["SFR"]))
    fig_int_raw.suptitle("{0:}, v_c = {1:.1f} km/h, SFR = {2:.1f}".format(data.params_obs["name"], data.params_obs["v_c"], data.params_obs["SFR"]))
    fig_int_conv.suptitle("{0:}, v_c = {1:.1f} km/h, SFR = {2:.1f}".format(data.params_obs["name"], data.params_obs["v_c"], data.params_obs["SFR"]))


    for beta_el in betas:
    
        params.update(beta = beta_el)
    
        if load_sol_from_file == False:
        
            from sol_modules import get_profiles
            profiles = get_profiles(params, resol=1000)
        
            profiles.to_file()
            profiles.plot(ax=pltc.axs_sol)
            
            if profiles.check_nans() == True:
                string_nans = "Integration error: presence of nans"
            else:
                string_nans = "Integration successful"
                
            print("beta=", beta_el, "\t", string_nans)
        
        else:
        
            profiles = load_from_file(params, type_class = "sol_profiles")
         
            if profiles.check_nans() == True:
                string_nans = "Integration error: presence of nans"
            else:
                string_nans = "Integration successful"
                
            print("beta=", beta_el, "\t", string_nans)
      
            ionization_state = get_ionization_states(profiles, params)
        
            sigma_CII = get_surface_density(profiles, ionization_state, params, rmax=16)
        
            intensity_raw = get_intensity_raw(sigma_CII, params, data.params_obs)
        
            intensity_conv = get_intensity_convolved(intensity_raw, params, data.params_obs, data)
                     
            #profiles.to_file()
            profiles.plot(ax=axs_sol, label=r"$\beta$={:.1f}".format(beta_el))
            
            
            #ionization_state.to_file()
            ionization_state.plot(ax=axs_ion,  label=r"$\beta$={:.1f}".format(beta_el))
            
            sigma_CII.to_file()
            #sigma_CII.plot(ax=ax_sigma)            

            #intensity_raw.to_file()
            intensity_raw.plot(ax=ax_int_raw,  label=r"$\beta$={:.1f}".format(beta_el))
            
            #intensity_conv.to_file()
            intensity_conv.plot(ax=ax_int_conv,  label=r"$\beta$={:.1f}".format(beta_el))
            
            
    ax_int_conv.errorbar(data.x/(1000*nc.pc), data.data, yerr=data.err, \
                         markerfacecolor='C3',markeredgecolor='C3', marker='o',\
                         linestyle='', ecolor = 'C3')
        
    
    fig_sol.legend(loc="lower center", ncol=8, fontsize="small")
    fig_ion.legend(loc="lower center", ncol=8, fontsize="small")
    fig_int_raw.legend(loc="lower center", ncol=8, fontsize="small")
    fig_int_conv.legend(loc="lower center", ncol=8, fontsize="small")

    time_elapsed = (time.perf_counter() - time_start)

    print("total time elapsed (s)=", time_elapsed)


