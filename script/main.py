"""
Created on Sun Mar  7 17:16:03 2021

@author: anna
"""



import os
import numpy as np
import mydir
import natconst as nc

from post_sol_modules import get_ionization_states, get_surface_density, get_intensity_raw, get_intensity_convolved, get_chi2

from model_classes import load_from_file
from load_data import obs_data_list, names, names_CII_halo, names_wo_CII_halo
import plot_config as pltc


import time




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

load_sol_from_file = False

to_file = True

plot_hydro = False

plot_emission = False


betas = np.asarray([1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6])
#betas = [3.0,3.2,3.4]


datas = obs_data_list

data_counter = 0

for data in datas:
    
    if data.params_obs["name"] not in names_CII_halo:
        pass
    else:
        print("#####################")
        print(data.params_obs["name"], "(number {})".format(data_counter) )
        print("#####################")
    
        params = dict([("beta", 1.0), 
                   ("SFR", data.params_obs["SFR"]),
                   ("f_esc", 0.), 
                   ("v_c", data.params_obs["v_c"]),
                   ("redshift", data.params_obs["redshift"]),
                   ("Zeta", 1.0),
                   ("alfa", 1.0),
                   ("R_in", 0.3)])    
        
        time_start = time.perf_counter()
    
        if plot_hydro:
            fig_sol, axs_sol = pltc.plot_configurator(plot_type="sol")    
            fig_ion, axs_ion = pltc.plot_configurator(plot_type="ion")
    
            fig_sol.suptitle("{0:}, v_c = {1:.1f} km/h, SFR = {2:.1f}".format(data.params_obs["name"], data.params_obs["v_c"], data.params_obs["SFR"]))
            fig_ion.suptitle("{0:}, v_c = {1:.1f} km/h, SFR = {2:.1f}".format(data.params_obs["name"], data.params_obs["v_c"], data.params_obs["SFR"]))
        
        if plot_emission:
            fig_int_raw, ax_int_raw = pltc.plot_configurator(plot_type="int") 
            fig_int_conv, ax_int_conv = pltc.plot_configurator(plot_type="int")
        
            ax_int_raw.set_ylim((1e-3,1e2))
            ax_int_conv.set_ylim((1e-3,1e2))
 
            fig_int_raw.suptitle("{0:}, v_c = {1:.1f} km/h, SFR = {2:.1f}".format(data.params_obs["name"], data.params_obs["v_c"], data.params_obs["SFR"]))
            fig_int_conv.suptitle("{0:}, v_c = {1:.1f} km/h, SFR = {2:.1f}".format(data.params_obs["name"], data.params_obs["v_c"], data.params_obs["SFR"]))
    
        beta_counter = 0
        
        for beta_el in betas:
        
            params.update(beta = beta_el)
        
            if load_sol_from_file == False:
            
                from sol_modules import get_profiles
                profiles = get_profiles(params, resol=1000)
            
                if to_file:
                    profiles.to_file()

                if plot_hydro:
                    profiles.plot(ax=axs_sol)
                
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
                              
                ionization_state = get_ionization_states(profiles, params)
            
                sigma_CII = get_surface_density(profiles, ionization_state, params, rmax=18, h_resol=100, add_CMB_suppression=True)
            
                intensity_raw = get_intensity_raw(sigma_CII, params, data.params_obs)
            
                intensity_conv = get_intensity_convolved(intensity_raw, params, data.params_obs, data, add_central_contribution=False)

                if to_file:                         
                    profiles.to_file()
                    ionization_state.to_file()
                    sigma_CII.to_file()
                    intensity_raw.to_file()
                    intensity_conv.to_file()

                
                if plot_hydro:
                    profiles.plot(ax=axs_sol, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))
                    ionization_state.plot(ax=axs_ion,  label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))

                if plot_emission:                
                    #sigma_CII.plot(ax=ax_sigma)            
    
                    intensity_raw.plot(ax=ax_int_raw,  label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))
                
                    intensity_conv.plot(ax=ax_int_conv,  label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))
                
                    sigma_CII_no_CMB = get_surface_density(profiles, ionization_state, params, rmax=18, h_resol=100, add_CMB_suppression=False)
                    intensity_raw_no_CMB = get_intensity_raw(sigma_CII_no_CMB, params, data.params_obs)
                    intensity_conv_no_CMB = get_intensity_convolved(intensity_raw_no_CMB, params, data.params_obs, data, add_central_contribution=False)
                        
                    intensity_conv_no_CMB.plot(ax=ax_int_conv, color="C{}".format(beta_counter), linestyle='--')

            beta_counter+=1
                
            if plot_emission:
                # data   

                ax_int_conv.errorbar(data.x/(1000*nc.pc), data.data, yerr=data.err, \
                                 markerfacecolor='C3',markeredgecolor='C3', marker='o',\
                                 linestyle='', ecolor = 'C3')
            
                # central contribution
                
                luminosity_central = nc.ls * 1e7 * data.params_obs["SFR"]
        
                factor = luminosity_central / np.trapz(2*np.pi*data.x_beam*data.beam, data.x_beam)
                
                ax_int_conv.plot(data.x_beam/1e3/nc.pc, data.beam*factor, linestyle="--", color="gray")
                
        
            if plot_hydro:
                fig_sol.legend(loc="lower center", ncol=8, fontsize="small")
                fig_ion.legend(loc="lower center", ncol=8, fontsize="small")
                
            if plot_emission:
                fig_int_raw.legend(loc="lower center", ncol=8, fontsize="small")
                fig_int_conv.legend(loc="lower center", ncol=8, fontsize="small")
    
        
            chi2 = get_chi2(intensity_conv, data)
        
            print("beta=", beta_el, "\t", string_nans, "\t chi2=", chi2)
   
        time_elapsed = (time.perf_counter() - time_start)
    
        print("total time elapsed (s)=", time_elapsed)

        data_counter+=1

