# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 11:37:50 2021

@author: anna
"""

import os
import numpy as np
import mydir
import natconst as nc

from post_sol_modules import get_ionization_states, get_surface_density, get_intensity_raw, \
                             get_intensity_convolved, get_chi2

from model_classes import load_from_file

from load_data import obs_data_list, names, names_CII_halo, names_wo_CII_halo, observational_data_fuji

import plot_config as pltc

import time



load_sol_from_file = True

to_file = False

plot_hydro = True

plot_emission = True

plot_eta = True

f_esc_ion = 0.

f_esc_FUV = 0.

betas = np.asarray([1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,\
                    3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0])
#betas = np.asarray([1.8,2.1,2.4,2.7,3.0,3.3,3.6,3.9,4.2,4.5])
#betas = [3.0]
#betas = np.asarray([2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0])
betas = np.asarray([1.3,1.6,1.9,2.2,2.5,2.8,3.1,3.4,3.7,4.0])
betas = np.asarray([1.6,4.0])

data = observational_data_fuji

params = dict([("DM_model", None),
           ("beta", 1.0), 
           ("SFR", 50.),
           ("f_esc_ion", f_esc_ion), 
           ("f_esc_FUV", f_esc_FUV), 
           ("v_c", 175.),
           ("redshift", data.params_obs["redshift"]),
           ("Zeta", 1.0),
           ("alfa", 1.0),
           ("R_in", 0.3)])    

time_start = time.perf_counter()

if plot_hydro:
    fig_sol, axs_sol = pltc.plot_configurator(plot_type="sol")    
    fig_ion, axs_ion = pltc.plot_configurator(plot_type="ion")

    fig_sol.suptitle("v_c = {:.1f} km/s, SFR = {:.1f}, f_escFUV = {:.1f}, f_escION = {:.1f}"\
                     .format(params["v_c"], params["SFR"], params["f_esc_FUV"], params["f_esc_ion"]))
    fig_ion.suptitle("v_c = {:.1f} km/s, SFR = {:.1f}, f_escFUV = {:.1f}, f_escION = {:.1f}"\
                     .format(params["v_c"], params["SFR"], params["f_esc_FUV"], params["f_esc_ion"]))

if plot_emission:
    fig_int_raw, ax_int_raw = pltc.plot_configurator(plot_type="int", xlim=10) 
    fig_int_conv, ax_int_conv = pltc.plot_configurator(plot_type="int", xlim=10)

    ax_int_raw.set_ylim((1e-3,1e2))
    ax_int_conv.set_ylim((1e-3,1e2))
 
    fig_int_raw.suptitle("v_c = {:.1f} km/s, SFR = {:.1f}, f_escFUV = {:.1f}, f_escION = {:.1f}"\
                     .format(params["v_c"], params["SFR"], params["f_esc_FUV"], params["f_esc_ion"]))
    fig_int_conv.suptitle("v_c = {:.1f} km/s, SFR = {:.1f}, f_escFUV = {:.1f}, f_escION = {:.1f}"\
                     .format(params["v_c"], params["SFR"], params["f_esc_FUV"], params["f_esc_ion"]))

if plot_eta:
    
    fig_eta, ax_eta = pltc.plot_configurator(plot_type="eta", xlim=10) 
    fig_eta.suptitle("v_c = {:.1f} km/s, SFR = {:.1f}, f_escFUV = {:.1f}, f_escION = {:.1f}"\
                     .format(params["v_c"], params["SFR"], params["f_esc_FUV"], params["f_esc_ion"]))
    
    
beta_counter = 0

chi2_betas = []

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
    
        profiles = load_from_file(params, class_type = "profiles")
     
        if profiles.check_nans() == True:
            string_nans = "Integration error: presence of nans"
        else:
            string_nans = "Integration successful"
                      
        ionization_state = get_ionization_states(profiles, params)
    
        sigma_CII = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500, add_CMB_suppression=True)
        
        intensity_raw = get_intensity_raw(sigma_CII, params, data.params_obs)

        intensity_conv = get_intensity_convolved(intensity_raw, params, data.params_obs, data, add_central_contribution=False)

        print(np.trapz(1./profiles.v, profiles.r)/1e6/nc.year)
        if to_file:                         
            profiles.to_file()
            ionization_state.to_file()
            sigma_CII.to_file()
            intensity_raw.to_file()
            intensity_conv.to_file()

        
        if plot_hydro:
            profiles.plot(ax=axs_sol, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))
            ionization_state.plot(ax=axs_ion,  label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))

            try_radius=np.linspace(.3,30.,1000)
            axs_sol[0].plot(np.log10(try_radius), np.sqrt(395.1**2-2*(175*1.51)**2*(np.log10(try_radius)+0.484))/1e3, color="black")
            axs_sol[0].plot(np.log10(try_radius), np.sqrt(664**2-2*(175*1.51)**2*(np.log10(try_radius)+0.375))/1e3, color="black")


            axs_sol[1].plot(np.log10(try_radius), np.log10(10**1.38*(10**(-0.484)/try_radius)**2*395.1/np.sqrt(395.1**2-2*(175*1.51)**2*(np.log10(try_radius)+0.484))), color="black")
            axs_sol[1].plot(np.log10(try_radius), np.log10(10**0.553*(10**(-0.375)/try_radius)**2*664/np.sqrt(664**2-2*(175*1.51)**2*(np.log10(try_radius)+0.375))), color="black")
            
        if plot_emission:                
            #sigma_CII.plot(ax=ax_sigma)            

            intensity_raw.plot(ax=ax_int_raw,  label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))
        
            intensity_conv.plot(ax=ax_int_conv,  label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))
        
            sigma_CII_no_CMB = get_surface_density(profiles, ionization_state, params, rmax=18, h_resol=100, add_CMB_suppression=False)
            intensity_raw_no_CMB = get_intensity_raw(sigma_CII_no_CMB, params, data.params_obs)
            intensity_conv_no_CMB = get_intensity_convolved(intensity_raw_no_CMB, params, data.params_obs, data, add_central_contribution=False)
                
            #intensity_conv_no_CMB.plot(ax=ax_int_conv, color="C{}".format(beta_counter), linestyle='--')
            #intensity_raw_no_CMB.plot(ax=ax_int_raw, color="C{}".format(beta_counter), linestyle='--')
            
        if plot_eta: 
            ax_eta.plot(intensity_conv.h/1e3/nc.pc, intensity_conv.eta, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))

    beta_counter+=1
        
    
    if load_sol_from_file == True:
        
        if plot_emission:
                # data   

                fuji = ax_int_conv.errorbar(data.x/(1000*nc.pc), data.data, yerr=data.err, \
                             markerfacecolor='maroon',markeredgecolor='maroon', marker='o',\
                             linestyle='', ecolor = 'maroon')
                
                ax_int_conv.legend([fuji], ["Fujimoto+19"])##

                # central contribution
            
                luminosity_central = nc.ls * 1e7 * params["SFR"]
    
                factor_lum = luminosity_central / np.trapz(2*np.pi*data.x_beam*data.beam, data.x_beam)
                factor_data = data.data[0]/data.beam[0]
                
                ax_int_conv.plot(data.x_beam/1e3/nc.pc, data.beam*factor_data, linestyle="--", color="gray")
                
                
        ncol = 5
        if plot_hydro:
                fig_sol.legend(loc="lower center", fontsize="small", ncol=ncol)
                fig_ion.legend(loc="lower center", fontsize="small", ncol=ncol)
                
        if plot_emission:
                fig_int_raw.legend(loc="lower center", fontsize="small", ncol=ncol)
                fig_int_conv.legend(loc="lower center", fontsize="small", ncol=ncol)
        
        if plot_eta:
                fig_eta.legend(loc="lower center", fontsize="small", ncol=ncol)
    
    
        chi2 = get_chi2(intensity_conv, data)
        
        chi2_betas.append(chi2)
            
        print("beta=", beta_el, "\t", string_nans, "\t chi2/ndof=", "{:.1f}/{:.0f}".format(chi2, len(data.data)))
   
time_elapsed = (time.perf_counter() - time_start)

print("total time elapsed (s)=", time_elapsed)

