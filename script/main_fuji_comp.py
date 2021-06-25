# -*- coding: utf-8 -*-
"""
Created on Sat Apr 24 09:59:08 2021

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



alpha = 0.3

fesc_comp = False

NFW_comp = True

f_esc_ion = 0.

f_esc_FUV = 0.

betas = np.asarray([1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,\
                    3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0])
#betas = np.asarray([1.8,2.1,2.4,2.7,3.0,3.3,3.6,3.9,4.2,4.5])
#betas = [3.0]
#betas = np.asarray([2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0])
#betas = np.asarray([3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0])
betas = np.asarray([1.3,1.6,1.9,2.2,2.5,2.8,3.1,3.4,3.7,4.0])

data = observational_data_fuji
  

time_start = time.perf_counter()

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

fig_sol, axs_sol = pltc.plot_configurator(plot_type="sol")    
fig_ion, axs_ion = pltc.plot_configurator(plot_type="ion")

fig_sol.suptitle("v_c = {:.1f} km/h, SFR = {:.1f}, f_escFUV = 0.0, f_escION = 0.0; NFW model"\
                     .format(params["v_c"], params["SFR"]))
fig_ion.suptitle("v_c = {:.1f} km/h, SFR = {:.1f}, f_escFUV = 0.0, f_escION = 0.0; NFW model"\
                     .format(params["v_c"], params["SFR"]))

fig_int_raw, ax_int_raw = pltc.plot_configurator(plot_type="int", xlim=10) 
fig_int_conv, ax_int_conv = pltc.plot_configurator(plot_type="int", xlim=10)

ax_int_raw.set_ylim((1e-3,1e2))
ax_int_conv.set_ylim((1e-3,1e2))
 
fig_int_raw.suptitle("v_c = {:.1f} km/h, SFR = {:.1f}, f_escFUV = 0.0, f_escION = 0.0; NFW model"\
                     .format(params["v_c"], params["SFR"]))
fig_int_conv.suptitle("v_c = {:.1f} km/h, SFR = {:.1f}, f_escFUV = 0.0, f_escION = 0.0; NFW model"\
                     .format(params["v_c"], params["SFR"]))

fig_eta, ax_eta = pltc.plot_configurator(plot_type="eta", xlim=10) 
fig_eta.suptitle("v_c = {:.1f} km/h, SFR = {:.1f}, f_escFUV = 0.0, f_escION = 0.0; NFW model"\
                     .format(params["v_c"], params["SFR"]))
    
    
beta_counter = 0

chi2_betas = []

for beta_el in betas: 
    
    params.update(beta = beta_el)
    
    profiles = load_from_file(params, class_type = "profiles")                  
    ionization_state = get_ionization_states(profiles, params)
    sigma_CII = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500, add_CMB_suppression=True)    
    intensity_raw = get_intensity_raw(sigma_CII, params, data.params_obs)
    intensity_conv = get_intensity_convolved(intensity_raw, params, data.params_obs, data, add_central_contribution=False)

    profiles.plot(ax=axs_sol, color="C{}".format(beta_counter), alpha=alpha)
    ionization_state.plot(ax=axs_ion, color="C{}".format(beta_counter), alpha=alpha)

    intensity_raw.plot(ax=ax_int_raw, color="C{}".format(beta_counter), alpha=alpha)
    intensity_conv.plot(ax=ax_int_conv, color="C{}".format(beta_counter),alpha=alpha)
        
    ax_eta.plot(intensity_conv.h/1e3/nc.pc, intensity_conv.eta, color="C{}".format(beta_counter),alpha=alpha)
    
    if fesc_comp:
        params.update(f_esc_FUV = 0.2)
        #params.update(f_esc_ion = 0.2)
    
        profiles_FUV = load_from_file(params, class_type = "profiles")                  
        ionization_state_FUV = get_ionization_states(profiles_FUV, params)
        sigma_CII_FUV = get_surface_density(profiles_FUV, ionization_state_FUV, params, rmax=30, h_resol=500, add_CMB_suppression=True)    
        intensity_raw_FUV = get_intensity_raw(sigma_CII_FUV, params, data.params_obs)
        intensity_conv_FUV = get_intensity_convolved(intensity_raw_FUV, params, data.params_obs, data, add_central_contribution=False)
    
        profiles_FUV.plot(ax=axs_sol, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))
        ionization_state_FUV.plot(ax=axs_ion, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))
    
        intensity_raw_FUV.plot(ax=ax_int_raw, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))
        intensity_conv_FUV.plot(ax=ax_int_conv, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))
            
        ax_eta.plot(intensity_conv_FUV.h/1e3/nc.pc, intensity_conv_FUV.eta, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))
      
        params.update(f_esc_FUV = f_esc_FUV)
        params.update(f_esc_ion = f_esc_ion)
        
    if NFW_comp:
        params.update(DM_model = "NFW")
        
        profiles_NFW = load_from_file(params, class_type = "profiles")                  
        ionization_state_NFW = get_ionization_states(profiles_NFW, params)
        sigma_CII_NFW = get_surface_density(profiles_NFW, ionization_state_NFW, params, rmax=30, h_resol=500, add_CMB_suppression=True)    
        intensity_raw_NFW = get_intensity_raw(sigma_CII_NFW, params, data.params_obs)
        intensity_conv_NFW = get_intensity_convolved(intensity_raw_NFW, params, data.params_obs, data, add_central_contribution=False)
    
        profiles_NFW.plot(ax=axs_sol, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))
        ionization_state_NFW.plot(ax=axs_ion, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))
    
        intensity_raw_NFW.plot(ax=ax_int_raw, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))
        intensity_conv_NFW.plot(ax=ax_int_conv, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))
            
        ax_eta.plot(intensity_conv_NFW.h/1e3/nc.pc, intensity_conv_NFW.eta, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))
      
        params.update(DM_model = None)
        

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

    fig_sol.legend(loc="lower center", fontsize="small", ncol=ncol)
    fig_ion.legend(loc="lower center", fontsize="small", ncol=ncol)
    
    fig_int_raw.legend(loc="lower center", fontsize="small", ncol=ncol)
    fig_int_conv.legend(loc="lower center", fontsize="small", ncol=ncol)

    fig_eta.legend(loc="lower center", fontsize="small", ncol=ncol)


    chi2 = get_chi2(intensity_conv, data)
    
    chi2_betas.append(chi2)
        
    print("beta=", beta_el, "\t chi2/ndof=", "{:.1f}/{:.0f}".format(chi2, len(data.data)))
   
    beta_counter+=1

time_elapsed = (time.perf_counter() - time_start)

print("total time elapsed (s)=", time_elapsed)

