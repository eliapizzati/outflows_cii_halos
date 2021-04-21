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
from load_data import obs_data_list, names, names_CII_halo, names_wo_CII_halo, names_other,  observational_data_fuji
import plot_config as pltc


import time



load_sol_from_file = True

to_file = False

plot_hydro = False

plot_emission = True

save_chi2 = False

plot_eta = False

f_esc_ion = 0.0

f_esc_FUV = 0.0

betas = np.asarray([1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,\
                    3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0])
#betas = np.asarray([1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1])
#betas = np.asarray([1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8])
betas = np.asarray([4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9])


datas = obs_data_list

data_counter = 0

chi2_names = []

datas_real = []

data_container_name = "other"

for data in datas:
    
    if data.params_obs["name"] not in names_CII_halo: #or data.params_obs["name"] != "DEIMOS_COSMOS_881725":
    #if data.params_obs["name"] in names_wo_CII_halo or data.params_obs["name"] in names_CII_halo:#names_wo_CII_halodata.params_obs["name"] != "DEIMOS_COSMOS_881725":
    #if data.params_obs["name"] != "vuds_cosmos_5110377875":
        pass
    else:
        datas_real.append(data)
        print("#####################")
        print(data.params_obs["name"], "(number {})".format(data_counter) )
        print("#####################")
    
        params = dict([("DM_model", "NFW"),
                   ("beta", 1.0), 
                   ("SFR", data.params_obs["SFR"]),
                   ("f_esc_ion", f_esc_ion), 
                   ("f_esc_FUV", f_esc_FUV), 
                   ("v_c", data.params_obs["v_c"]),
                   ("M_vir", data.params_obs["M_vir"]),
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
            fig_int_raw, ax_int_raw = pltc.plot_configurator(plot_type="int", xlim=15) 
            fig_int_conv, ax_int_conv = pltc.plot_configurator(plot_type="int", xlim=15)
        
            ax_int_raw.set_ylim((1e-3,1e2))
            ax_int_conv.set_ylim((1e-3,1e2))
 
            fig_int_raw.suptitle("{0:}, v_c = {1:.1f} km/h, SFR = {2:.1f}".format(data.params_obs["name"], data.params_obs["v_c"], data.params_obs["SFR"]))
            fig_int_conv.suptitle("{0:}, v_c = {1:.1f} km/h, SFR = {2:.1f}".format(data.params_obs["name"], data.params_obs["v_c"], data.params_obs["SFR"]))
    
        if plot_eta:
    
            fig_eta, ax_eta = pltc.plot_configurator(plot_type="eta", xlim=15) 

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
                
                    sigma_CII_no_CMB = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500, add_CMB_suppression=False)
                    intensity_raw_no_CMB = get_intensity_raw(sigma_CII_no_CMB, params, data.params_obs)
                    intensity_conv_no_CMB = get_intensity_convolved(intensity_raw_no_CMB, params, data.params_obs, data, add_central_contribution=False)
                        
                    #intensity_conv_no_CMB.plot(ax=ax_int_conv, color="C{}".format(beta_counter), linestyle='--')

                if plot_eta: 
                    ax_eta.plot(intensity_conv.h/1e3/nc.pc, intensity_conv.eta, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))

            beta_counter+=1
                
            
            if load_sol_from_file == True:
                
                if plot_emission:
                        # data   
    
                        ax_int_conv.errorbar(data.x/(1000*nc.pc), data.data, yerr=data.err, \
                                     markerfacecolor='maroon',markeredgecolor='maroon', marker='o',\
                                     linestyle='', ecolor = 'maroon')
                        
                        ax_int_conv.errorbar(observational_data_fuji.x/(1000*nc.pc), observational_data_fuji.data, yerr=observational_data_fuji.err, \
                                     markerfacecolor='navy',markeredgecolor='navy', marker='d',\
                                     linestyle='', ecolor = 'navy')
                        
                        # central contribution
                    
                        luminosity_central = nc.ls * 1e7 * data.params_obs["SFR"]
            

                        factor_lum = luminosity_central / np.trapz(2*np.pi*data.x_beam*data.beam, data.x_beam)
                        factor_data = data.data[0]/data.beam[0]
                                            
                        ax_int_conv.plot(data.x_beam/1e3/nc.pc, data.beam*factor_data, linestyle="--", color="gray")
                        

            
                if plot_hydro:
                        fig_sol.legend(loc="lower center", ncol=8, fontsize="small")
                        fig_ion.legend(loc="lower center", ncol=8, fontsize="small")
                        
                if plot_emission:
                        fig_int_raw.legend(loc="lower center", ncol=8, fontsize="small")
                        fig_int_conv.legend(loc="lower center", ncol=8, fontsize="small")
        
            
                if plot_eta:
                    fig_eta.legend(loc="lower center", fontsize="small", ncol=8)

                chi2 = get_chi2(intensity_conv, data)
                
                likelihood = np.exp(-chi2/data.data.shape[0])  

                chi2_betas.append(likelihood)
                
                    
                print("beta=", beta_el, "\t", string_nans, "\t chi2/ndof=", "{:.1f}/{:.0f}".format(chi2, len(data.data)))
        
        if load_sol_from_file == True:

            chi2_betas = np.asarray(chi2_betas)
            
            beta_best_fit = np.nan
            
            if chi2_betas.max() >= 0.1:
                beta_best_fit = betas[chi2_betas == chi2_betas.max()][0]
            
            print(chi2_betas)
            print("beta best fit = ", beta_best_fit)
                        
            data.params_obs.update(beta_best_fit=beta_best_fit)
            
            chi2_names.append(chi2_betas)
            
            
        time_elapsed = (time.perf_counter() - time_start)
        
        print("total time elapsed (s)=", time_elapsed)

        data_counter+=1
       
if save_chi2 == True:     
    
    chi2_names = np.asarray(chi2_names)
    
    out_filename = os.path.join(mydir.data_dir, "data_chi2", data_container_name)
    
    np.save(out_filename, chi2_names)


    