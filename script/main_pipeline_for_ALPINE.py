"""
This script runs the entire pipeline for (selected) ALPINE systems and show some plots/saves data in files depending on
the switches
"""



import os
import numpy as np
import my_dir
import natconst as nc

import matplotlib.pyplot as plt
from post_sol_modules import get_ionization_states, get_surface_density, get_intensity_raw, \
                      get_intensity_convolved, get_chi2

from model_classes import load_from_file
from load_data import obs_data_list, observational_data_fuji
import plot_config as pltc
from my_utils import get_vc_from_virial_mass


import time



load_sol_from_file = True

to_file = False

plot_hydro = True

plot_emission = False

save_chi2 = False

plot_eta = False

plot_vc_uncertainty = False

plot_SFR_uncertainty = False

f_esc_ion = 0.0

f_esc_FUV = 0.0

betas = np.asarray([4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,\
                    5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,\
                    6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,\
                    7.0,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9])
#betas = np.asarray([1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1])
#betas = np.asarray([1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8])
betas = np.asarray([6.80])
#betas = np.asarray([7.9,8.0,8.1,8.2,8.3,8.4,8.5,8.6,8.7])
#betas = np.asarray([1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1])
#betas = np.asarray([4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9])
#betas = np.asarray([2.6,2.7,2.8,2.9,3.0])


datas = obs_data_list

data_counter = 0

chi2_names = []

datas_real = []

data_container_name = "CII_halo_NFW"

for data in datas:
    
    #if data.params_obs["name"] not in names_CII_halo:
    #if data.params_obs["name"] in names_wo_CII_halo or data.params_obs["name"] in names_CII_halo:#names_wo_CII_halodata.params_obs["name"] != "DEIMOS_COSMOS_881725":
    if data.params_obs["name_short"] != "VC_5110377875":
        pass
    else:
        datas_real.append(data)
        print("#####################")
        print(data.params_obs["name"], "(number {})".format(data_counter) )
        print("#####################")
        params = dict([("DM_model", "NFW+disk"),
                   ("beta", 1.0),
                   ("SFR", data.params_obs["SFR"]),
                   ("f_esc_ion", f_esc_ion),
                   ("f_esc_FUV", f_esc_FUV),
                   ("v_c", data.params_obs["v_c"]),
                   #("M_vir", data.params_obs["M_vir"]),
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
            #fig_int_raw, ax_int_raw = pltc.plot_configurator(plot_type="int", xlim=17)
            fig_int_conv, ax_int_conv = pltc.plot_configurator(plot_type="int", xlim=17)
        
            #ax_int_raw.set_ylim((1e-3,1e2))
            ax_int_conv.set_ylim((1e-3,1e2))
 
            #fig_int_raw.suptitle("{0:}, v_c = {1:.1f} km/h, SFR = {2:.1f}".format(data.params_obs["name"], data.params_obs["v_c"], data.params_obs["SFR"]))
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

                #intensity_raw.plot(ax=ax_int_raw,  label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))

                intensity_conv.plot(ax=ax_int_conv,  label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))

                sigma_CII_no_CMB = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500, add_CMB_suppression=False)
                intensity_raw_no_CMB = get_intensity_raw(sigma_CII_no_CMB, params, data.params_obs)
                intensity_conv_no_CMB = get_intensity_convolved(intensity_raw_no_CMB, params, data.params_obs, data, add_central_contribution=False)

                #intensity_conv_no_CMB.plot(ax=ax_int_conv, color="C{}".format(beta_counter), linestyle='--')

            if plot_eta:
                ax_eta.plot(intensity_conv.h/1e3/nc.pc, intensity_conv.eta, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter))

            if plot_vc_uncertainty:

                M_vir_up = data.params_obs["M_vir"]+data.params_obs["M_vir_err_up"]
                v_c_up = get_vc_from_virial_mass(M_vir_up, params["redshift"])/1e5

                M_vir_down = data.params_obs["M_vir"]-data.params_obs["M_vir_err_down"]
                v_c_down = get_vc_from_virial_mass(M_vir_down, params["redshift"])/1e5

                params.update(M_vir = M_vir_up)
                params.update(v_c = v_c_up)
                print("{:.1f}".format(v_c_up))

                if load_sol_from_file == False:
                    profiles_up = get_profiles(params, resol=1000)
                else:
                    profiles_up = load_from_file(params, class_type = "profiles")
                ionization_state_up = get_ionization_states(profiles_up, params)
                sigma_CII_up = get_surface_density(profiles_up, ionization_state_up, params, rmax=30, h_resol=500, add_CMB_suppression=True)
                intensity_raw_up = get_intensity_raw(sigma_CII_up, params, data.params_obs)
                intensity_conv_up = get_intensity_convolved(intensity_raw_up, params, data.params_obs, data, add_central_contribution=False)

                params.update(M_vir = M_vir_down)
                params.update(v_c = v_c_down)
                print("{:.1f}".format(v_c_down))

                if load_sol_from_file == False:
                    profiles_down = get_profiles(params, resol=1000)
                else:
                    profiles_down = load_from_file(params, class_type = "profiles")

                ionization_state_down = get_ionization_states(profiles_down, params)
                sigma_CII_down = get_surface_density(profiles_down, ionization_state_down, params, rmax=30, h_resol=500, add_CMB_suppression=True)
                intensity_raw_down = get_intensity_raw(sigma_CII_down, params, data.params_obs)
                intensity_conv_down = get_intensity_convolved(intensity_raw_down, params, data.params_obs, data, add_central_contribution=False)

                params.update(M_vir = data.params_obs["M_vir"])
                params.update(v_c = data.params_obs["v_c"])

                ax_int_conv.plot(intensity_conv_up.h/(1000*nc.pc), intensity_conv_up.var, color="C{}".format(beta_counter))
                ax_int_conv.plot(intensity_conv_down.h/(1000*nc.pc), intensity_conv_down.var, color="C{}".format(beta_counter))
                ax_int_conv.fill_between(intensity_conv.h/(1000*nc.pc), intensity_conv_down.var, intensity_conv_up.var, color="C{}".format(beta_counter), alpha=0.2)

                #ax_int_raw.plot(intensity_raw_up.h/(1000*nc.pc), intensity_raw_up.var, color="C{}".format(beta_counter))
                #ax_int_raw.plot(intensity_raw_down.h/(1000*nc.pc), intensity_raw_down.var, color="C{}".format(beta_counter))
                #ax_int_raw.fill_between(intensity_raw.h/(1000*nc.pc), intensity_raw_down.var, intensity_raw_up.var, color="C{}".format(beta_counter), alpha=0.2)

                if plot_hydro:
                    profiles_up.plot(ax=axs_sol,color="C{}".format(beta_counter))
                    profiles_down.plot(ax=axs_sol, color="C{}".format(beta_counter))

                    profiles_up_extended_v = np.interp(profiles_down.r, profiles_up.r, profiles_up.v) #right=0.)
                    profiles_up_extended_n = np.interp(profiles_down.r, profiles_up.r, profiles_up.n)# right=0.)
                    profiles_up_extended_T = np.interp(profiles_down.r, profiles_up.r, profiles_up.T)#right=0.)

                    axs_sol[0].fill_between(np.log10(profiles_down.r/(1000*nc.pc)),profiles_down.v/10**8,profiles_up_extended_v/10**8,\
                           color="C{}".format(beta_counter), alpha=0.2)
                    axs_sol[1].fill_between(np.log10(profiles_down.r/(1000*nc.pc)),np.log10(profiles_down.n), np.log10(profiles_up_extended_n),\
                           color="C{}".format(beta_counter), alpha=0.2)
                    axs_sol[2].fill_between(np.log10(profiles_down.r/(1000*nc.pc)),np.log10(profiles_down.T),np.log10(profiles_up_extended_T),\
                           color="C{}".format(beta_counter), alpha=0.2)

                    ionization_state_up.plot(ax=axs_ion,color="C{}".format(beta_counter))
                    ionization_state_down.plot(ax=axs_ion, color="C{}".format(beta_counter))

                    ion_up_extended_xe = np.interp(ionization_state_down.r, ionization_state_up.r, ionization_state_up.x_e) #right=0.)
                    ion_up_extended_xCII = np.interp(ionization_state_down.r, ionization_state_up.r, ionization_state_up.x_CII)# right=0.)

                    axs_ion[0].fill_between(np.log10(ionization_state_down.r/(1000*nc.pc)),1.-ionization_state_down.x_e,1.-ion_up_extended_xe,\
                           color="C{}".format(beta_counter), alpha=0.2)
                    axs_ion[1].fill_between(np.log10(ionization_state_down.r/(1000*nc.pc)),ionization_state_down.x_CII,ion_up_extended_xCII,\
                           color="C{}".format(beta_counter), alpha=0.2)


            if plot_SFR_uncertainty:

                SFR_up = data.params_obs["SFR"]+data.params_obs["SFR_err_up"]

                SFR_down = data.params_obs["SFR"]-data.params_obs["SFR_err_down"]

                params.update(SFR = SFR_up)

                if load_sol_from_file == False:
                    profiles_up = get_profiles(params, resol=1000)
                else:
                    profiles_up = load_from_file(params, class_type = "profiles")

                ionization_state_up = get_ionization_states(profiles_up, params)
                sigma_CII_up = get_surface_density(profiles_up, ionization_state_up, params, rmax=30, h_resol=500, add_CMB_suppression=True)
                intensity_raw_up = get_intensity_raw(sigma_CII_up, params, data.params_obs)
                intensity_conv_up = get_intensity_convolved(intensity_raw_up, params, data.params_obs, data, add_central_contribution=False)

                params.update(SFR = SFR_down)

                if load_sol_from_file == False:
                    profiles_down = get_profiles(params, resol=1000)
                else:
                    profiles_down = load_from_file(params, class_type = "profiles")

                ionization_state_down = get_ionization_states(profiles_down, params)
                sigma_CII_down = get_surface_density(profiles_down, ionization_state_down, params, rmax=30, h_resol=500, add_CMB_suppression=True)
                intensity_raw_down = get_intensity_raw(sigma_CII_down, params, data.params_obs)
                intensity_conv_down = get_intensity_convolved(intensity_raw_down, params, data.params_obs, data, add_central_contribution=False)

                params.update(SFR = data.params_obs["SFR"])

                ax_int_conv.plot(intensity_conv_up.h/(1000*nc.pc), intensity_conv_up.var, color="C{}".format(beta_counter))
                ax_int_conv.plot(intensity_conv_down.h/(1000*nc.pc), intensity_conv_down.var, color="C{}".format(beta_counter))
                ax_int_conv.fill_between(intensity_conv.h/(1000*nc.pc), intensity_conv_down.var, intensity_conv_up.var, color="C{}".format(beta_counter), alpha=0.2)

                #ax_int_raw.plot(intensity_raw_up.h/(1000*nc.pc), intensity_raw_up.var, color="C{}".format(beta_counter))
                #ax_int_raw.plot(intensity_raw_down.h/(1000*nc.pc), intensity_raw_down.var, color="C{}".format(beta_counter))
                #ax_int_raw.fill_between(intensity_raw.h/(1000*nc.pc), intensity_raw_down.var, intensity_raw_up.var, color="C{}".format(beta_counter), alpha=0.2)


                if plot_hydro:
                    profiles_up.plot(ax=axs_sol,color="C{}".format(beta_counter))
                    profiles_down.plot(ax=axs_sol, color="C{}".format(beta_counter))

                    profiles_up_extended_v = np.interp(profiles_down.r, profiles_up.r, profiles_up.v) #right=0.)
                    profiles_up_extended_n = np.interp(profiles_down.r, profiles_up.r, profiles_up.n)# right=0.)
                    profiles_up_extended_T = np.interp(profiles_down.r, profiles_up.r, profiles_up.T)#right=0.)

                    axs_sol[0].fill_between(np.log10(profiles_down.r/(1000*nc.pc)),profiles_down.v/10**8,profiles_up_extended_v/10**8,\
                           color="C{}".format(beta_counter), alpha=0.2)
                    axs_sol[1].fill_between(np.log10(profiles_down.r/(1000*nc.pc)),np.log10(profiles_down.n), np.log10(profiles_up_extended_n),\
                           color="C{}".format(beta_counter), alpha=0.2)
                    axs_sol[2].fill_between(np.log10(profiles_down.r/(1000*nc.pc)),np.log10(profiles_down.T),np.log10(profiles_up_extended_T),\
                           color="C{}".format(beta_counter), alpha=0.2)

                    ionization_state_up.plot(ax=axs_ion,color="C{}".format(beta_counter))
                    ionization_state_down.plot(ax=axs_ion, color="C{}".format(beta_counter))

                    ion_up_extended_xe = np.interp(ionization_state_down.r, ionization_state_up.r, ionization_state_up.x_e) #right=0.)
                    ion_up_extended_xCII = np.interp(ionization_state_down.r, ionization_state_up.r, ionization_state_up.x_CII)# right=0.)

                    axs_ion[0].fill_between(np.log10(ionization_state_down.r/(1000*nc.pc)),1.-ionization_state_down.x_e,1.-ion_up_extended_xe,\
                           color="C{}".format(beta_counter), alpha=0.2)
                    axs_ion[1].fill_between(np.log10(ionization_state_down.r/(1000*nc.pc)),ionization_state_down.x_CII,ion_up_extended_xCII,\
                           color="C{}".format(beta_counter), alpha=0.2)


            beta_counter+=1
                
            
            if True == True:
                
                if plot_emission:
                        # data   

                        alpine = ax_int_conv.errorbar(data.x/(1000*nc.pc), data.data, yerr=data.err, \
                                     markerfacecolor='maroon',markeredgecolor='maroon', marker='o',\
                                     linestyle='', ecolor = 'maroon')

                        fuji = ax_int_conv.errorbar(observational_data_fuji.x/(1000*nc.pc), observational_data_fuji.data, yerr=observational_data_fuji.err, \
                                     markerfacecolor='navy',markeredgecolor='navy', marker='d',\
                                     linestyle='', ecolor = 'navy')

                        ax_int_conv.legend([alpine, fuji], [data.params_obs["name_short"], "Fujimoto+19"])##
                        #ax_int_conv.legend([alpine], [data.params_obs["name_short"]])  ##

                        # central contribution

                        luminosity_central = nc.ls * 1e7 * data.params_obs["SFR"]
            

                        factor_lum = luminosity_central / np.trapz(2*np.pi*data.x_beam*data.beam, data.x_beam)
                        factor_data = data.data[0]/data.beam[0]
                                            
                        ax_int_conv.plot(data.x_beam/1e3/nc.pc, data.beam*factor_data, linestyle="--", color="gray")
                        

            
                if plot_hydro:
                        fig_sol.legend(loc="lower center", ncol=8, fontsize="small")
                        fig_ion.legend(loc="lower center", ncol=8, fontsize="small")
                        
                if plot_emission:
                        #fig_int_raw.legend(loc="lower center", ncol=8, fontsize="small")
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
    
    out_filename = os.path.join(my_dir.data_dir, "data_chi2", data_container_name)
    
    np.save(out_filename, chi2_names)

plt.show()
