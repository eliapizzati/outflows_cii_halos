# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 16:40:06 2021

@author: anna
"""


import numpy as np
import matplotlib.pyplot as plt
import natconst as nc


import matplotlib



matplotlib.rcParams.update({
        "font.size": 16.0,
        "font.family": 'sans-serif',
#        "font.sans-serif": ['Helvetica'],
        "axes.titlesize": 16.,
        "axes.labelsize": 16.,
        "xtick.labelsize": 16.,
        "ytick.labelsize": 16.,
        "xtick.major.size": 6.0,
        "ytick.major.size": 6.0,
        "xtick.minor.size": 4.0,
        "ytick.minor.size": 4.0,
        "legend.fontsize": 13.0,
        "legend.frameon": False
 #       "figure.dpi": 200,
 #       "savefig.dpi": 200,
#        "text.usetex": True
})



plot1 = True  #single profiles results 
plot2 = False      # SFR and vc uncertainties


params = dict([("DM_model", "NFW"),
           ("beta", 1.0), 
           ("SFR", 50.),
           ("f_esc_ion", 0.), 
           ("f_esc_FUV", 0.), 
           ("v_c", 200.),
           ("redshift", 5.),
           ("Zeta", 1.0),
           ("alfa", 1.0),
           ("R_in", 0.3)])    


if plot1:
    """
    #single profiles results 
    """
    from sol_modules import get_profiles, get_ionization_states, get_surface_density, get_intensity_raw, get_intensity_convolved

    from load_data import obs_data_list

    betas = np.arange(1.0,9.0,0.8)

    cmap_rend_col = matplotlib.cm.get_cmap('tab10')
    
    norm = matplotlib.colors.Normalize(vmin=betas.min()-0.4, vmax=betas.max()+0.4)


    fig, axs = plt.subplots(4, 2, sharey=True, sharex=True, figsize=(1.3*8.27,1.3*12.))

    for ax in axs:
        ax.set_xlabel("b [kpc]")
        ax.set_ylabel(r"flux [mJy/arcsec$^2$]")
        ax.set_yscale('log')
        ax.set_ylim((0.01,12))    
        ax.set_xlim((0.3,10))
    


    for data,  data_counter in zip(obs_data_list, range(len(obs_data_list))):
    
        for i in range(len(betas)):
        
        
            params.update(beta = betas[i])
            
            profiles = get_profiles(params, resol=1000)
            ionization_state = get_ionization_states(profiles, params)
            sigma_CII = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500, add_CMB_suppression=True)    
            intensity_raw = get_intensity_raw(sigma_CII, params, data.params_obs)
            intensity_conv = get_intensity_convolved(intensity_raw, params, data.params_obs, data, add_central_contribution=False)
         
            
            axs[data_counter].plot(intensity_conv.h/(1e3*nc.pc),intensity_conv.var, color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))
            
            
        alpine = axs[data_counter].errorbar(data.x/(1000*nc.pc), data.data, yerr=data.err, \
                 markerfacecolor='maroon',markeredgecolor='maroon', marker='o',\
                 linestyle='', ecolor = 'maroon')


        factor_data = data.data[0]/data.beam[0]        
        axs[data_counter].plot(data.x_beam/1e3/nc.pc, data.beam*factor_data, linestyle="--", color="gray")

        axs[data_counter].legend([alpine], ["{}".format(data["name_short"])])##
    
    plt.subplots_adjust(left = 0.1,  # the left side of the subplots of the figure
    right = 0.95,   # the right side of the subplots of the figure
    bottom = 0.3,  # the bottom of the subplots of the figure
    top = 0.95,     # the top of the subplots of the figure
    wspace = 0.05,  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    hspace = 0.05)  # the amount of height reserved for space between subplots,
                  # expressed as a fraction of the average axis height
    
    cax = plt.axes([0.15, 0.13,  0.72,0.03])
    
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    #betas_ticks = betas[10*betas%2 ==0]
    betas_ticks = betas
    
    
    cb = fig.colorbar(cmap, ticks=betas_ticks.round(1), cax=cax, orientation='horizontal')
    cb.set_label(r'$\eta$', rotation=0.)

   

if plot2:
    """
    #SFR and vc uncertainties
    """
    from sol_modules import get_profiles, get_ionization_states, get_surface_density, get_intensity_raw, get_intensity_convolved

    from load_data import obs_data_list

    from my_utils import get_vc_from_virial_mass

    betas = np.arange(1.0,9.0,0.8)

    cmap_rend_col = matplotlib.cm.get_cmap('tab10')
    
    norm = matplotlib.colors.Normalize(vmin=betas.min()-0.05, vmax=betas.max()+0.05)

    fig, (ax_vc, ax_sfr) = plt.subplots(1, 2, sharey=True, figsize=(1.3*8.27,1.3*4.))


    data = obs_data_list[0]
    
    
    M_vir_up = data.params_obs["M_vir"]+data.params_obs["M_vir_err_up"]
    v_c_up = get_vc_from_virial_mass(M_vir_up, params["redshift"])/1e5

    M_vir_down = data.params_obs["M_vir"]-data.params_obs["M_vir_err_down"]
    v_c_down = get_vc_from_virial_mass(M_vir_down, params["redshift"])/1e5

    SFR_up = data.params_obs["SFR"]+data.params_obs["SFR_err_up"]
    SFR_down = data.params_obs["SFR"]-data.params_obs["SFR_err_down"]

    
    for i in range(len(betas)):
        
        
            params.update(beta = betas[i])
            
            profiles = get_profiles(params, resol=1000)
            ionization_state = get_ionization_states(profiles, params)
            sigma_CII = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500, add_CMB_suppression=True)    
            intensity_raw = get_intensity_raw(sigma_CII, params, data.params_obs)
            intensity_conv = get_intensity_convolved(intensity_raw, params, data.params_obs, data, add_central_contribution=False)
         
            
            params.update(M_vir = M_vir_up)
            params.update(v_c = v_c_up)
                
            profiles_up = get_profiles(params, resol=1000)
            ionization_state_up = get_ionization_states(profiles_up, params)
            sigma_CII_up = get_surface_density(profiles_up, ionization_state_up, params, rmax=30, h_resol=500, add_CMB_suppression=True)
            intensity_raw_up = get_intensity_raw(sigma_CII_up, params, data.params_obs)
            intensity_conv_up = get_intensity_convolved(intensity_raw_up, params, data.params_obs, data, add_central_contribution=False)

            params.update(M_vir = M_vir_down)
            params.update(v_c = v_c_down)

            profiles_down = get_profiles(params, resol=1000)
            ionization_state_down = get_ionization_states(profiles_down, params)
            sigma_CII_down = get_surface_density(profiles_down, ionization_state_down, params, rmax=30, h_resol=500, add_CMB_suppression=True)
            intensity_raw_down = get_intensity_raw(sigma_CII_down, params, data.params_obs)
            intensity_conv_down = get_intensity_convolved(intensity_raw_down, params, data.params_obs, data, add_central_contribution=False)

            params.update(M_vir = data.params_obs["M_vir"])
            params.update(v_c = data.params_obs["v_c"])

            ax_vc.plot(intensity_conv_up.h/(1000*nc.pc), intensity_conv_up.var, color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))
            ax_vc.plot(intensity_conv_down.h/(1000*nc.pc), intensity_conv_down.var, color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))
            ax_vc.fill_between(intensity_conv.h/(1000*nc.pc), intensity_conv_down.var, intensity_conv_up.var,  alpha=0.2, \
                               color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))

                
            params.update(SFR = SFR_up)

            profiles_up = get_profiles(params, resol=1000)
            ionization_state_up = get_ionization_states(profiles_up, params)
            sigma_CII_up = get_surface_density(profiles_up, ionization_state_up, params, rmax=30, h_resol=500, add_CMB_suppression=True)
            intensity_raw_up = get_intensity_raw(sigma_CII_up, params, data.params_obs)
            intensity_conv_up = get_intensity_convolved(intensity_raw_up, params, data.params_obs, data, add_central_contribution=False)

            params.update(SFR = SFR_down)
            
            profiles_down = get_profiles(params, resol=1000)
            ionization_state_down = get_ionization_states(profiles_down, params)
            sigma_CII_down = get_surface_density(profiles_down, ionization_state_down, params, rmax=30, h_resol=500, add_CMB_suppression=True)
            intensity_raw_down = get_intensity_raw(sigma_CII_down, params, data.params_obs)
            intensity_conv_down = get_intensity_convolved(intensity_raw_down, params, data.params_obs, data, add_central_contribution=False)

            params.update(SFR = data.params_obs["SFR"])

            ax_sfr.plot(intensity_conv_up.h/(1000*nc.pc), intensity_conv_up.var, color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))
            ax_sfr.plot(intensity_conv_down.h/(1000*nc.pc), intensity_conv_down.var, color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))
            ax_sfr.fill_between(intensity_conv.h/(1000*nc.pc), intensity_conv_down.var, intensity_conv_up.var,  alpha=0.2, \
                               color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))


            
    alpine = axs[data_counter].errorbar(data.x/(1000*nc.pc), data.data, yerr=data.err, \
                 markerfacecolor='maroon',markeredgecolor='maroon', marker='o',\
                 linestyle='', ecolor = 'maroon')


    factor_data = data.data[0]/data.beam[0]        
    ax_vc.plot(data.x_beam/1e3/nc.pc, data.beam*factor_data, linestyle="--", color="gray")

    ax_vc.legend([alpine], ["{}".format(data["name_short"])])##
    
    ax_sfr.plot(data.x_beam/1e3/nc.pc, data.beam*factor_data, linestyle="--", color="gray")

    ax_sfr.legend([alpine], ["{}".format(data["name_short"])])##
    
    
    plt.subplots_adjust(left=0.1,  # the left side of the subplots of the figure
                        right=0.98,  # the right side of the subplots of the figure
                        bottom=0.29,  # the bottom of the subplots of the figure
                        top=0.97,  # the top of the subplots of the figure
                        wspace=0.05,  # the amount of width reserved for space between subplots,
                        # expressed as a fraction of the average axis width
                        hspace=0.15)  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height

    cax = plt.axes([0.15, 0.12,  0.72,0.03])
    
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    #betas_ticks = betas[10*betas%2 ==0]
    betas_ticks = betas
    
    
    cb = fig.colorbar(cmap, ticks=betas_ticks.round(1), cax=cax, orientation='horizontal')
    cb.set_label(r'$\eta$', rotation=0.)
    
    plt.show()
    