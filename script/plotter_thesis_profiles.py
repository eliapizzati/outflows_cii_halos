


"""
This script creates some basic figures for the thesis (model chapter) such as v,n,T profiles, ionization states, cmb emission,
final intensity and surface density profiles
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
        "xtick.minor.size": 3.0,
        "ytick.minor.size": 3.0,
        "legend.fontsize": 13.0,
        "legend.frameon": False
 #       "figure.dpi": 200,
 #       "savefig.dpi": 200,
#        "text.usetex": True
})



plot1 = False  # profiles v,n,T
plot2 = False  # ionization state
plot3 = True   # cmb emission
plot4 = False  # final emission



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
    #profiles v,n,T
    """
    from sol_modules import get_profiles

     
    betas = np.arange(0.2,8,0.8)
    betas = np.append(betas,  [0.4,1.,1.8])
    cmap_rend_col = matplotlib.cm.get_cmap('viridis_r')
    
    norm = matplotlib.colors.Normalize(vmin=betas.min(), vmax=betas.max())


    fig, axs = plt.subplots(3, 2, sharex=True,sharey=False, figsize=(1.3*8.27,1.3*10.))


    ax_v_gal = axs[0,0] 
    ax_n_gal = axs[1,0]
    ax_T_gal = axs[2,0]
    ax_v = axs[0,1]
    ax_n = axs[1,1]
    ax_T = axs[2,1]


 
    # plot preparation
      
    #ax_v_gal.set_xlabel("log (r [kpc])", size=size)
    ax_v_gal.set_ylabel("v [1000 km/s] ")
    ax_v_gal.set_title(r'$f\,_\mathrm{esc}\,=\,0.2$', pad = 7)
    #ax_v_gal.tight_layout()
    ax_v_gal.set_xlim((np.log10(0.3),np.log10(30)))
    
    #ax_v.set_xlabel("log (r [kpc])", size=size)
    #ax_v.set_ylabel("v [1000 km/s] ", size=size)
    ax_v.set_title(r'$f\,_\mathrm{esc}\,=\,0.0$', pad = 7)
    #ax_v.tight_layout()
    ax_v.set_xlim((np.log10(0.3),np.log10(30)))
    
    
    #ax_n_gal.set_xlabel("log (r [kpc])", size=size)
    ax_n_gal.set_ylabel("log (n [cm$^{-3}$]) ")
    #ax_n_gal.tight_layout()
    

    #ax_n.set_xlabel("log (r [kpc])", size=size)
    #ax_n.set_ylabel("log (n [cm$^{-3}$]) ", size=size)
    #ax_n.tight_layout()
    
    ax_T_gal.set_xlabel("log (r [kpc])")
    ax_T_gal.set_ylabel("log (T [K])")
    ax_T_gal.set_ylim((1.7,8))
    #ax_T_gal.tight_layout()
    
    ax_T.set_xlabel("log (r [kpc])")
    #ax_T.set_ylabel("log (T [K])", size=size)
    ax_T.set_ylim((1.7,8))
    #ax_T.tight_layout()


    for i in range(len(betas)):
        print("beta=", betas[i])
        params.update(f_esc_ion = 0.)
        params.update(f_esc_FUV = 0.)

        params.update(beta = betas[i])

        
        profiles = get_profiles(params, resol=1000)

        r = profiles.r
        n = profiles.n
        T = profiles.T
        v = profiles.v
        
        ax_v.plot(np.log10(r/(1000*nc.pc)),v/10**8, color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))

        ax_n.plot(np.log10(r/(1000*nc.pc)),np.log10(n), color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))

        ax_T.plot(np.log10(r/(1000*nc.pc)),np.log10(T), color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))

        params.update(f_esc_ion = 0.2)
        params.update(f_esc_FUV = 0.2)
        
        profiles = get_profiles(params, resol=1000)

        r = profiles.r
        n_gal = profiles.n
        T_gal = profiles.T
        v_gal = profiles.v

        
        ax_v_gal.plot(np.log10(r/(1000*nc.pc)),v_gal/10**8, color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))
       
        ax_n_gal.plot(np.log10(r/(1000*nc.pc)),np.log10(n_gal), color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))
           
        ax_T_gal.plot(np.log10(r/(1000*nc.pc)),np.log10(T_gal), color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))
        


    plt.subplots_adjust(left = 0.08,  # the left side of the subplots of the figure
            right = 0.97,   # the right side of the subplots of the figure
            bottom = 0.22,  # the bottom of the subplots of the figure
            top = 0.95,     # the top of the subplots of the figure
            wspace = 0.12,  # the amount of width reserved for space between subplots,
            # expressed as a fraction of the average axis width
            hspace = 0.15)  # the amount of height reserved for space between subplots,
                          # expressed as a fraction of the average axis height


    cax = plt.axes([0.15, 0.1,  0.72,0.015])
    
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, cax=cax, orientation='horizontal')
    cb.set_label(r'$\eta$', rotation=0.)
    cb.set_ticks(np.arange(1.0,10.,1.0))
    cb.set_ticklabels(np.arange(1.0,10.,1.0))
    plt.show()


if plot2:
    """
    #ionization state
    """
    from sol_modules import get_profiles
    from post_sol_modules import get_ionization_states

        
    betas = np.arange(0.2,8,0.8)
    betas = np.append(betas,  [0.4,1.,1.8])
    cmap_rend_col = matplotlib.cm.get_cmap('viridis_r')
    
    norm = matplotlib.colors.Normalize(vmin=betas.min(), vmax=betas.max())

        
    fig, axs = plt.subplots(2, 2, sharex=True, figsize=(1.3*8.27,1.3*5.5))
    
    ax_xe_gal = axs[0,0] 
    ax_xCII_gal = axs[1,0]
    ax_xe = axs[0,1]
    ax_xCII = axs[1,1]
    
    
    
    #ax_xe.set_xlabel("log (r [kpc])", size=size)
    ax_xe.set_ylabel("n$_\\mathrm{HI}$/n$_\\mathrm{H}$")
    ax_xe.set_title(r'$f\,_\mathrm{esc}\,=\,0.0$', pad = 7)
    ax_xe.set_ylim((0,1))
    ax_xe.set_xlim((np.log10(0.3),np.log10(30)))
    
    #ax_xe_gal.set_xlabel("log (r [kpc])", size=size)
    ax_xe_gal.set_ylabel("log (n$_\\mathrm{HI}$/n$_\\mathrm{H})$")
    ax_xe_gal.set_title(r'$f\,_\mathrm{esc}\,=\,0.2$', pad = 7)
    ax_xe_gal.set_ylim((-9.6,0))
    
    
    ax_xCII_gal.set_xlabel("log (r [kpc])")
    ax_xCII_gal.set_ylabel("log (n$_\\mathrm{CII}$/n$_\\mathrm{C})$")
    ax_xCII_gal.ticklabel_format(axis='y', style='plain')
    ax_xCII_gal.set_ylim((-5.2,0))
    
    ax_xCII.set_xlabel("log (r [kpc])")
    ax_xCII.set_ylabel("n$_\\mathrm{CII}$/n$_\\mathrm{C}$")
    ax_xCII.ticklabel_format(axis='y', style='plain')
    ax_xCII.set_ylim((0,1))
    
    for i in range(len(betas)):
        print("beta=", betas[i])

        params.update(f_esc_ion = 0.)
        params.update(f_esc_FUV = 0.)

        params.update(beta = betas[i])

        
        profiles = get_profiles(params, resol=1000)
        ionization_state = get_ionization_states(profiles, params)

        r = ionization_state.r
        x_e = ionization_state.x_e
        x_CII = ionization_state.x_CII
        
        ax_xCII.plot(np.log10(r/(1000*nc.pc)),x_CII, color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))

        ax_xe.plot(np.log10(r/(1000*nc.pc)),1.-x_e, color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))

        params.update(f_esc_ion = 0.2)
        params.update(f_esc_FUV = 0.2)
        
        profiles = get_profiles(params, resol=1000)
        ionization_state = get_ionization_states(profiles, params)

        r = ionization_state.r
        x_e_gal = ionization_state.x_e
        x_CII_gal = ionization_state.x_CII
        
    
        ax_xe_gal.plot(np.log10(r/(1000*nc.pc)),np.log10(1.-x_e_gal), color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))
    
        ax_xCII_gal.plot(np.log10(r/(1000*nc.pc)),np.log10(x_CII_gal), color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))

    
    plt.subplots_adjust(left = 0.08,  # the left side of the subplots of the figure
    right = 0.97,   # the right side of the subplots of the figure
    bottom = 0.22,  # the bottom of the subplots of the figure
    top = 0.95,     # the top of the subplots of the figure
    wspace = 0.22,  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    hspace = 0.15)  # the amount of height reserved for space between subplots,
                  # expressed as a fraction of the average axis height
    
    
    cax = plt.axes([0.15, 0.09,  0.72,0.017])
    
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, ticks=betas, cax=cax, orientation='horizontal')
    cb.set_label(r'$\eta$',rotation=0.)
    cb.set_ticks(np.arange(1.0,10.,1.0))
    cb.set_ticklabels(np.arange(1.0,10.,1.0))

    plt.show()

if plot3:
    """
    #cmb emission
    """

    from sol_modules import get_profiles
    from post_sol_modules import get_ionization_states, get_surface_density

        
    betas = np.arange(0.2,8,0.8)
    betas = np.append(betas,  [0.4,1.3])
    cmap_rend_col = matplotlib.cm.get_cmap('viridis_r')
    
    norm = matplotlib.colors.Normalize(vmin=betas.min(), vmax=betas.max())

    fig, [ax_eta, ax_sigma] = plt.subplots(1,2,figsize=(1.3*8.27,1.3*4.))

    ax_eta.set_xlabel("b [kpc]")
    ax_eta.set_ylabel(r"$\zeta$",  labelpad=-8)
    ax_eta.set_ylim((3e-2,1))
    ax_eta.set_xlim((0.3,15.5))
    ax_eta.set_yscale("log")
    ax_eta.set_yticks([0.1,1])
    ax_eta.set_yticklabels([0.1,1])

    ax_sigma.set_xlabel("b [kpc]")
    ax_sigma.set_ylabel(r"log ($\Sigma_{\rm CII}$ [erg/cm s$^2$])", labelpad=-1)
    ax_sigma.set_yscale("log")
    ax_sigma.set_xlim((0.3, 15.5))
    ax_sigma.set_ylim(1e-8,2e-2)



    for i in range(len(betas)):
        print("beta=", betas[i])

        params.update(beta = betas[i])

        
        profiles = get_profiles(params, resol=1000)
        ionization_state = get_ionization_states(profiles, params)
        sigma_CII = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500, add_CMB_suppression=True)    
        sigma_CII_no_cmb = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500, add_CMB_suppression=False)


        ax_eta.plot(sigma_CII.h/1e3/nc.pc, sigma_CII.eta,  color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))
        ax_sigma.plot(sigma_CII.h/(1e3*nc.pc),sigma_CII.var, color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))
        ax_sigma.plot(sigma_CII_no_cmb.h/(1e3*nc.pc),sigma_CII_no_cmb.var, color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())), alpha=0.4)


    
    plt.subplots_adjust(left = 0.06,  # the left side of the subplots of the figure
    right = 0.98,   # the right side of the subplots of the figure
    bottom = 0.28,  # the bottom of the subplots of the figure
    top = 0.96,     # the top of the subplots of the figure
    wspace = 0.25,  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    hspace = 0.27)  # the amount of height reserved for space between subplots,
                  # expressed as a fraction of the average axis height
    
    
    cax = plt.axes([0.15, 0.12,  0.72,0.015*10/7])
    
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, ticks=betas, cax=cax, orientation='horizontal')
    cb.set_label(r'$\eta$',rotation=0.)
    cb.set_ticks(np.arange(1.0,10.,1.0))
    cb.set_ticklabels(np.arange(1.0,10.,1.0))


    plt.show()

if plot4:
    """
    #final emission
    """
    from sol_modules import get_profiles
    from post_sol_modules import get_ionization_states, get_surface_density, get_intensity_raw, get_intensity_convolved, get_luminosity_CII, get_halo_mass

    from load_data import observational_data_fuji

    betas = np.arange(1.0,9.0,0.8)

    data = observational_data_fuji

    cmap_rend_col = matplotlib.cm.get_cmap('tab10')
    
    norm = matplotlib.colors.Normalize(vmin=betas.min()-0.4, vmax=betas.max()+0.4)
    
    
    fig, (ax_em, ax_comp) = plt.subplots(1, 2, sharey=True, figsize=(1.3*8.27,1.3*4.))
    
    
    ax_em.set_xlabel("b [kpc]")
    ax_em.set_ylabel(r"flux [mJy/arcsec$^2$]")
    ax_em.set_yscale('log')
    ax_em.set_ylim((0.01,12))    
    ax_em.set_xlim((0.3,10))
    
    
    ax_comp.set_xlabel("b [kpc]")
    #ax_comp.set_ylabel("flux [mJy/arcsec$^2$]", size=size)
    ax_comp.set_yscale('log')
    ax_comp.set_xlim((0.3,10))
    
    
    
    for i in range(len(betas)):
        print("beta=", betas[i])

        params.update(f_esc_ion = 0.)
        params.update(f_esc_FUV = 0.)

        params.update(beta = betas[i])

        
        profiles = get_profiles(params, resol=1000)
        ionization_state = get_ionization_states(profiles, params)
        sigma_CII = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500, add_CMB_suppression=True)    
        intensity_raw = get_intensity_raw(sigma_CII, params, data.params_obs)
        intensity_conv = get_intensity_convolved(intensity_raw, params, data.params_obs, data, add_central_contribution=False)
     
        print("L_CII=", get_luminosity_CII(sigma_CII))
        print("M_halo=", get_halo_mass(profiles, params))
        ax_em.plot(intensity_raw.h/(1e3*nc.pc),intensity_raw.var, color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))
        
        ax_comp.plot(intensity_conv.h/(1e3*nc.pc),intensity_conv.var, color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))
                    
            
    fuji = ax_comp.errorbar(data.x/(1000*nc.pc), data.data, yerr=data.err, \
                 markerfacecolor='maroon',markeredgecolor='maroon', marker='o',\
                 linestyle='', ecolor = 'maroon')
    

    # central contribution

    factor_data = data.data[0]/data.beam[0]
    
    ax_comp.plot(data.x_beam/1e3/nc.pc, data.beam*factor_data, linestyle="--", color="gray")

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
    
    ax_comp.legend([fuji], ["Fujimoto+19"])##

    #ax_comp.legend(loc='upper right', fontsize = 'x-large', frameon=False )
    #ax_em.legend(loc='upper right', fontsize = 'large')
    
    plt.show()
    
