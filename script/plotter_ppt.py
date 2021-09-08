# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 11:38:20 2021

@author: anna
"""




import numpy as np
import matplotlib.pyplot as plt
import natconst as nc


from cc85 import cc85_profiles
from my_utils import get_circular_velocity_profile_NFW, get_vc_from_virial_mass
from cmb_suppression import eta_func, T_spin_func

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



plot1 = False  #cc85 
plot2 = False      # cooling
plot3 = False    # cmb suppression theory
plot4 = True    # profiles
plot5 = False   # ionization

if plot1:
    """
    #cc85
    """
    
    
    # CC85
    
    betas = np.arange(0.2,8,0.8)
    betas = np.append(betas,  [0.4,1.3])
    cmap_rend_col = matplotlib.cm.get_cmap('viridis_r')

    norm = matplotlib.colors.Normalize(vmin=betas.min(), vmax=betas.max())

    fig, ax_cc = plt.subplots(1,1,figsize=(1.3*4.93,1.3*3.8))

    ax_cc.set_xlabel("log ( r / R )")
    ax_cc.set_ylabel("log ( T [K] )")
    ax_cc.set_ylim((4.5,8))   
    ax_cc.set_xlim(xmin=-1.0)   

    
    for i in range(len(betas)):
        
        r, v, p, T = cc85_profiles(beta=betas[i], SFR=50)
        
        ax_cc.plot(np.log10(r), np.log10(T), color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))
    
    
    #y = np.linspace(4,8,1000)
    #ax_cc.plot([np.log10(10/0.3)]*len(y), y, linestyle='--', color='gray', alpha=0.5)
    
    
    cax_cc = plt.axes([0.12, 0.15, 0.025,0.8])
    
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, ticks=betas.round(1), cax=cax_cc, orientation='vertical')
    cb.set_label(r'$\eta$', rotation=90., labelpad=-3)
    cb.set_ticks(np.arange(0.5,7.5,1.0))
    cb.set_ticklabels(np.arange(0.5,7.5,1.0))
   
    cax_cc.yaxis.set_ticks_position('left')
    cax_cc.yaxis.set_label_position('left')
    
    

    plt.subplots_adjust(left = 0.28,  # the left side of the subplots of the figure
                        right = 0.95,   # the right side of the subplots of the figure
                        bottom = 0.15,  # the bottom of the subplots of the figure
                        top = 0.95,     # the top of the subplots of the figure
                        wspace = 0.28,  # the amount of width reserved for space between subplots,
                        # expressed as a fraction of the average axis width
                        hspace = 0.1)  # the amount of height reserved for space between subplots,
                      # expressed as a fraction of the average axis height


if plot2:
    """
    # cooling
    """
    from cooling import cooling_vectorized, heating_vectorized

    fesc = 0.0
    
    fig, ax = plt.subplots(1, 1, figsize=(1.3*4.93,1.3*3.8))
    
    cmap_rend_col = matplotlib.cm.get_cmap('Dark2_r')
    
    T = np.logspace(3,8,1000)
    n_s = np.logspace(-6,1, 8)
    
    n_min = min(n_s)
    n_max = max(n_s)
    
    norm = matplotlib.colors.Normalize(vmin=np.log10(n_min)-0.5, vmax=np.log10(n_max)+0.5)
    
    

    ax.set_xlim((3.2,8.2))
    #plt.yscale('log')
    #ax_gal.set_title(r'$f\,_\mathrm{esc}\,=\,0.2$', pad = 7)
    ax.set_xlabel("log (T [K])")
    ax.set_ylabel("log ($|\Lambda$(T)| $[\mathrm{erg}\,\mathrm{cm}^3\,\mathrm{s}^{-1}]$)")
    
    ax.set_xlim((3.2,8.2))
    ax.set_ylim((-26,-21.4))
    
    
    for n in n_s:

        # gal        
        
        
        if fesc == 0.2:
            
            alpha_gal = 1.
            alpha_uvb = 0.1
        
        elif fesc == 0.0:
             
            alpha_gal = 0.1
            alpha_uvb = 1.
                 
            
        Plw = 1.4229125141877616e-08 
        Ph1 =  5.48031935502901e-09 # s^-1
        Pg1 = 1.7687762344020628e-09 # s^-1
        
        T, cooling = cooling_vectorized(T, n, Plw, Ph1, Pg1, 0.)
        T, heating = heating_vectorized(T, n, Plw, Ph1, Pg1, 0.)
    
        net_cooling = np.abs(cooling - heating)
        T_min = np.argmin(net_cooling)
    
        ax.plot(np.log10(T[T > T[T_min-1]]), np.log10(net_cooling[T > T[T_min-1]]),\
                    color = cmap_rend_col((np.log10(n)-np.log10(n_min))/(np.log10(n_max)-np.log10(n_min))), \
                    alpha=alpha_gal)
    
        ax.plot(np.log10(T[T < T[T_min-1]]), np.log10(net_cooling[T < T[T_min-1]]), linestyle='--',\
                    color = cmap_rend_col((np.log10(n)-np.log10(n_min))/(np.log10(n_max)-np.log10(n_min))),\
                    alpha=alpha_gal)


        # UVB        
        
        Plw = 2.053552516714885e-13
        Ph1 = 2.3335173959061055e-13
        Pg1 = 1.348999532748642e-13
        
        T, cooling = cooling_vectorized(T, n, Plw, Ph1, Pg1, 0.)
        T, heating = heating_vectorized(T, n, Plw, Ph1, Pg1, 0.)

        net_cooling = np.abs(cooling - heating)

        T_min = np.argmin(net_cooling)

        ax.plot(np.log10(T[T > T[T_min]]), np.log10(net_cooling[T > T[T_min]]),\
                    color = cmap_rend_col((np.log10(n)-np.log10(n_min))/(np.log10(n_max)-np.log10(n_min))),\
                    alpha=alpha_uvb)

        ax.plot(np.log10(T[T < T[T_min]]), np.log10(net_cooling[T < T[T_min]]), linestyle="--",\
                    color = cmap_rend_col((np.log10(n)-np.log10(n_min))/(np.log10(n_max)-np.log10(n_min))),\
                    alpha=alpha_uvb)

    net_cooling_cie = np.abs(cooling_vectorized(T, 1.,0.,0.,0.,0.)[1]-heating_vectorized(T,1.,0.,0.,0.,0.)[1])
    print(net_cooling_cie)
    ax.plot(np.log10(T), np.log10(net_cooling_cie),color="gray", linestyle=":")
    
    ax.text(0.19,0.87,s=r'$f_{{esc}}$={:.1f}'.format(fesc), fontsize=16, bbox={'facecolor': 'lightgray', 'alpha': 0.8,\
                                                      'boxstyle':"round", 'edgecolor':'none','pad': 0.5},\
                 horizontalalignment='center',verticalalignment='center',\
                 transform=ax.transAxes)

    
 

    cax = plt.axes([0.82, 0.15, 0.025,0.8])
    
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, orientation='vertical',  cax=cax)
    cb.set_label(r'log( n [cm$^{-3}$] )', rotation=90.)
    cb.set_ticks(np.arange(-6.0,2.0,1.0))
    cb.set_ticklabels(np.arange(-6.0,2.0,1.0))
    cb.minorticks_off()
    
    plt.subplots_adjust(left = 0.17,  # the left side of the subplots of the figure
                        right = 0.8,   # the right side of the subplots of the figure
                        bottom = 0.15,  # the bottom of the subplots of the figure
                        top = 0.95,     # the top of the subplots of the figure
                        wspace = 0.28,  # the amount of width reserved for space between subplots,
                        # expressed as a fraction of the average axis width
                        hspace = 0.1)  # the amount of height reserved for space between subplots,
                      # expressed as a fraction of the average axis height

    plt.show()



if plot3:
    """
    cmb theory
    """
    
    pass


if plot4:
    """
    profiles
    """
#    from sol_modules import get_profiles

    fesc = 0.0
    
    if fesc == 0.2:
            
        alpha_gal = 1.
        alpha_uvb = 0.1
        
    elif fesc == 0.0:
             
        alpha_gal = 0.1
        alpha_uvb = 1.
        
        

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

     
    betas = np.arange(0.2,8,0.8)
    betas = np.append(betas,  [0.4,1.,1.8])
    cmap_rend_col = matplotlib.cm.get_cmap('viridis_r')
    
    norm = matplotlib.colors.Normalize(vmin=betas.min(), vmax=betas.max())

    
    fig, [ax_v,ax_n,ax_T] = plt.subplots(1, 3, sharex=True, figsize=(1.5*12.27,1.5*4.0))



    # plot preparation
  

    ax_v.set_xlabel("log (r [kpc])")
    ax_v.set_ylabel("v [1000 km/s] ")
    #ax_v.set_title(r'$f\,_\mathrm{esc}\,=\,0.0$', fontsize = 'x-large', pad = 7)
    #ax_v.tight_layout()
    ax_v.set_xlim((np.log10(0.3),np.log10(30)))

    ax_n.set_xlabel("log (r [kpc])")
    ax_n.set_ylabel("log (n [cm$^{-3}$]) ")
    
    ax_T.set_xlabel("log (r [kpc])")
    ax_T.set_ylabel("log (T [K])")
    ax_T.set_ylim((2,8))
    
    
    for i in range(len(betas)):
        print("beta=", betas[i])
        params.update(f_esc_ion = 0.)
        params.update(f_esc_FUV = 0.)

        params.update(beta = betas[i])
     
#        profiles = get_profiles(params, resol=1000)
#
#        r = profiles.r
#        n = profiles.n
#        T = profiles.T
#        v = profiles.v
#        
#        ax_v.plot(np.log10(r/(1000*nc.pc)),v/10**8, color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())), alpha=alpha_uvb)
#
#        ax_n.plot(np.log10(r/(1000*nc.pc)),np.log10(n), color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())), alpha=alpha_uvb)
#
#        ax_T.plot(np.log10(r/(1000*nc.pc)),np.log10(T), color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())), alpha=alpha_uvb)
#
#        params.update(f_esc_ion = 0.2)
#        params.update(f_esc_FUV = 0.2)
#        
#        profiles = get_profiles(params, resol=1000)
#
#        r = profiles.r
#        n_gal = profiles.n
#        T_gal = profiles.T
#        v_gal = profiles.v
#
#        
#        ax_v.plot(np.log10(r/(1000*nc.pc)),v_gal/10**8, color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())), alpha=alpha_gal)
#       
#        ax_n.plot(np.log10(r/(1000*nc.pc)),np.log10(n_gal), color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())), alpha=alpha_gal)
#           
#        ax_T.plot(np.log10(r/(1000*nc.pc)),np.log10(T_gal), color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())), alpha=alpha_gal)
        
    

    plt.subplots_adjust(left = 0.06,  # the left side of the subplots of the figure
                        right = 0.97,   # the right side of the subplots of the figure
                        bottom = 0.29,  # the bottom of the subplots of the figure
                        top = 0.96,     # the top of the subplots of the figure
                        wspace = 0.18,  # the amount of width reserved for space between subplots,
                        # expressed as a fraction of the average axis width
                        hspace = 0.05)  # the amount of height reserved for space between subplots,
                        # expressed as a fraction of the average axis height

    ax_T.text(0.79,0.87,s=r'$f_{{esc}}$={:.1f}'.format(fesc), fontsize=16, bbox={'facecolor': 'lightgray', 'alpha': 0.8,\
                                                      'boxstyle':"round", 'edgecolor':'none','pad': 0.5},\
                 horizontalalignment='center',verticalalignment='center',\
                 transform=ax_T.transAxes)


    cax = plt.axes([0.15, 0.12,  0.72,0.028])
    
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, cax=cax, orientation='horizontal')
    cb.set_label(r'$\eta$', rotation=0., labelpad=-2)
    cb.set_ticks(np.arange(1.0,10.,1.0))
    cb.set_ticklabels(np.arange(1.0,10.,1.0))
    plt.show()



if plot5:
    """
    ionization
    """
    
    
#    from sol_modules import get_profiles
#    from post_sol_modules import get_ionization_states

    fesc = 0.0
    
    if fesc == 0.2:
            
        alpha_gal = 1.
        alpha_uvb = 0.1
        
    elif fesc == 0.0:
             
        alpha_gal = 0.1
        alpha_uvb = 1.
        
        

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

     
    betas = np.arange(0.2,8,0.8)
    betas = np.append(betas,  [0.4,1.,1.8])
    cmap_rend_col = matplotlib.cm.get_cmap('viridis_r')
    
    norm = matplotlib.colors.Normalize(vmin=betas.min(), vmax=betas.max())

    

        
    fig, axs = plt.subplots(1, 2, sharex=True, figsize=(1.3*8.27,1.3*2.9))

    ax_xe = axs[0]
    ax_xCII = axs[1]
    
    
    
    ax_xe.set_xlabel("log (r [kpc])")
    ax_xe.set_ylabel("n$_\\mathrm{HI}$/n$_\\mathrm{H}$")
    ax_xe.set_ylim((0,1))
    ax_xe.set_xlim((np.log10(0.3),np.log10(30)))
        
    ax_xCII.set_xlabel("log (r [kpc])")
    ax_xCII.set_ylabel("n$_\\mathrm{CII}$/n$_\\mathrm{C}$")
    ax_xCII.ticklabel_format(axis='y', style='plain')
    ax_xCII.set_ylim((0,1))
    
    for i in range(len(betas)):
        print("beta=", betas[i])

        params.update(f_esc_ion = 0.)
        params.update(f_esc_FUV = 0.)

        params.update(beta = betas[i])

        
#        profiles = get_profiles(params, resol=1000)
#        ionization_state = get_ionization_states(profiles, params)
#
#        r = ionization_state.r
#        x_e = ionization_state.x_e
#        x_CII = ionization_state.x_CII
#        
#        ax_xCII.plot(np.log10(r/(1000*nc.pc)),x_CII, color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())), alpha=alpha_uvb)
#
#        ax_xe.plot(np.log10(r/(1000*nc.pc)),1.-x_e, color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())), alpha=alpha_uvb)
#
#        params.update(f_esc_ion = 0.2)
#        params.update(f_esc_FUV = 0.2)
#        
#        profiles = get_profiles(params, resol=1000)
#        ionization_state = get_ionization_states(profiles, params)
#
#        r = ionization_state.r
#        x_e_gal = ionization_state.x_e
#        x_CII_gal = ionization_state.x_CII
#        
#    
#        ax_xe_gal.plot(np.log10(r/(1000*nc.pc)),np.log10(1.-x_e_gal), color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())), alpha=alpha_gal)
#    
#        ax_xCII_gal.plot(np.log10(r/(1000*nc.pc)),np.log10(x_CII_gal), color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())), alpha=alpha_gal)

    
    plt.subplots_adjust(left = 0.08,  # the left side of the subplots of the figure
    right = 0.89,   # the right side of the subplots of the figure
    bottom = 0.17,  # the bottom of the subplots of the figure
    top = 0.95,     # the top of the subplots of the figure
    wspace = 0.24,  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    hspace = 0.15)  # the amount of height reserved for space between subplots,
                  # expressed as a fraction of the average axis height
    
    
    ax_xCII.text(0.79,0.87,s=r'$f_{{esc}}$={:.1f}'.format(fesc), fontsize=16, bbox={'facecolor': 'lightgray', 'alpha': 0.8,\
                                                      'boxstyle':"round", 'edgecolor':'none','pad': 0.5},\
                 horizontalalignment='center',verticalalignment='center',\
                 transform=ax_xCII.transAxes)

    cax = plt.axes([0.91,  0.17, 0.015,0.78])
    
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, ticks=betas, cax=cax, orientation='vertical')
    cb.set_label(r'$\eta$',rotation=90.)
    cb.set_ticks(np.arange(1.0,10.,1.0))
    cb.set_ticklabels(np.arange(1.0,10.,1.0))

    plt.show()
    
    