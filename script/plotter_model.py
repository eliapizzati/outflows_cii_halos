# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 19:36:13 2021

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



plot1 = True  #cc85 + NFW vc
plot2 = False      # cooling
plot3 = False    # cmb suppression theory


if plot1:
    """
    #cc85 + NFW vc
    """
    
    
    # CC85
    
    betas = np.arange(0.2,8,0.8)
    betas = np.append(betas,  [0.4,1.3])
    cmap_rend_col = matplotlib.cm.get_cmap('viridis_r')

    norm = matplotlib.colors.Normalize(vmin=betas.min(), vmax=betas.max())

    fig, [ax_cc, ax_nfw] = plt.subplots(1,2,figsize=(1.3*8.27,1.3*3.2))

    ax_cc.set_xlabel("log ( r / R )")
    ax_cc.set_ylabel("log ( T [K] )")
    ax_cc.set_ylim((4.,8))   
    ax_cc.set_xlim(xmin=-1.0)   

    
    for i in range(len(betas)):
        
        r, v, p, T = cc85_profiles(beta=betas[i], SFR=50)
        
        ax_cc.plot(np.log10(r), np.log10(T), color = cmap_rend_col((betas[i]-betas.min())/(betas.max()-betas.min())))
    
    
    #y = np.linspace(4,8,1000)
    #ax_cc.plot([np.log10(10/0.3)]*len(y), y, linestyle='--', color='gray', alpha=0.5)
    
    
    cax_cc = plt.axes([0.07, 0.15, 0.018,0.8])
    
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, ticks=betas.round(1), cax=cax_cc, orientation='vertical')
    cb.set_label(r'$\eta$', rotation=90., labelpad=-3)
    cb.set_ticks(np.arange(0.5,7.5,1.0))
    cb.set_ticklabels(np.arange(0.5,7.5,1.0))
   
    cax_cc.yaxis.set_ticks_position('left')
    cax_cc.yaxis.set_label_position('left')

    # NFW

    M_vir = np.linspace(1,10,11)
    
    cmap_rend_col = matplotlib.cm.get_cmap('inferno_r')

    norm = matplotlib.colors.Normalize(vmin=M_vir.min(), vmax=M_vir.max())

    radius = np.linspace(1,15,1000)

    for i in range(len(M_vir)):
        v_c = get_circular_velocity_profile_NFW(radius*nc.pc*1e3, M_vir[i]*1e11, z=6)
        v_c_global = get_vc_from_virial_mass( M_vir[i]*1e11, z=6)
        
        ax_nfw.plot(radius, v_c/1e5,  color = cmap_rend_col((M_vir[i]-M_vir.min())/(M_vir.max()-M_vir.min())))
        ax_nfw.plot(radius, [v_c_global/1e5]*len(radius),  color = cmap_rend_col((M_vir[i]-M_vir.min())/(M_vir.max()-M_vir.min())),\
                    linestyle="--")

    ax_nfw.set_ylabel(r"$v_c$ [km/s]")
    ax_nfw.set_xlabel("r [kpc]")

    cax_nfw = plt.axes([0.9, 0.15, 0.018,0.8])
    
    
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, ticks=M_vir, cax=cax_nfw, orientation='vertical')
    cb.set_label(r'$M_{vir}$ [10$^{11}$ M$_\odot$]', rotation=90., labelpad=5)
    cb.set_ticks(np.arange(1.0,10.,1.0))
    cb.set_ticklabels(np.arange(1.0,10.,1.0))



    plt.subplots_adjust(left = 0.15,  # the left side of the subplots of the figure
                        right = 0.885,   # the right side of the subplots of the figure
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

    fig, (ax_gal, ax_no) = plt.subplots(1, 2, sharey=True,figsize=(1.3*8.27,1.3*3.2))
    
    cmap_rend_col = matplotlib.cm.get_cmap('Dark2')
    
    T = np.logspace(3,8)
    n_s = np.logspace(-6,1, 8)
    
    n_min = min(n_s)
    n_max = max(n_s)
    
    norm = matplotlib.colors.LogNorm(vmin=n_min/10**0.5, vmax=n_max*10**0.5)
    
    

    ax_gal.set_xlim((3.2,8.2))
    #plt.yscale('log')
    ax_gal.set_title(r'$f\,_\mathrm{esc}\,=\,0.2$', fontsize = 'x-large', pad = 7)
    ax_gal.set_xlabel("log (T [K])")
    ax_gal.set_ylabel("log ($|\Lambda$(T)| $[\mathrm{erg}\,\mathrm{cm}^3\,\mathrm{s}^{-1}]$)")
    
    ax_no.set_xlim((3.2,8.2))
    ax_no.set_ylim((-26,-21.4))
    
    #plt.yscale('log')
    ax_no.set_title(r'$f\,_\mathrm{esc}\,=\,0.0$', fontsize = 'x-large', pad = 7)
    ax_no.set_xlabel("log (T [K])")
    #ax_no.set_ylabel("log ($|\Lambda$(T)| $[\mathrm{erg}\,\mathrm{cm}^3\,\mathrm{s}^{-1}]$)" , size=18)
    
    for n in n_s:

        # gal        
        
        Plw = 1.4229125141877616e-08 
        Ph1 =  5.48031935502901e-09 # s^-1
        Pg1 = 1.7687762344020628e-09 # s^-1
        
        cooling = cooling_vectorized(T, n, Plw, Ph1, Pg1, 0.) - heating_vectorized(T, n, Plw, Ph1, Pg1, 0.)

        T_min = np.argmin(cooling)

        ax_gal.plot(np.log10(T[T > T[T_min]]), cooling[T > T[T_min]],\
                    color = cmap_rend_col((np.log10(n)-np.log10(n_min))/(np.log10(n_max)-np.log10(n_min))))

        ax_gal.plot(np.log10(T[T < T[T_min]]), cooling[T < T[T_min]], linestyle='--',\
                    color = cmap_rend_col((np.log10(n)-np.log10(n_min))/(np.log10(n_max)-np.log10(n_min))))
    
        # UVB        
        
        Plw = 2.053552516714885e-13
        Ph1 = 2.3335173959061055e-13
        Pg1 = 1.348999532748642e-13
        
        cooling = cooling_vectorized(T, n, Plw, Ph1, Pg1, 0.) - heating_vectorized(T, n, Plw, Ph1, Pg1, 0.)

        T_min = np.argmin(cooling)

        ax_no.plot(np.log10(T[T > T[T_min]]), cooling[T > T[T_min]],\
                    color = cmap_rend_col((np.log10(n)-np.log10(n_min))/(np.log10(n_max)-np.log10(n_min))))

        ax_no.plot(np.log10(T[T < T[T_min]]), cooling[T < T[T_min]], linestyle='--',\
                    color = cmap_rend_col((np.log10(n)-np.log10(n_min))/(np.log10(n_max)-np.log10(n_min))))
    
    ax_gal.plot(np.log10(T), cooling_vectorized(T, 1.,0.,0.,0.,0.)-heating_vectorized(T,1.,0.,0.,0.,0.),color="gray")
    ax_no.plot(np.log10(T), cooling_vectorized(T, 1.,0.,0.,0.,0.)-heating_vectorized(T,1.,0.,0.,0.,0.),color="gray")
    
    
    plt.subplots_adjust(left = 0.1,  # the left side of the subplots of the figure
    right = 0.95,   # the right side of the subplots of the figure
    bottom = 0.28,  # the bottom of the subplots of the figure
    top = 0.9,     # the top of the subplots of the figure
    wspace = 0.05,  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    hspace = 0.05)  # the amount of height reserved for space between subplots,
                  # expressed as a fraction of the average axis height
    
    cax = plt.axes([0.15, 0.13,  0.72,0.03])
    
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col, cax=cax)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, orientation='horizontal')
    cb.set_label(r'n [cm$^{-3}$]', rotation=0.)


if plot3:
    """
    # cmb suppression theory
    """
    
    fig, [ax_eta, ax_texc] = plt.subplots(1,2,figsize=(1.3*8.27,1.3*3.2))
    
    for z in [3,4,5,6]:
        
        redshift_CII_temp = z
        T_CMB_z = nc.CMB_temperature*(1+z)

        T_spin_list = np.linspace(T_CMB_z, 200., 1000)
    
        eta_list = eta_func(T_spin_list, z)
        

        
        ax_eta.plot(np.log10(T_spin_list/T_CMB_z), eta_list, label="z={}".format(z))
    
    #ax_eta.set_xscale("log")
    ax_eta.set_ylim(0.,1.)
    ax_eta.set_xlabel(r"log( T$_{exc}$ / T$_{\rm CMB}$(z) )")
    ax_eta.set_ylabel(r"$\eta$", labelpad = -5)
    ax_eta.axvline(0., linestyle='--', color='gray')
    ax_eta.legend(ncol=2)
    

    ax_texc.set_xlabel(r"log( $n_e$ [cm$^{-3}$] )")
    ax_texc.set_ylabel(r"log( I$_{UV}$ [erg s$^{-1}$cm$^2$Hz$^{-1}$sr$^{-1}$] )")
    
    n_es = np.logspace(-5,3, 10000)
    I_UVs = np.logspace(-12,-19,1000)
    
    nn, II = np.meshgrid(n_es, I_UVs)
    
    T_spin = T_spin_func(n=nn/0.5, T=1e3, I_UV=II, x_e=0.5, z=z)

    z = 6.
    T_CMB_z = nc.CMB_temperature*(1+z)
                    

    im = plt.contourf(np.log10(nn),np.log10(II),np.log10(T_spin/T_CMB_z), cmap=matplotlib.cm.inferno, \
                      norm=matplotlib.colors.PowerNorm(gamma=0.4),\
                      vmin=0., vmax=5,\
                      levels = [0.0,0.05,0.1,0.15,0.2,0.3,0.6,0.8,1.0, 1.4,1.8,2.4,3.0,4.0])
    # Create colorbar
    
    cax = plt.axes([0.9, 0.17, 0.018,0.78])

    cbar = ax_texc.figure.colorbar(im, cax=cax)
    cbar.ax.set_ylabel(r"log( T$_{exc}$ / T$_{\rm CMB}$(z) )", rotation=90., labelpad=0)
    
    #ax_texc.set_title("T_spin dependence on n_e, I_UV; (T = 10^3 K, z = 5.)")

    ax2 = ax_texc.twinx()
    
    habing_g0 = 4*np.pi*I_UVs*(13.6-6.)*nc.ev/nc.hh/nc.cc/5.29e-14
    
    ax2.set_yticks(np.around(np.log10(habing_g0[5::100]),1))
    ax2.set_ylabel(r"log( G$_0 )$")

    ax_texc.text(0.17,0.87,s=r'z=6', fontsize=16, bbox={'facecolor': 'lightgray', 'alpha': 0.8,\
                                                      'boxstyle':"round", 'edgecolor':'none','pad': 0.5},\
                 horizontalalignment='center',verticalalignment='center',\
                 transform=ax_texc.transAxes)

    plt.subplots_adjust(left = 0.07,  # the left side of the subplots of the figure
    right = 0.82,   # the right side of the subplots of the figure
    bottom = 0.16,  # the bottom of the subplots of the figure
    top = 0.95,     # the top of the subplots of the figure
    wspace = 0.3,  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    hspace = 0.05)  # the amount of height reserved for space between subplots,
                  # expressed as a fraction of the average axis height

