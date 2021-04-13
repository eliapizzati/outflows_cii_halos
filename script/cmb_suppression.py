# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 15:48:30 2021

@author: anna
"""


import os
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import ticker, cm, colors


import natconst as nc
from model_classes import load_from_file
from radiation_fields import UVB_rates



#functions

def gamma_ul(T):
    
    if T < 1.5e2:
        return 1.58
    
    elif T >= 1.5e2 and T < 2.5e3:
        return 1.60
    
    elif T >= 2.5e3:
        return 2.11

def gamma_ul_vec(T_vec):
    
    gamma_vec = []
    
    for T in T_vec:
        gamma_vec.append(gamma_ul(T))
            
    return np.asarray(gamma_vec)

def gamma_H_ul(T):
    
    if T < 1.5e2:
        return 10.
    
    elif T >= 1.5e2 and T < 5e2:
        return 15
    
    elif T >= 5e2:
        return 20
    

def gamma_H_ul_vec(T_vec):
    
    gamma_vec = []
    
    for T in T_vec:
        gamma_vec.append(gamma_H_ul(T))
            
    return np.asarray(gamma_vec)


def T_spin_func(n, T, I, x_e, z):
    
    
    
    I_CMB = (2*nc.hh*nc.line_CII_rest_frame**3) /  \
            (nc.cc**2 * (np.exp((nc.hh*nc.line_CII_rest_frame) / (nc.kk*(1.+z)*nc.CMB_temperature)) - 1.))
        
    T_star = 91.7 # K

    A_coeff = 2.36e-6 # s^-1

    B_coeff = A_coeff * nc.cc**2 / (2*nc.hh*nc.line_CII_rest_frame**3)

    C_ul_e = 8.63e-6 / 2. / np.sqrt(T) * gamma_ul(T) 
    
    C_ul_H = gamma_H_ul(T) * 1e-10 / 1.3 

    n_e = x_e * nc.A_H * n
    n_H = (1. - x_e) * nc.A_H * n
    
    num = B_coeff * (I + I_CMB) + A_coeff + n_e * C_ul_e + n_H * C_ul_H
    
    den =  B_coeff * (I + I_CMB) + n_e * C_ul_e * np.exp(-T_star/T) + n_H * C_ul_H * np.exp(-T_star/T)

    T_spin = T_star / np.log(num/den)

    return T_spin


    
def T_spin_func_vec(n_vec, T_vec, I_vec, x_e_vec, z):
    
     
    
    I_CMB = (2*nc.hh*nc.line_CII_rest_frame**3) /  \
            (nc.cc**2 * (np.exp((nc.hh*nc.line_CII_rest_frame) / (nc.kk*(1.+z)*nc.CMB_temperature)) - 1.))
        
    T_star = 91.7 # K

    A_coeff = 2.36e-6 # s^-1

    B_coeff = A_coeff * nc.cc**2 / (2*nc.hh*nc.line_CII_rest_frame**3)

    C_ul_e = 8.63e-6 / 2. / np.sqrt(T_vec) * gamma_ul_vec(T_vec) 
    
    C_ul_H = gamma_H_ul_vec(T_vec) * 1e-10 / 1.3

    n_e = x_e_vec * nc.A_H * n_vec
    n_H = (1. - x_e_vec) * nc.A_H * n_vec
    
    num = B_coeff * (I_vec + I_CMB) + A_coeff + n_e * C_ul_e + n_H * C_ul_H
    
    den =  B_coeff * (I_vec + I_CMB) + n_e * C_ul_e * np.exp(-T_star/T_vec) + n_H * C_ul_H * np.exp(-T_star/T_vec)

    T_spin = T_star / np.log(num/den)

    return T_spin


def eta_func(T_spin, z):
       
    return 1. - (np.exp(nc.hh*nc.line_CII_rest_frame / (nc.kk*(1.+z)*T_spin)) - 1.) /\
           (np.exp((nc.hh*nc.line_CII_rest_frame) / (nc.kk*(1.+z)**2*nc.CMB_temperature)) - 1.)



if __name__ == "__main__":

    
#    params = dict([("SFR", 20.),
#                   ("beta", 1.0), 
#                   ("redshift", 5.),
#                   ("f_esc", 0.), 
#                   ("v_c", 175.),
#                   ("Zeta", 1.0),
#                   ("alfa", 1.0),
#                   ("R_in", 0.3)])
#
#    params_obs = dict([#("name", name),
#                       ("redshift", 5.),
#                       #("line_FWHM", CII_FWHM*nc.km),
#                       #("M_star", 10**log_M_star),
#                       #("SFR", 10**log_SFR),
#                       #("v_c", np.sqrt(nc.gg*M_halo)),
#                       ("sersic_effective_radius", 1.1),
#                       ("sersic_index", 1.)])
#       
#    
#    
#    profiles = load_from_file(params, type_class = "sol_profiles")
#    
#    ionization_state = get_ionization_states(profiles, params)    
#          
#    intensity_tot = intensity_CII_UVB + nc.intensity_CII_1000 * (1000*nc.pc/profiles.r)**2 * params["SFR"]*params["f_esc"]
#
#    redshift_CII = params_obs["redshift"]
#        
#    T_spin = T_spin_func_vec(n_vec=profiles.n, T_vec=profiles.T, I_vec= intensity_tot, x_e_vec = ionization_state.x_e, z=redshift_CII)
#    
#    print(T_spin)
#    
#    eta = eta_func_vec(T_spin_vec=T_spin, z= redshift_CII)
#    
#    print(eta)
    
    
    # FIGURE 1
    
    fig,ax= plt.subplots()
    
    
    for z in [3,4,5,6]:
        
        redshift_CII_temp = z
        T_spin_list = np.linspace(10., 105., 1000)
    
        eta_list = eta_func(T_spin_list, z)
        
        ax.plot(T_spin_list, eta_list, label="z={}".format(z))
        ax.set_ylim(0.,1.)
        ax.set_xlabel("spin temperature [K]")
        ax.set_ylabel("eta")
        ax.axvline(nc.CMB_temperature*(1+redshift_CII_temp), linestyle='--', color='gray')
        ax.legend()
    
    fig.tight_layout()
        
    # FIGURE 2
    
    fig,ax= plt.subplots()

    T_temp = np.logspace(0., 4, 1000)
    
    ax.plot(T_temp, np.exp(-nc.T_star/T_temp))
    ax.set_xlabel("temperature [K]")
    ax.set_ylabel(r"$\exp(-T_*/T)$")
    ax.set_xscale("log")
    #ax.set_xlim(10,15)
    
    fig.tight_layout()
    
    # FIGURE 3     
    
    fig, ax = plt.subplots(figsize=(8.27,5.))

    cmap_rend_col = cm.get_cmap('viridis')
    
    n_temp = np.logspace(-5,3, 10000)
    
    
    temperatures = np.logspace(1.25,3.,10)
    
    T_min = min(temperatures)
    T_max = max(temperatures)
    
    norm = colors.Normalize(vmin=T_min, vmax=T_max)
    
    
    for T in temperatures:

        T_spin_temp = T_spin_func(n_temp, T=T, I=0.,x_e=0.5, z=5.)
    
        eta_temp = eta_func(T_spin_temp, z=5.)
    
        ax.plot(n_temp, eta_temp,  color = cmap_rend_col((T-T_min)/(T_max-T_min)))

    ax.set_xlabel("n [cm-3]")
    ax.set_ylabel("eta")
    ax.set_xscale("log")
    ax.set_yscale("log")
    
    ax.axhline(1., linestyle='--', color='gray')

    ax.set_title("Eta dependence on n, T (I=0, x_e=0.5, z=5.)")    
    
    plt.subplots_adjust(left = 0.15,  # the left side of the subplots of the figure
        right = 0.98,   # the right side of the subplots of the figure
        bottom = 0.15,  # the bottom of the subplots of the figure
        top = 0.9,     # the top of the subplots of the figure
        wspace = 0.1,  # the amount of width reserved for space between subplots,
        # expressed as a fraction of the average axis width
        hspace = 0.1)  # the amount of height reserved for space between subplots,
                      # expressed as a fraction of the average axis height

    
    cmap = cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, orientation='vertical')
    cb.set_label(r'$T$', rotation=0.)

    fig.tight_layout()
    
    # FIGURE 4
    
    fig, ax = plt.subplots(figsize=(8.27,5.))

    cmap_rend_col = cm.get_cmap('Dark2')
    
    n_temp = np.logspace(-5,3, 10000)
    
    
    redshifts = np.linspace(1,8,8)
    
    z_min = min(redshifts)
    z_max = max(redshifts)
    
    norm = colors.Normalize(vmin=z_min-0.5, vmax=z_max+0.5)
    
    
    for z in redshifts:

        T_spin_temp = T_spin_func(n_temp, T=1e3, I=0.,x_e=0.5, z=z)
    
        eta_temp = eta_func(T_spin_temp, z=z)
    
        ax.plot(n_temp, eta_temp,  color = cmap_rend_col((z-z_min)/(z_max-z_min)))

    ax.set_xlabel("n [cm-3]")
    ax.set_ylabel("eta")
    ax.set_xscale("log")
    ax.set_yscale("log")
    
    ax.axhline(1., linestyle='--', color='gray')

    ax.set_title("Eta dependence on n, z (T=1000, I=0, x_e=0.5)")    
    
    plt.subplots_adjust(left = 0.15,  # the left side of the subplots of the figure
        right = 0.98,   # the right side of the subplots of the figure
        bottom = 0.15,  # the bottom of the subplots of the figure
        top = 0.9,     # the top of the subplots of the figure
        wspace = 0.1,  # the amount of width reserved for space between subplots,
        # expressed as a fraction of the average axis width
        hspace = 0.1)  # the amount of height reserved for space between subplots,
                      # expressed as a fraction of the average axis height

    
    cmap = cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, orientation='vertical')
    cb.set_label(r'  $z$', rotation=0.)

    fig.tight_layout()
    
    # FIGURE 5
    
    fig, ax = plt.subplots(figsize=(8.27,5.))

    cmap_rend_col = cm.get_cmap('magma_r')
    
    n_temp = np.logspace(-5,3, 10000)
    
    
    x_es = np.logspace(-2,0,8)
    
    x_e_min = min(x_es)
    x_e_max = max(x_es)
    
    norm = colors.Normalize(vmin=x_e_min, vmax=x_e_max)
    
    
    for x_e in x_es:

        T_spin_temp = T_spin_func(n_temp, T=1e3, I=0.,x_e=x_e, z=5.)
    
        eta_temp = eta_func(T_spin_temp, z=5.)
    
        ax.plot(n_temp, eta_temp,  color = cmap_rend_col((x_e-x_e_min)/(x_e_max-x_e_min)))

    ax.set_xlabel("n [cm-3]")
    ax.set_ylabel("eta")
    ax.set_xscale("log")
    ax.set_yscale("log")
    
    ax.axhline(1., linestyle='--', color='gray')

    ax.set_title("Eta dependence on n, x_e (T=1000, I=0, z=5)")    
    
    plt.subplots_adjust(left = 0.15,  # the left side of the subplots of the figure
        right = 0.98,   # the right side of the subplots of the figure
        bottom = 0.15,  # the bottom of the subplots of the figure
        top = 0.9,     # the top of the subplots of the figure
        wspace = 0.1,  # the amount of width reserved for space between subplots,
        # expressed as a fraction of the average axis width
        hspace = 0.1)  # the amount of height reserved for space between subplots,
                      # expressed as a fraction of the average axis height

    
    cmap = cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, orientation='vertical')
    cb.set_label(r'  $x_e$', rotation=0.)

    fig.tight_layout()
    
    
    # FIGURE 6
    
    fig, ax = plt.subplots(figsize=(8.27,5.))

    cmap_rend_col = cm.get_cmap('tab10_r')
    
    n_temp = np.logspace(-5,3, 10000)
    
    
    I_s = np.logspace(-13,-22,10)
    
    I_min = min(I_s)
    I_max = max(I_s)
    
    norm = colors.LogNorm(vmin=I_min/10**0.5, vmax=I_max*10**0.5)
    
    
    for I in I_s:
        print(I)

        T_spin_temp = T_spin_func(n_temp, T=1e3, I=I,x_e=0.5, z=5.)
    
        eta_temp = eta_func(T_spin_temp, z=5.)
    
        ax.plot(n_temp, eta_temp,  color = cmap_rend_col((np.log10(I)-np.log10(I_min))/(np.log10(I_max)-np.log10(I_min))))

    ax.set_xlabel("n [cm-3]")
    ax.set_ylabel("eta")
    ax.set_xscale("log")
    ax.set_yscale("log")
    
    ax.axhline(1., linestyle='--', color='gray')

    ax.set_title("Eta dependence on n, I (T=1000, x_e=0.5, z=5)")    
    
    plt.subplots_adjust(left = 0.15,  # the left side of the subplots of the figure
        right = 0.98,   # the right side of the subplots of the figure
        bottom = 0.15,  # the bottom of the subplots of the figure
        top = 0.9,     # the top of the subplots of the figure
        wspace = 0.1,  # the amount of width reserved for space between subplots,
        # expressed as a fraction of the average axis width
        hspace = 0.1)  # the amount of height reserved for space between subplots,
                      # expressed as a fraction of the average axis height

    
    cmap = cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, orientation='vertical')
    cb.set_label(r'  $I$', rotation=0.)

    fig.tight_layout()
    








