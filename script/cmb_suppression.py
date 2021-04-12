# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 15:48:30 2021

@author: anna
"""


import os
import numpy as np
import matplotlib.pyplot as plt 

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
    


def T_spin_func(n, T, I, x_e, z):
    
    
    line_CII_rest_frame = 1900*1e9 
    
    I_CMB = (2*nc.hh*line_CII_rest_frame**3) /  \
            (nc.cc**2 * (np.exp((nc.hh*line_CII_rest_frame) / (nc.kk*(1.+z)*nc.CMB_temperature)) - 1.))
        
    T_star = 91.7 # K

    A_coeff = 2.36e-6 # s^-1

    B_coeff = A_coeff * nc.cc**2 / (2*nc.hh*line_CII_rest_frame**3)

    C_ul_e = 8.63e-6 / 2. / np.sqrt(T) * gamma_ul(T) 
    
    C_ul_H = 1e-6  

    n_e = x_e * nc.A_H * n
    n_H = (1. - x_e) * nc.A_H * n
    
    num = B_coeff * (I + I_CMB) + A_coeff + n_e * C_ul_e + n_H * C_ul_H
    
    den =  B_coeff * (I + I_CMB) + n_e * C_ul_e * np.exp(-T_star/T) + n_H * C_ul_H * np.exp(-T_star/T)

    T_spin = T_star / np.log(num/den)

    return T_spin


    
def T_spin_func_vec(n_vec, T_vec, I_vec, x_e_vec, z):
    
    
    line_CII_rest_frame = 1900*1e9 
    
    I_CMB = (2*nc.hh*line_CII_rest_frame**3) /  \
            (nc.cc**2 * (np.exp((nc.hh*line_CII_rest_frame) / (nc.kk*(1.+z)*nc.CMB_temperature)) - 1.))
        
    T_star = 91.7 # K

    A_coeff = 2.36e-6 # s^-1

    B_coeff = A_coeff * nc.cc**2 / (2*nc.hh*line_CII_rest_frame**3)

    C_ul_e = 8.63e-6 / 2. / np.sqrt(T_vec) * gamma_ul_vec(T_vec) 
    
    C_ul_H = 1e-6  

    n_e = x_e_vec * nc.A_H * n_vec
    n_H = (1. - x_e_vec) * nc.A_H * n_vec
    
    num = B_coeff * (I_vec + I_CMB) + A_coeff + n_e * C_ul_e + n_H * C_ul_H
    
    den =  B_coeff * (I_vec + I_CMB) + n_e * C_ul_e * np.exp(-T_star/T_vec) + n_H * C_ul_H * np.exp(-T_star/T_vec)

    T_spin = T_star / np.log(num/den)

    return T_spin


def eta_func(T_spin, z):
       
    line_CII_rest_frame = 1900*1e9 

    return 1. - (np.exp(nc.hh*line_CII_rest_frame / (nc.kk*(1.+z)*T_spin)) - 1.) /\
           (np.exp((nc.hh*line_CII_rest_frame) / (nc.kk*(1.+z)**2*nc.CMB_temperature)) - 1.)



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
    
    redshift_CII = 5.

    
    fig,ax= plt.subplots()
    
    n_temp = np.logspace(-5,3, 10000)
    
    T_spin_temp = T_spin_func(n=n_temp, T=1e3, I=1e-20, x_e= 0.5, z=5)
    
    eta_temp = eta_func(T_spin_temp, z=5)
    
    ax.plot(n_temp, eta_temp)
    ax.set_xscale("log")
    
    
    fig,ax= plt.subplots()
    
    
    for z in [3,4,5,6]:
        
        redshift_CII_temp = z
        T_spin_list = np.linspace(10., 105., 1000)
    
        eta_list = 1. - (np.exp((nc.hh*nc.line_CII_rest_frame) / (nc.kk*T_spin_list*(1.+redshift_CII_temp))) - 1.) /  \
                    (np.exp((nc.hh*nc.line_CII_rest_frame) / (nc.kk*(1.+redshift_CII_temp)**2*nc.CMB_temperature)) - 1.)
    
        ax.plot(T_spin_list, eta_list, label="z={}".format(z))
        ax.set_ylim(0.,1.)
        ax.set_xlabel("spin temperature [K]")
        ax.set_ylabel("eta")
        ax.axvline(nc.CMB_temperature*(1+redshift_CII_temp), linestyle='--', color='gray')
        ax.legend()
    
    
    T_temp = np.linspace(10., 10000, 10000)
    
    fig,ax= plt.subplots()
    
    ax.plot(T_temp, np.exp(-nc.T_star/T_temp))
    ax.set_xlabel("temperature [K]")
    ax.set_ylabel("exp(-T_star/T)")
    #ax.set_xlim(10,15)
    
    
    fig,ax= plt.subplots()
    
    ax.plot(T_temp, 8.63e-6 / np.sqrt(T_temp))
    ax.set_xlabel("temperature [K]")
    ax.set_ylabel("C_e(T)")
    #ax.set_xlim(10,15)
    
    
    
    fig,ax= plt.subplots()
    
    n_temp = np.logspace(-5,3, 10000)
    alpha_temp = 1.0
    beta_temp = 1.0
    
    T_spin_temp = nc.T_star /  np.log((1.+ 1e-3 +n_temp*alpha_temp)/(1e-3 + n_temp*alpha_temp*beta_temp))
    
    eta_temp = 1. - (np.exp((nc.hh*nc.line_CII_rest_frame) / (nc.kk*T_spin_temp*(1.+redshift_CII))) - 1.) /  \
                    (np.exp((nc.hh*nc.line_CII_rest_frame) / (nc.kk*(1.+redshift_CII)**2*nc.CMB_temperature)) - 1.)
    
    ax.plot(n_temp, eta_temp)
    ax.set_xlabel("n [cm-3]")
    ax.set_ylabel("eta")
    ax.set_xscale("log")
    
    
    #print(eta)











