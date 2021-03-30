# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 15:48:30 2021

@author: anna
"""


import os
import numpy as np
import matplotlib.pyplot as plt 

import natconst as nc
from load_data import params_obs
from post_sol_modules import get_ionization_states
from model_classes import load_from_file
import mydir



# loading profiles for varibales & ionization states

params = dict([("class", "sol"),
               ("type", "n"),
               ("SFR", 20.),
               ("beta", 1.0), 
               ("f_esc", 0.), 
               ("v_c", 175.),
               ("Zeta", 1.0),
               ("alfa", 1.0),
               ("R_in", 0.3)])    


profiles = load_from_file(params, type_class = "sol_profiles")

ionization_state = get_ionization_states(profiles, params)

# loading params from obs and radation fields

line_CII_rest_frame = params_obs["line_frequency"] # in Hz
redshift_CII = params_obs["redshift"]

line_CII_redshifted = line_CII_rest_frame / (1.+ redshift_CII) # in Hz

line_CII_wav_rest_frame =  nc.cc / line_CII_rest_frame / nc.angs  # in Angs

line_CII_wav_redshifted =  nc.cc / line_CII_redshifted / nc.angs  # in Angs

data_flux_UVB =  np.loadtxt(os.path.join(mydir.script_dir, "input_data", "krome_HMflux.dat"), unpack=True)
data_flux_gal =  np.loadtxt(os.path.join(mydir.script_dir, "input_data", "flux.dat"), unpack=True, usecols = (0,36))

# UVB

redshift_list, energy_exp, flux_exp_UVB = data_flux_UVB
    
energy_eV = 10**energy_exp #eV
    
energy = energy_eV * nc.ev #erg
    
wav_UVB = nc.hh * nc.cc / energy / nc.angs #Angstrom

flux_eV_UVB = 10**flux_exp_UVB #eV/s/cm2/Hz/sr

flux_UVB = flux_eV_UVB * nc.ev # erg/s/cm2/Hz/sr

index_redshift = np.searchsorted(redshift_list, redshift_CII) 
redshift_value = redshift_list[index_redshift-1]

wav_z = []
flux_UVB_z = []

for i in range(len(redshift_list)):
    if redshift_list[i] == redshift_value:
            wav_z.append(wav_UVB[i])
            flux_UVB_z.append(flux_UVB[i])

wav_z = np.asarray(wav_z)
intensity_UVB_z = np.asarray(flux_UVB_z)

wav_z_rev = np.flip(wav_z)

index_UVB = np.searchsorted(wav_z, line_CII_wav_rest_frame)

intensity_UVB_z_rev = np.flip(intensity_UVB_z)

intensity_UVB_CII = intensity_UVB_z_rev[index_UVB-1]



# gal flux

wav_gal, exp_lum_gal = data_flux_gal  # wav in A

lum_wav_gal = 10**(exp_lum_gal) #erg/(s*A); missing SFR and f_esc, will add them later

R_sample = 1000. * nc.pc

intensity_gal = lum_wav_gal * (wav_gal*nc.angs)**2 / (4**2*np.pi**2*(R_sample)**2*nc.cc*nc.angs)

index_gal = np.searchsorted(wav_gal, line_CII_wav_rest_frame)

intensity_gal_CII = intensity_gal[index_gal-1]


print("####################")
print("intensity CII gal at 1 kpc (erg/cm^2/s/Hz/sr) =", intensity_gal_CII) 
print("intensity CII UVB (erg/cm^2/s/Hz/sr) =", intensity_UVB_CII) 
print("####################")

      
intensity_tot = intensity_UVB_CII + intensity_gal_CII * (1000*nc.pc/profiles.r)**2 * params["SFR"]#*params["f_esc"]


# computing the spin temperature for CII


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
    

    
def T_spin_func_vec(n_vec, T_vec, I_vec, x_e_vec, z):
    
    
    line_CII_rest_frame = 1900*1e9 
    
    I_CMB = (2*nc.hh*line_CII_rest_frame**3) /  \
            (nc.cc**2 * (np.exp((nc.hh*line_CII_rest_frame) / (nc.kk*(1.+z)*nc.CMB_temperature)) - 1.))
    
    #I_CMB = 0.
    
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

def eta_func_vec(T_spin_vec, z):
       
    line_CII_rest_frame = 1900*1e9 

    return 1. - (np.exp(nc.hh*line_CII_rest_frame / (nc.kk*(1.+z)*T_spin_vec)) - 1.) /\
           (np.exp((nc.hh*line_CII_rest_frame) / (nc.kk*(1.+z)**2*nc.CMB_temperature)) - 1.)
    
T_spin = T_spin_func_vec(n_vec=profiles.n, T_vec=profiles.T, I_vec= intensity_tot, x_e_vec = ionization_state.x_e, z=redshift_CII)

print(T_spin)

eta = eta_func_vec(T_spin_vec=T_spin, z= redshift_CII)

print(eta)

#fig,ax= plt.subplots()
#
#n_temp = np.logspace(-5,3, 10000)
#
#T_spin_temp = T_spin_func(n=n_temp, T=1e3, I=1e-20, x_e = 0.5, z=5)
#
#eta_temp = eta_func(T_spin_temp, z=5)
#
#ax.plot(n_temp, eta_temp)
#ax.set_xscale("log")

#
#fig,ax= plt.subplots()
#
#
#for z in [3,4,5,6]:
#    
#    redshift_CII_temp = z
#    T_spin_list = np.linspace(10., 105., 1000)
#
#    eta_list = 1. - (np.exp((nc.hh*line_CII_rest_frame) / (nc.kk*T_spin_list*(1.+redshift_CII_temp))) - 1.) /  \
#                (np.exp((nc.hh*line_CII_rest_frame) / (nc.kk*(1.+redshift_CII_temp)**2*nc.CMB_temperature)) - 1.)
#
#    ax.plot(T_spin_list, eta_list, label="z={}".format(z))
#    ax.set_ylim(0.,1.)
#    ax.set_xlabel("spin temperature [K]")
#    ax.set_ylabel("eta")
#    ax.axvline(nc.CMB_temperature*(1+redshift_CII_temp), linestyle='--', color='gray')
#    ax.legend()
#
#
#T_temp = np.linspace(10., 10000, 10000)
#
#fig,ax= plt.subplots()
#
#ax.plot(T_temp, np.exp(-T_star/T_temp))
#ax.set_xlabel("temperature [K]")
#ax.set_ylabel("exp(-T_star/T)")
##ax.set_xlim(10,15)
#
#
#fig,ax= plt.subplots()
#
#ax.plot(T_temp, 8.63e-6 / np.sqrt(T_temp))
#ax.set_xlabel("temperature [K]")
#ax.set_ylabel("C_e(T)")
##ax.set_xlim(10,15)
#
#
#
#fig,ax= plt.subplots()
#
#n_temp = np.logspace(-5,3, 10000)
#alpha_temp = 1.0
#beta_temp = 1.0
#
#T_spin_temp = T_star /  np.log((1.+ 1e-3 +n_temp*alpha_temp)/(1e-3 + n_temp*alpha_temp*beta_temp))
#
#eta_temp = 1. - (np.exp((nc.hh*line_CII_rest_frame) / (nc.kk*T_spin_temp*(1.+redshift_CII))) - 1.) /  \
#                (np.exp((nc.hh*line_CII_rest_frame) / (nc.kk*(1.+redshift_CII)**2*nc.CMB_temperature)) - 1.)
#
#ax.plot(n_temp, eta_temp)
#ax.set_xlabel("n [cm-3]")
#ax.set_ylabel("eta")
#ax.set_xscale("log")


#print(eta)











