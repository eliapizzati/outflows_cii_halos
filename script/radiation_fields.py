# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 14:02:23 2021

@author: anna
"""



import numpy as np

import os

import natconst as nc
import mydir


#photoionization params (THE ORDER IS H, He, CI, CII)

nu_thr = np.asarray([1.097,1.983,0.909,1.97])*1e5 #cm^-1

threshold = nu_thr**(-1) / nc.angs

alfa_T = np.asarray([6.3,7.83,12.2,4.60])*1e-18 #cm^2

a = [2.99,2.05,2.0,3.0]

b = [1.34,1.66,3.32,1.95]


#UV BACKGROUND
    
data_flux_UVB =  np.loadtxt(os.path.join(mydir.script_dir, "input_data", "krome_HMflux.dat"), unpack=True)
    
redshift, energy, UVB_flux = data_flux_UVB
    
energy_exp = []
flux_exp_UVB = []
    
for i in range(len(redshift)):
    if redshift[i] == 5.940:
            energy_exp.append(energy[i])
            flux_exp_UVB.append(UVB_flux[i])
    
energy_exp = np.asarray(energy_exp)        
flux_exp_UVB = np.asarray(flux_exp_UVB)
    
energy_eV = 10**energy_exp #eV
    
energy = energy_eV * nc.ev #erg
    
wav = nc.hh * nc.cc / energy / nc.angs #Angstrom
    
    
flux_eV_UVB = 10**flux_exp_UVB #eV/s/cm2/Hz/sr

flux_UVB = flux_eV_UVB * nc.ev # erg/s/cm2/Hz/sr
    
# photoionization coefficients
    
# CALCULATIONS
    
photoionization_rates_UVB = []
    
for i in range(0,4):
        
        
    wav_thr = wav[wav<=threshold[i]]

    flux_thr_UVB = flux_UVB[wav<=threshold[i]]

    alfa_crosssec = alfa_T[i] * ( b[i] * (threshold[i]/wav_thr)**(-a[i]) + (1-b[i]) * (threshold[i]/wav_thr)**(-1.-a[i]))
    
    N_tot = np.trapz(4*np.pi*flux_thr_UVB*alfa_crosssec/ (nc.hh*wav_thr), - wav_thr)
    
    photoionization_rates_UVB.append(N_tot)  


Gamma_H_UVB, Gamma_He_UVB, Gamma_CI_UVB, Gamma_CII_UVB = photoionization_rates_UVB


# LW photodestruction coefficient
    
index = np.searchsorted(energy_eV, 12.87)
    
Gamma_LW_UVB= 1.38e9 * flux_UVB[index] #s^-1 --> value ca 2.05e-13
    

print("######### UVB ###########")
print("Gamma_H_UVB (s^-1) =", Gamma_H_UVB) 
print("Gamma_He_UVB (s^-1) =", Gamma_He_UVB) 
print("Gamma_CI_UVB (s^-1) =", Gamma_CI_UVB) 
print("Gamma_CII_UVB (s^-1) =", Gamma_CII_UVB) 
print("Gamma_LW_UVB (s^-1) =", Gamma_LW_UVB)     
print("#########     ###########")


    
# GALACTIC FLUX
    
# integral of the spectral luminosity
    
    
data_flux_gal = np.loadtxt(os.path.join(mydir.script_dir, "input_data", "flux.dat"), unpack=True, usecols = (0,36))
    
wav = data_flux_gal[0]  #Armstrong
    
exp_lum = data_flux_gal[1]
    
lum_wav = 10**(exp_lum) #erg/(s*A); missing SFR and f_esc, will add them later

# CALCULATIONS

R_sample = 1000. * nc.pc

photoionization_rates_gal = []
    
for i in range(0,4):
    
    wav_thr = wav[wav<=threshold[i]]

    lum_wav_thr = lum_wav[wav<=threshold[i]]

    alfa_crosssec = alfa_T[i] * ( b[i] * (threshold[i]/wav_thr)**(-a[i]) + (1-b[i]) * (threshold[i]/wav_thr)**(-1.-a[i]))

    N_tot = np.trapz(lum_wav_thr*(wav_thr*nc.angs)*alfa_crosssec/ (nc.hh*nc.cc), wav_thr)
            
    photoionization_rates_gal.append(N_tot / (4*np.pi*R_sample**2))  


Gamma_H_1000, Gamma_He_1000, Gamma_CI_1000, Gamma_CII_1000 = photoionization_rates_gal
   
# LW photodestruction coefficient

flux_LW = 10**40.7

Gamma_LW_1000 = 1.38e9 * flux_LW * (963*nc.angs)**2 / (4**2*np.pi**2*(R_sample)**2*nc.cc*nc.angs) 



print("######### gal ###########")
print("Gamma_H_1000 (s^-1) =", Gamma_H_1000) 
print("Gamma_He_1000 (s^-1) =", Gamma_He_1000) 
print("Gamma_CI_1000 (s^-1) =", Gamma_CI_1000) 
print("Gamma_CII_1000 (s^-1) =", Gamma_CII_1000) 
print("Gamma_LW_1000 (s^-1) =", Gamma_LW_1000)     
print("#########     ###########")
    