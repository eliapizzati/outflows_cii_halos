# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 14:02:23 2021

@author: anna
"""



import numpy as np

import os

import natconst as nc
import mydir
from my_utils import to_habing_flux 


def UVB_rates(redshift, quantity):
    """
    computes many quantities loading the UVB flux (use the Haardt&Madau2012 UVB)
    
    Parameters
    ==========
    redshift: float
    
    quantity: string
        quantities can be "UV intensity", "H rate", "He rate", "CII rate", "CI rate"

    
    Returns
    =======
    UVB rate selected

    """  

    #photoionization params (THE ORDER IS H, He, CI, CII)
    
    nu_thr = np.asarray([1.097,1.983,0.909,1.97, 3.81])*1e5 # in cm^-1
    
    threshold = nu_thr**(-1) # in cm
    
    alfa_T = np.asarray([6.3,7.83,12.2,4.60,1.0])*1e-18 # in cm^2
    
    a = [2.99,2.05,2.0,3.0,1.0]
    
    b = [1.34,1.66,3.32,1.95,1.0]
    
    
    #UV BACKGROUND
        
    data_flux_UVB =  np.loadtxt(os.path.join(mydir.script_dir, "input_data", "krome_HMflux.dat"), unpack=True)
        
    redshift_list, energy_exp, flux_exp_UVB = data_flux_UVB
        
    energy_eV = 10**energy_exp # in eV
        
    energy = energy_eV * nc.ev # in erg
        
    wav_UVB = nc.hh * nc.cc / energy  # in cm 
    
    flux_eV_UVB = 10**flux_exp_UVB # in eV/s/cm2/Hz/sr
    
    flux_UVB = flux_eV_UVB * nc.ev # in erg/s/cm2/Hz/sr
    
    
    index_redshift = np.searchsorted(redshift_list, redshift) 
    redshift_value = redshift_list[index_redshift-1]
        
    energy_eV_z =[]
    wav_z = []
    flux_UVB_z = []
        
    for i in range(len(redshift_list)):
        if redshift_list[i] == redshift_value:
                energy_eV_z.append(energy_eV[i])
                wav_z.append(wav_UVB[i])
                flux_UVB_z.append(flux_UVB[i])
    
    energy_eV_z = np.asarray(energy_eV_z)
    wav_z = np.asarray(wav_z)
    flux_UVB_z = np.asarray(flux_UVB_z)
    

    
    if quantity == "UV intensity":
        
        wav_z_rev = np.flip(wav_z)
  
        line_UV_wav_rest_frame =  nc.cc / nc.line_UV_mixing_rest_frame  # in cm
    
        index_UVB = np.searchsorted(wav_z_rev, line_UV_wav_rest_frame)
    
        intensity_UVB_z_rev = np.flip(flux_UVB_z)
        
        intensity_UVB_UV = intensity_UVB_z_rev[index_UVB-1]

        return intensity_UVB_UV   
    
        
    elif quantity == "H rate" or quantity == "He rate" or quantity == "CII rate" or quantity == "CI rate":
        
        # photoionization coefficients
            
        # CALCULATIONS
            
        photoionization_rates_UVB = []
            
        for i in range(0,5):
                
                
            wav_thr = wav_z[wav_z<=threshold[i]]
        
            flux_thr_UVB = flux_UVB_z[wav_z<=threshold[i]]
        
            alfa_crosssec = alfa_T[i] * ( b[i] * (threshold[i]/wav_thr)**(-a[i]) + (1-b[i]) * (threshold[i]/wav_thr)**(-1.-a[i]))
            
            N_tot = np.trapz(4*np.pi*flux_thr_UVB*alfa_crosssec/ (nc.hh*wav_thr), - wav_thr)
            
            photoionization_rates_UVB.append(N_tot)  
        
        
        Gamma_H_UVB, Gamma_He_UVB, Gamma_CI_UVB, Gamma_CII_UVB, Gamma_CVI_UVB = photoionization_rates_UVB
        
        if quantity == "H rate":
            return Gamma_H_UVB
        elif quantity == "He rate":
            return Gamma_He_UVB
        elif quantity == "CII rate":
            return Gamma_CII_UVB
        elif quantity == "CI rate":
            return Gamma_CI_UVB
    
    
    elif quantity == "LW rate":
        
        # LW photodestruction coefficient
            
        index = np.searchsorted(energy_eV_z, 12.87)
            
        Gamma_LW_UVB= 1.38e9 * flux_UVB_z[index] #s^-1 --> value ca 2.05e-13
        
        return Gamma_LW_UVB
    
    else:
        raise ValueError("No correct quantity in input (quantities are UV intensity, H rate, He rate, CII rate, CI rate, LW rate)")
    
        

if __name__ == "__main__":
    
    redshift = 5
    
    # photoionization params (THE ORDER IS H, He, CI, CII, CVI?)
    
    nu_thr = np.asarray([1.097,1.983,0.909,1.97, 3.81])*1e5 # in cm^-1
    
    threshold = nu_thr**(-1) # in cm
        
    alfa_T = np.asarray([6.3,7.83,12.2,4.60,1.0])*1e-18 # in cm^2
    
    a = [2.99,2.05,2.0,3.0, 1.0]
    
    b = [1.34,1.66,3.32,1.95, 1.0]
    
    
    # UV BACKGROUND
        
    data_flux_UVB =  np.loadtxt(os.path.join(mydir.script_dir, "input_data", "krome_HMflux.dat"), unpack=True)
        
    redshift_list, energy_exp, flux_exp_UVB = data_flux_UVB
        
    energy_eV = 10**energy_exp # eV
        
    energy = energy_eV * nc.ev # erg
        
    wav_UVB = nc.hh * nc.cc / energy # in cm
    
    flux_eV_UVB = 10**flux_exp_UVB # eV/s/cm2/Hz/sr
    
    flux_UVB = flux_eV_UVB * nc.ev # erg/s/cm2/Hz/sr
    
    
    index_redshift = np.searchsorted(redshift_list, redshift) 
    redshift_value = redshift_list[index_redshift-1]
        
    energy_eV_z =[]
    wav_z = []
    flux_UVB_z = []
        
    for i in range(len(redshift_list)):
        if redshift_list[i] == redshift_value:
                energy_eV_z.append(energy_eV[i])
                wav_z.append(wav_UVB[i])
                flux_UVB_z.append(flux_UVB[i])
    
    energy_eV_z = np.asarray(energy_eV_z)
    wav_z = np.asarray(wav_z)
    flux_UVB_z = np.asarray(flux_UVB_z)
    
    
    wav_z_rev = np.flip(wav_z)
    
    line_UV_wav_rest_frame =  nc.cc / nc.line_UV_mixing_rest_frame # in cm

    index_UVB = np.searchsorted(wav_z_rev, line_UV_wav_rest_frame)

    intensity_UVB_z_rev = np.flip(flux_UVB_z)
    
    intensity_UVB_UV = intensity_UVB_z_rev[index_UVB-1]

    
        
        
    # photoionization coefficients
        
    # CALCULATIONS
        
    photoionization_rates_UVB = []
        
    for i in range(0,5):
                 
        wav_thr = wav_z[wav_z<=threshold[i]]
    
        flux_thr_UVB = flux_UVB_z[wav_z<=threshold[i]]
        
        alfa_crosssec = alfa_T[i] * ( b[i] * (threshold[i]/wav_thr)**(-a[i]) + (1-b[i]) * (threshold[i]/wav_thr)**(-1.-a[i]))
        
        N_tot = np.trapz(4*np.pi*flux_thr_UVB*alfa_crosssec/ (nc.hh*wav_thr), - wav_thr)
        
        photoionization_rates_UVB.append(N_tot)  
    
    
    Gamma_H_UVB, Gamma_He_UVB, Gamma_CI_UVB, Gamma_CII_UVB, Gamma_CVI_UVB = photoionization_rates_UVB
    
 
    # LW photodestruction coefficient
        
    index = np.searchsorted(energy_eV_z, 12.87)

    intensity_LW_UVB = flux_UVB_z[index]
    Gamma_LW_UVB = 1.38e9 * intensity_LW_UVB  #s^-1 --> value ca 2.05e-13
    

    print("######### UVB ###########")
    print("Gamma_H_UVB (s^-1) =", Gamma_H_UVB) 
    print("Gamma_He_UVB (s^-1) =", Gamma_He_UVB) 
    print("Gamma_CI_UVB (s^-1) =", Gamma_CI_UVB) 
    print("Gamma_CII_UVB (s^-1) =", Gamma_CII_UVB) 
    print("Gamma_CVI_UVB (s^-1) =", Gamma_CVI_UVB) 
    print("Gamma_LW_UVB (s^-1) =", Gamma_LW_UVB)  
    print("intensity_UV_UVB (erg/cm^2/s/Hz/sr) =", intensity_UVB_UV) 
    print("intensity_UV_UVB (G_0) =", to_habing_flux(intensity_UVB_UV))          
    print("intensity_LW_UVB (erg/cm^2/s/Hz/sr) =", intensity_LW_UVB)
    print("intensity_LW_UVB (G_0) =", to_habing_flux(intensity_LW_UVB))
    print("#########     ###########")
    
    
        
    # GALACTIC FLUX
        
    # integral of the spectral luminosity
        
        
    data_flux_gal = np.loadtxt(os.path.join(mydir.script_dir, "input_data", "flux.dat"), unpack=True, usecols = (0,36))
        
    wav = data_flux_gal[0] * nc.angs  # in cm
        
    exp_lum = data_flux_gal[1]
        
    lum_wav = 10**(exp_lum) / nc.angs # in erg/s/cm; missing SFR and f_esc, will add them later
    
    # CALCULATIONS
    
    R_sample = 1000. * nc.pc
    
    photoionization_rates_gal = []
        
    for i in range(0,5):
        
        wav_thr = wav[wav<=threshold[i]]
    
        lum_wav_thr = lum_wav[wav<=threshold[i]]
    
        alfa_crosssec = alfa_T[i] * ( b[i] * (threshold[i]/wav_thr)**(-a[i]) + (1-b[i]) * (threshold[i]/wav_thr)**(-1.-a[i]))
    
        N_tot = np.trapz(lum_wav_thr* wav_thr * alfa_crosssec/ (nc.hh*nc.cc), wav_thr)
                
        photoionization_rates_gal.append(N_tot / (4*np.pi*R_sample**2))  
    
    
    Gamma_H_1000, Gamma_He_1000, Gamma_CI_1000, Gamma_CII_1000, Gamma_CVI_1000 = photoionization_rates_gal
       
    threshold_H = threshold[0]
    threshold_CI = threshold[2]
    
    mask_FUV = np.logical_and(wav <= threshold_CI, wav >= threshold_H)
    mask_EUV = wav <= threshold_H
    
    
    wav_FUV = wav[mask_FUV]
    lum_wav_FUV = lum_wav[mask_FUV]
    alfa_crosssec_FUV = alfa_T[2] * ( b[2] * (threshold[2]/wav_FUV)**(-a[2]) + (1-b[2]) * (threshold[2]/wav_FUV)**(-1.-a[2]))
    N_tot_FUV = np.trapz(lum_wav_FUV* wav_FUV * alfa_crosssec_FUV/ (nc.hh*nc.cc), wav_FUV)


    Gamma_CI_FUV_1000 = N_tot_FUV / (4*np.pi*R_sample**2)
    
    wav_EUV = wav[mask_EUV]
    lum_wav_EUV = lum_wav[mask_EUV]
    alfa_crosssec_EUV = alfa_T[2] * ( b[2] * (threshold[2]/wav_EUV)**(-a[2]) + (1-b[2]) * (threshold[2]/wav_EUV)**(-1.-a[2]))
    N_tot_EUV = np.trapz(lum_wav_EUV* wav_EUV * alfa_crosssec_EUV/ (nc.hh*nc.cc), wav_EUV)

    Gamma_CI_EUV_1000 = N_tot_EUV / (4*np.pi*R_sample**2)
    
    
    # LW photodestruction coefficient
    
    flux_LW = 10**40.7 / nc.angs # in erg/s/cm
    
    Gamma_LW_1000 = 1.38e9 * flux_LW * (963*nc.angs)**2 / (4**2*np.pi**2*(R_sample)**2*nc.cc)  # in s^-1

    # UV intensity (for the CMB suppression)
    
    intensity_gal = lum_wav * wav**2 / (4**2*np.pi**2*(R_sample)**2*nc.cc) # in erg/s/cm^2/Hz/sr
        
    line_UV_wav_rest_frame =  nc.cc / nc.line_UV_mixing_rest_frame  # in cm
    
    index_gal = np.searchsorted(wav, line_UV_wav_rest_frame)

    intensity_UV_1000 = intensity_gal[index_gal-1]
    
    
    print("######### gal ###########")
    print("Gamma_H_1000 (s^-1) =", Gamma_H_1000) 
    print("Gamma_He_1000 (s^-1) =", Gamma_He_1000) 
    print("Gamma_CI_1000 (s^-1) =", Gamma_CI_1000) 
    print("Gamma_CI_FUV_1000 (s^-1) =", Gamma_CI_FUV_1000) 
    print("Gamma_CI_EUV_1000 (s^-1) =", Gamma_CI_EUV_1000) 
    print("Gamma_CII_1000 (s^-1) =", Gamma_CII_1000) 
    print("Gamma_CVI_1000 (s^-1) =", Gamma_CVI_1000) 
    print("Gamma_LW_1000 (s^-1) =", Gamma_LW_1000)  
    print("intensity_UV_1000 (erg/s/cm^2/Hz/sr) =", intensity_UV_1000)   
    print("intensity_UV_1000 (G_0) =", to_habing_flux(intensity_UV_1000))   
    print("#########     ###########")
          
          
