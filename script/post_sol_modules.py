# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 15:10:15 2021

@author: anna
"""

import os
import natconst as nc
import mydir 


import numpy as np

from model_classes import sol_profiles, ion_profiles, lum_profile
from utils import twod_making



def get_ionization_states(profiles, params):
    """
    computes the ionization state for hydrogen (x_HI) and Carbon (x_CII)
    
    Parameters
    ==========    
    profiles: sol_profiles class
    
    params: dict
        parameters to be passed to all functions
    
    Returns
    =======
    ionization_state: ion_profiles class element

    """    
    
    # params definition
    
    if "SFR" in params:
        SFR_pure = params["SFR"]
    else: 
        raise ValueError("No SFR given")
            
    if "f_esc" in params:
        f_esc = params["f_esc"]
    else:
        f_esc = 0.
    
    # unpacking profiles 
    
    v = profiles.v
    n = profiles.n
    T = profiles.T
    
    #UV BACKGROUND
        
    gamma_H = nc.Gamma_H_UVB
    gamma_CI = nc.Gamma_CI_UVB               
    gamma_CII = nc.Gamma_CII_UVB
                
    if f_esc != 0.0:
    
        gamma_H += nc.Gamma_H_1000 * (1000*nc.pc/profiles.r)**2 * SFR_pure*f_esc              
        gamma_CI = nc.Gamma_CI_1000 * (1000*nc.pc/profiles.r)**2 * SFR_pure*f_esc
        gamma_CII = nc.Gamma_CII_1000 * (1000*nc.pc/profiles.r)**2 * SFR_pure*f_esc
            
        
    beta_H = 4.18e-13 * (T/1e4)**(-0.75) #cm^3 s^-1
    beta_CII = 4.66e-13 * (T/1e4)**(-0.62) + 1.84e-13 #cm^3 s^-1
    beta_CIII = 24.5e-13 * (T/1e4)**(-0.65) + 60.6e-13 #cm^3 s^-1
    
    T_eV = 25.7 * T / 1000. / 298

    kappa_H = np.exp(-32.71396786 + 13.5365560*np.log(T_eV)\
                     - 5.73932875*np.log(T_eV)**2 + 1.56315498*np.log(T_eV)**3 \
                     - 0.28770560*np.log(T_eV)**4 + 3.48255977e-2*np.log(T_eV)**5 \
                     - 2.63197617e-3*np.log(T_eV)**6 + 1.11954395e-4*np.log(T_eV)**7 \
                     - 2.03914985e-6*np.log(T_eV)**8)
    
    #kappa_H = 0.
    kappa_CI = 6.85e-8 * (0.193 + 11.26/T_eV)**(-1) * (11.26/T_eV)**0.25 * np.exp(-11.26/T_eV) 
    kappa_CII = 1.86e-8 * (1. + 24.4/T_eV**0.5)*(0.286 + 24.4/T_eV)**(-1) * (24.4/T_eV)**0.24 * np.exp(-24.4/T_eV) 
    
    ratio =  gamma_H / nc.A_H / n
    
    x_e = ( np.sqrt( (ratio-kappa_H)**2 + 4*ratio*(beta_H+kappa_H)) - (ratio-kappa_H) ) / (2*(beta_H+kappa_H))
                
    n_e = x_e * nc.A_H * n
                
    x_CII = 1. / (1. + (gamma_CII)/(beta_CIII * n_e) + (beta_CII * n_e)/(gamma_CI + kappa_CI * n_e) + kappa_CII/beta_CIII)
    
    ionization_states = []
    
    ionization_states.append(x_e)
    ionization_states.append(x_CII)
    
    return ion_profiles(radius=profiles.r, variables=ionization_states, params=params)
            


def get_surface_density(profiles, ionization_states, params, central_contribution=False, h_resol = 1000):
    """
    computes the surface density predicted by the model; 
    
    possibility to add the central contribution from the galaxy
    
    Parameters
    ==========    
    profiles: sol_profiles class
    
    ionization_states: ion_profiles class
    
    params: dict
        parameters to be passed to all functions
        
    central_contribution: Boolean, optional
    
    h_resol: int, optional
    
    Returns
    =======
    sigma_CII: lum_profile class element

    """    
    
    # params definition
    
    if "Zeta" in params:
        Zeta = params["Zeta"]
    else: 
        Zeta = 1.0
    
    # unpacking profiles 
    
    v = profiles.v
    n = profiles.n
    T = profiles.T
    
    x_e = ionization_states.x_e
    x_CII = ionization_states.x_CII
    
    # computing the emissivity 
    
    #n_CII = nc.A_C * Zeta * x_CII * n
      
    epsilon = 3e-27 * n**2 * (nc.A_C * Zeta/1.4e-4) * (1+0.42*x_CII/1e-3) * np.exp(-92./T)
    
    h = np.linspace(min(profiles.r),max(max(profiles.r),10.*1000*nc.pc), h_resol)
   
    sigma_CII = np.zeros_like(h)
        
#    sigma_CII[:] = 2* np.trapz(epsilon[profiles.r>h[:]] * profiles.r[profiles.r>h[:]] / np.sqrt((profiles.r[profiles.r>h[:]])**2 - h[:]**2),\
#                 profiles.r[profiles.r>h[:]])        

    sigma_CII = []
    for el_h in h:
        
        integral = np.trapz(epsilon[profiles.r>el_h] * profiles.r[profiles.r>el_h] / np.sqrt((profiles.r[profiles.r>el_h])**2 - el_h**2), profiles.r[profiles.r>el_h])        
        sigma_CII.append(2.*integral)
                
    sigma_CII = np.asarray(sigma_CII)

    if central_contribution == True:    
        pass

    return lum_profile(radius=h, variable=sigma_CII, params=params, category="sigma")
    
    


def get_intensity_raw(sigma_CII, params, params_obs):
    """
    computes raw intensity (not convolved with the beam) from the surface density
    
    Parameters
    ==========
    sigma_CII: lum_profile class element
    
    params_obs: dict
        parameters from observations

    
    Returns
    =======
    intensity: lum_profile class element

    """    
    
    # params definition
    
    if "line_frequency" in params_obs:
        nu_0 = params_obs["line_frequency"]
    else: 
        raise ValueError("No line_frequency given")

    if "redshift" in params_obs:
        zeta = params_obs["redshift"]
    else: 
        raise ValueError("No redshift given")

    if "line_FWHM" in params_obs:
        FWHM_vel = params_obs["line_FWHM"]
    else: 
        raise ValueError("No line_FWHM given")


    # transforming sigma to the intensity
    
    intensity = sigma_CII.var / (nu_0*4*np.pi*(1.+zeta)**3*FWHM_vel/nc.cc)
        
    # changing the units
            
    intensity *= 1e26 #transformation to mJy 
    intensity /= 4.2e10 #transforming sr to arcsec^2
    
    return lum_profile(radius=sigma_CII.h, variable=intensity, params=params, category="int_raw")



def get_intensity_convolved(intensity_raw, params, params_obs, obs_data):
    """
    computes the ionization state for hydrogen (x_HI) and Carbon (x_CII)
    
    Parameters
    ==========
    r: array
        
    intensity_raw: array
    
    params: dict
    
    params_obs: dict
        parameters from observations
        
    obs_data: obs_data class element
    
    Returns
    =======
    intensity_convolved: lum_profile class element

    """    
    
        
    # creates the 2d profiles
    
    x, y, profile_2d = twod_making(intensity_raw.var, intensity_raw.h, nimage=1000)
    x_beam, y_beam, beam_2d = twod_making(obs_data.beam, obs_data.x_beam, nimage=1000)

    
    #makes the 2dfft
    
    f_beam = np.fft.fft2(beam_2d)
    f_imag = np.fft.fft2(profile_2d)
                        
    #convolution
    f_cimag = f_beam * f_imag
    cimage = np.real(np.fft.ifftshift(np.fft.ifft2(f_cimag)))
    cprofile_2d = cimage[:, cimage.shape[1]//2]
    
    intensity_convolved = np.interp(intensity_raw.h, x, cprofile_2d)
    
    #normalizes the convolved intensity
    
    norm_intensity = np.trapz(2*np.pi * intensity_raw.h * intensity_raw.var, intensity_raw.h)\
                                / np.trapz(2*np.pi * intensity_raw.h * intensity_convolved, intensity_raw.h)
    
    intensity_convolved *= norm_intensity

    return lum_profile(radius=intensity_raw.h, variable=intensity_convolved, params=params, category="int_conv")



def add_central_contribution():
    pass

def add_CMB_suppression():
    pass



    

    
#    #CALCULATION OF THE CARBON AND IONIZED CARBON HALO MASSES
#    
#    p = n*nc.mp*mu
#
#    M_C = np.trapz(4*np.pi*r**2*p*A_C, r)
#
#    M_CII = np.trapz(4*np.pi*r**2*p*A_C*x_CII, r)
#    print('M_C  ',beta, M_C/nc.ms/1e6)
#    print('M_CII',beta, M_CII/nc.ms/1e6)
        
    
    
    # adding to sigma the contribution of the inner galaxy
    
#    luminosity_central = nc.ls * 1e7 * SFR_pure
#    
#    beam_h = fuji.h_data * 1e3 * nc.pc
#    
#    normalization_beam = np.trapz(2*np.pi * beam_h * fuji.beam, beam_h)
#
#    factor = luminosity_central / np.trapz(2*np.pi*h_data*fuji.beam, h_data)
#    
#    
#    sigma_CII_plus_central = sigma_CII + factor * fuji.beam

