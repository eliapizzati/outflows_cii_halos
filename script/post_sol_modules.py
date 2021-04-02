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
from cmb_suppression import T_spin_func_vec, eta_func_vec
from radiation_fields import UVB_rates



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
        
    if "redshift" in params:
        redshift = params["redshift"]
    else: 
        raise ValueError("No redshift given")
            
    if "f_esc" in params:
        f_esc = params["f_esc"]
    else:
        f_esc = 0.
    
    # unpacking profiles 
    
    v = profiles.v
    n = profiles.n
    T = profiles.T
    
    #UV BACKGROUND
        
    gamma_H = UVB_rates(redshift, quantity="H rate")
    gamma_CI = UVB_rates(redshift, quantity="CI rate")
    gamma_CII = UVB_rates(redshift, quantity="CII rate")
                
    if f_esc != 0.0:
    
        gamma_H += nc.Gamma_H_1000 * (1000*nc.pc/profiles.r)**2 * SFR_pure*f_esc              
        gamma_CI += nc.Gamma_CI_1000 * (1000*nc.pc/profiles.r)**2 * SFR_pure*f_esc
        gamma_CII += nc.Gamma_CII_1000 * (1000*nc.pc/profiles.r)**2 * SFR_pure*f_esc
            
        
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
            


def get_surface_density(profiles, ionization_states, params, add_CMB_suppression=False, h_resol = 1000, rmax=10.):
    """
    computes the surface density predicted by the model; 
    
    possibility to add the central contribution from the galaxy
    
    Parameters
    ==========    
    profiles: sol_profiles class
    
    ionization_states: ion_profiles class
    
    params: dict
        parameters to be passed to all functions
        
    add_CMB_suppression: Boolean, optional
    
    h_resol: int, optional
    
    rmax: int, optional
    
    Returns
    =======
    sigma_CII: lum_profile class element

    """    
    
    # params definition
    
    if "redshift" in params:
        redshift = params["redshift"]
    else: 
        raise ValueError("No redshift given")

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
    
    if add_CMB_suppression == True:
        
        intensity_tot = UVB_rates(redshift, quantity="CII intensity") +\
                        nc.intensity_CII_1000 * (1000*nc.pc/profiles.r)**2 * params["SFR"]*params["f_esc"]
        
        T_spin = T_spin_func_vec(n_vec=profiles.n, T_vec=profiles.T, I_vec= intensity_tot, x_e_vec = ionization_states.x_e, z=redshift)

        eta = eta_func_vec(T_spin_vec=T_spin, z=redshift)
        
        epsilon *= eta

    # intergrating along the line of sight
    
    h = np.linspace(min(profiles.r),max(max(profiles.r),rmax*1000*nc.pc), h_resol)
   
    sigma_CII = np.zeros_like(h)
        
#    sigma_CII[:] = 2* np.trapz(epsilon[profiles.r>h[:]] * profiles.r[profiles.r>h[:]] / np.sqrt((profiles.r[profiles.r>h[:]])**2 - h[:]**2),\
#                 profiles.r[profiles.r>h[:]])        

    sigma_CII = []
    for el_h in h:
        
        integral = np.trapz(epsilon[profiles.r>el_h] * profiles.r[profiles.r>el_h] / np.sqrt((profiles.r[profiles.r>el_h])**2 - el_h**2), profiles.r[profiles.r>el_h])        
        sigma_CII.append(2.*integral)
                
    sigma_CII = np.asarray(sigma_CII)


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
    
    if "redshift" in params_obs:
        redshift = params_obs["redshift"]
    else: 
        raise ValueError("No redshift given")

    if "line_FWHM" in params_obs:
        FWHM_vel = params_obs["line_FWHM"]
    else: 
        raise ValueError("No line_FWHM given")

    nu_0 = nc.line_CII_rest_frame
    
    # transforming sigma to the intensity
    
    intensity = sigma_CII.var / (nu_0*4*np.pi*(1.+redshift)**3*FWHM_vel/nc.cc)
        
    # changing the units
            
    intensity *= 1e26 #transformation to mJy 
    intensity /= 4.2e10 #transforming sr to arcsec^2
    
    return lum_profile(radius=sigma_CII.h, variable=intensity, params=params, category="int_raw")



def get_intensity_convolved(intensity_raw, params, params_obs, obs_data, add_central_contribution=False):
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
    
    add_CMB_suppression: Boolean, optional
    
    Returns
    =======
    intensity_convolved: lum_profile class element

    """    
    
    if "SFR" in params:
        SFR_pure = params["SFR"]
    else: 
        raise ValueError("No SFR given")

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
    
    if add_central_contribution == True:
        
        beam_interp = np.interp(intensity_raw.h, obs_data.x_beam, obs_data.beam)
        
        luminosity_central = nc.ls * 1e7 * SFR_pure

        factor = luminosity_central / np.trapz(2*np.pi*intensity_raw.h*beam_interp, intensity_raw.h)

        intensity_convolved += factor * beam_interp


    return lum_profile(radius=intensity_raw.h, variable=intensity_convolved, params=params, category="int_conv")



def get_chi2(intensity_convolved, obs_data):
    """
    computes the chi squared given the data and the intensity profiles
    
    Parameters
    ==========
    r: array
        
    intensity_raw: array
    
    params: dict
    
    params_obs: dict
        parameters from observations
        
    obs_data: obs_data class element
    
    add_CMB_suppression: Boolean, optional
    
    Returns
    =======
    intensity_convolved: lum_profile class element

    """    
    
    from scipy.interpolate import interp1d
    
    emission_profile = interp1d(obs_data.x, intensity_convolved)

    res = emission_profile(obs_data.x) - obs_data.data
    
    residuals = 2*res / (obs_data.err_down + obs_data.err_up) 
    
    chi2 = np.sum(residuals**2)
    
    return chi2

    
    

    
#    #CALCULATION OF THE CARBON AND IONIZED CARBON HALO MASSES
#    
#    p = n*nc.mp*mu
#
#    M_C = np.trapz(4*np.pi*r**2*p*A_C, r)
#
#    M_CII = np.trapz(4*np.pi*r**2*p*A_C*x_CII, r)
#    print('M_C  ',beta, M_C/nc.ms/1e6)
#    print('M_CII',beta, M_CII/nc.ms/1e6)
        
    

