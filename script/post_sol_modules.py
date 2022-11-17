"""
This script contains the functions to compute the ionization state of the gas,
the surface density of the CII line, the surface brightness of the CII line (raw and convolved) and
the total CII luminosity.
"""

import os
import natconst as nc
import mydir 


import numpy as np

from model_classes import sol_profiles, ion_profiles, lum_profile

from my_utils import twod_making

from cmb_suppression import T_spin_func_vec, eta_func
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
            
    if "f_esc_ion" in params:
        f_esc_ion = params["f_esc_ion"]
    else:
        f_esc_ion = 0.
        
    if "f_esc_FUV" in params:
        f_esc_FUV = params["f_esc_FUV"]
    else:
        f_esc_FUV = 0.

    
    # unpacking profiles 
    
    v = profiles.v
    n = profiles.n
    T = abs(profiles.T)
    
    #UV BACKGROUND
        
    gamma_H = UVB_rates(redshift, quantity="H rate")
    gamma_CI = UVB_rates(redshift, quantity="CI rate")
    gamma_CII = UVB_rates(redshift, quantity="CII rate")
                
    if f_esc_FUV != 0.0:
    
        gamma_CI += nc.Gamma_CI_FUV_1000 * (1000*nc.pc/profiles.r)**2 * SFR_pure*f_esc_FUV
    
    if f_esc_ion != 0.0:
  
        gamma_H += nc.Gamma_H_1000 * (1000*nc.pc/profiles.r)**2 * SFR_pure*f_esc_ion        
        gamma_CI += nc.Gamma_CI_EUV_1000 * (1000*nc.pc/profiles.r)**2 * SFR_pure*f_esc_ion
        gamma_CII += nc.Gamma_CII_1000 * (1000*nc.pc/profiles.r)**2 * SFR_pure*f_esc_ion
            
        
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
            


def get_surface_density(profiles, ionization_states, params, add_CMB_suppression=True, h_resol = 500, rmax=30.):
    """
    computes the surface density predicted by the model; 

    Parameters
    ==========    
    profiles: sol_profiles class
    
    ionization_states: ion_profiles class

    params: dict
        parameters to be passed to all functions
        
    add_CMB_suppression: Boolean, optional
        if True, the CMB suppression in accounted for
    h_resol: int, optional
        number of points for the surface density profile
    rmax: int, optional
        maximum radius for the surface density profile
    Returns
    =======
    sigma_CII: lum_profile class element

    """    
    
    # params definition
    
    if "redshift" in params:
        redshift = params["redshift"]
    else: 
        raise ValueError("No redshift given")

    if "SFR" in params:
        SFR_pure = params["SFR"]
    else: 
        raise ValueError("No SFR given")
        
    if "redshift" in params:
        redshift = params["redshift"]
    else: 
        raise ValueError("No redshift given")
            
    if "f_esc_ion" in params:
        f_esc_ion = params["f_esc_ion"]
    else:
        f_esc_ion = 0.
        
    if "f_esc_FUV" in params:
        f_esc_FUV = params["f_esc_FUV"]
    else:
        f_esc_FUV = 0.
    
    if "Zeta" in params:
        Zeta = params["Zeta"]
    else: 
        Zeta = 1.0

    
    # unpacking profiles 
    
    v = profiles.v
    n = profiles.n
    T = abs(profiles.T)
    
    x_e = ionization_states.x_e
    x_CII = ionization_states.x_CII
    
    # computing the emissivity 
    
    #n_CII = nc.A_C * Zeta * x_CII * n

    #epsilon = 7.9e-20 * n**2 * (nc.A_C * Zeta) * nc.A_H * x_e * x_CII * T**(-0.5) * np.exp(-92./T)
    epsilon = 7.9e-20 *  n**2 * (nc.A_C * Zeta) * nc.A_H * x_e * x_CII * np.exp(-92./T) /  92.**0.5
    
    if add_CMB_suppression == True:
        
        intensity_tot = UVB_rates(redshift, quantity="UV intensity") +\
                        nc.intensity_UV_1000 * (1000*nc.pc/profiles.r)**2 * SFR_pure*f_esc_FUV
        
        T_spin = T_spin_func_vec(n_vec=profiles.n, T_vec=profiles.T, I_UV_vec= intensity_tot, x_e_vec = ionization_states.x_e, z=redshift)

        eta = eta_func(T_spin=T_spin, z=redshift)
        
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

    if add_CMB_suppression == True:
        
        eta = np.interp(h, profiles.r, eta, right=0.)
        
        return lum_profile(radius=h, variable=sigma_CII, params=params, category="sigma", eta=eta)

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
    
    return lum_profile(radius=sigma_CII.h, variable=intensity, params=params, category="int_raw", eta=sigma_CII.eta)



def get_intensity_convolved(intensity_raw, params, params_obs, obs_data, add_central_contribution=False):
    """
    computes convolved intensity from the raw intensity

    Parameters
    ==========
    intensity_raw: lum_profile class element
    params: dict
    params_obs: dict
    obs_data: obs_data class element
    add_central_contribution: bool
        if True, adds the central contribution of the galaxy to the intensity

    Returns
    =======
    intensity_convolved: lum_profile class element

    """    
    
    if "SFR" in params:
        SFR_pure = params["SFR"]
    else: 
        raise ValueError("No SFR given")

    # creates the 2d profiles
    
    x, profile_2d = twod_making(intensity_raw.var, intensity_raw.h, nimage=1000)
    
    beam_interp = np.interp(intensity_raw.h, obs_data.x_beam, obs_data.beam, right=0.)
    
    beam_interp[beam_interp<0.] = 0.

    x_beam, beam_2d = twod_making(beam_interp, intensity_raw.h, nimage=1000)

#    import matplotlib.pyplot as plt
#    
#    fig,ax = plt.subplots()
#    ax.contourf(x/1e3/nc.pc,x/1e3/nc.pc,profile_2d, alpha=0.5)
#    ax.contourf(x_beam/1e3/nc.pc,x_beam/1e3/nc.pc,beam_2d, alpha=0.2)
#    
#    print("x", x[-1]/1e3/nc.pc)
#    print("x beam", x[-1]/1e3/nc.pc)
    
    #makes the 2dfft
    
    f_beam = np.fft.fft2(beam_2d)
    f_imag = np.fft.fft2(profile_2d)
                        
    #convolution
    f_cimag = f_beam * f_imag
    cimage = np.real(np.fft.ifftshift(np.fft.ifft2(f_cimag)))
    cprofile_2d = cimage[:, cimage.shape[1]//2]
    
    intensity_convolved = np.interp(intensity_raw.h, x, cprofile_2d, right=0.)
    
    #normalizes the convolved intensity
    
    norm_intensity = np.trapz(2*np.pi * intensity_raw.h * intensity_raw.var, intensity_raw.h)\
                                / np.trapz(2*np.pi * intensity_raw.h * intensity_convolved, intensity_raw.h)
    
    intensity_convolved *= norm_intensity
    
    if add_central_contribution == True:
        
        beam_interp = np.interp(intensity_raw.h, obs_data.x_beam, obs_data.beam)
        
        luminosity_central = nc.ls * 1e7 * SFR_pure

        factor = luminosity_central / np.trapz(2*np.pi*intensity_raw.h*beam_interp, intensity_raw.h)

        intensity_convolved += factor * beam_interp


    return lum_profile(radius=intensity_raw.h, variable=intensity_convolved, params=params, category="int_conv", eta=intensity_raw.eta)



def get_chi2(intensity_convolved, obs_data):
    """
    computes the chi squared given the data and the intensity profiles
    
    Parameters
    ==========
    intensity_convolved: lum_profile class element
    obs_data: obs_data class element

    Returns
    =======
    chi2: float

    """    
    
    from scipy.interpolate import interp1d
    
    emission_profile = interp1d(intensity_convolved.h, intensity_convolved.var)
    

    res = emission_profile(obs_data.x) - obs_data.data
    
    residuals = 2*res / (obs_data.err_down + obs_data.err_up) 
    
    chi2 = np.sum(residuals**2)
    
    return chi2

    


def get_luminosity_CII(sigma_CII, r_min = None, r_max = None):
    """
    computes the total CII luminosity between two radii (from the CII sutface density)
    
    Parameters
    ==========
    sigma_CII: lum_profile class element
    
    r_min: float, optional
        minimum radius in cm (if None, the integral is computed for all elements in sigma_CII)    
    
    r_max: float, optional
        maximum radius in cm (if None, the integral is computed for all elements in sigma_CII)
    
    Returns
    =======
    lum_CII: float representing the luminosity of CII in erg/s

    """    

    if r_max != None and r_min != None:
        mask = np.logical_and(sigma_CII.h > r_min, sigma_CII.h < r_max)
    else:
        mask = True
    
    lum_CII = np.trapz(sigma_CII.var[mask]*2*np.pi*sigma_CII.h[mask], sigma_CII.h[mask])
    
    return lum_CII
    

    
    
    
    
def get_halo_mass(profiles, params, r_min = None, r_max = None):
    """
    computes the total Carbon (and CII) mass enclosed between two radii 
    
    Parameters
    ==========
    sol_profiles: sol_profiles class element
    
    r_min: float, optional
        minimum radius in cm (if None, the integral is computed for all elements in profiles)    
    
    r_max: float, optional
        maximum radius in cm (if None, the integral is computed for all elements in profiles)
    
    Returns
    =======
    C_mass, CII_mass: floats representing the total Carbon mass and the total CII mass (in solar masses)

    """ 

    if r_max != None and r_min != None:
        mask = np.logical_and(profiles.r > r_min, profiles.r < r_max)
    else:
        mask = True

    rho = profiles.n*nc.mp*nc.mus

    M = np.trapz(4*np.pi*profiles.r[mask]**2*rho[mask]*nc.A_C, profiles.r[mask])

    M_C = np.trapz(4*np.pi*profiles.r[mask]**2*rho[mask]*nc.A_C, profiles.r[mask])

    ionizations = get_ionization_states(profiles, params)
    
    M_CII = np.trapz(4*np.pi*profiles.r[mask]**2*rho[mask]*nc.A_C*ionizations.x_CII[mask], profiles.r[mask])
    
    return M_C/nc.ms/1e6, M_CII/nc.ms/1e6
        
    

