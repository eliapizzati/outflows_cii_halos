"""
Created on Fri Mar  5 16:12:54 2021

@author: anna
"""

import os
import natconst as nc
import mydir 
import gnedincooling as gc

import scipy.integrate as si

import numpy as np

from model_classes import sol_profiles, ion_profiles, lum_profile
from utils import twod_making


# FUNCTIONS
        
# cooling function

gc.frtinitcf(0, os.path.join(mydir.script_dir, "input_data", "cf_table.I2.dat"))

def lamda(T, n, r, params):
    """
    Cooling function (calling gnedincooling from Fortran functions)
    
    Parameters
    ==========
    T: float
        temperature
    n: float
        density
    r: float
        radius
    params: dict
        parameters to be passed to all functions

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
        
    if "Zeta" in params:
        Zeta = params["Zeta"]
    else:
        Zeta = 1.
    
    # loading photoionization rates
    #
    Pc6 = 0.

    # UV background (values around 1e-13 s^-1)
        
    Plw = nc.Gamma_LW_UVB
    Ph1 = nc.Gamma_H_UVB
    Pg1 = nc.Gamma_He_UVB         
        
    # ADDING THE GALACTIC FLUX CONTRIBUTION 
        
    if f_esc != 0.0:
        Plw += nc.Gamma_LW_1000 * (1000*nc.pc/r)**2  * SFR_pure*f_esc
        Ph1 += nc.Gamma_H_1000 * (1000*nc.pc/r)**2 * SFR_pure*f_esc
        Pg1 += nc.Gamma_He_1000 * (1000*nc.pc/r)**2 * SFR_pure*f_esc

    return gc.frtgetcf_cool(T, n, Zeta, Plw, Ph1, Pg1, Pc6) - gc.frtgetcf_heat(T, n, Zeta, Plw, Ph1, Pg1, Pc6)



# SYSTEM OF EQUATIONS

def diff_system(r, y, params):
    """
    differential system for the Euler equations
    
    Parameters
    ==========
    r: float
        temperature
    y: array
        array of the variables (v, n, T)
    params: dict
        parameters to be passed to all functions

    Returns
    =======
    output from the differential system: array

    """    
    
    # params definition
     
    if "v_c" in params:
        v_c_pure = params["v_c"]
    else: 
        raise ValueError("No v_c given")
    
    # defining the equations 
    
    assert len(y.shape) == 1, 'Dimensionality of the input array is wrong'

    v = y[0] # velocity
    rho = y[1] # density
    T = y[2] # temperature
        
    c_S2 = nc.gamma*nc.knorm * T
    c_T2 = nc.knorm * T
    
    n = rho / (nc.mus * nc.mp) #in cgs unit
    
    q = rho*lamda(T, n, r, params) / (nc.mus * nc.mp)**2
    
    v_c = v_c_pure * 1e5 #cm/s 
    v_e = v_c * np.sqrt(2) #cm/s 

    
    output_a = ((2*v/r)*(c_S2 - v_e**2/4) + (nc.gamma-1)*q) / (v**2-c_S2)
    
    output_b = ((2*rho/r)*(v_e**2/4-v**2) - (nc.gamma-1)*q*rho/v) / (v**2-c_S2)

    output_c = ((nc.gamma-1)*(2*T/r)*(v_e**2/4-v**2) - (nc.gamma-1)*q*(v**2-c_T2)/(nc.knorm*v)) / (v**2-c_S2)        
    
    output = np.asarray([output_a,output_b,output_c])
    
    assert len(output.shape) == 1  , 'Dimensionality of the output array is wrong'
    
    return output



def get_profiles(params, resol=1000):
    """
    computes the profiles for v, n, T as a function of r
    
    Parameters
    ==========
    params: dict
        parameters to be passed to all functions
    resol: int, optional
        number of r-steps
    
    Returns
    =======
    profiles: sol_profiles class element

    """    
    
    # params definition
    
    if "beta" in params:
        beta = params["beta"]
    else: 
        raise ValueError("No beta given")
    
    if "SFR" in params:
        SFR_pure = params["SFR"]
    else: 
        raise ValueError("No SFR given")
        
    if "v_c" in params:
        v_c_pure = params["v_c"]
    else: 
        raise ValueError("No v_c given")
    
    if "f_esc" in params:
        f_esc = params["f_esc"]
    else:
        f_esc = 0.
        
    if "Zeta" in params:
        Zeta = params["Zeta"]
    else:
        Zeta = 1.
    
    if "alfa" in params:
        alfa = params["alfa"]
    else:
        alfa = 1.
        
    if "R_in" in params:
        R_in_pure = params["R_in"]
    else:
        alfa = 0.3
    
    # getting the BC
    
    SFR = SFR_pure/nc.year #1/s
    
    E_SN = 1e51*(SFR)/100 #erg (energy from supernovae-- energy of a single SN \
                        #per 100 M_sun of Star Formation)

    M_dot = beta*SFR*nc.ms  #mass from SN
    
    E_dot = alfa*E_SN #erg/s
  
    #M0 = 1. 
    v0 = np.sqrt(E_dot/M_dot)/np.sqrt(2)  #m/s
    #c0 = v0/M0

    R = R_in_pure*1000*nc.pc #cm

    rho0 = 0.1125395*np.sqrt(M_dot**3/E_dot) / R**2 #g/cm^3
    P0 = 0.0337618*np.sqrt(M_dot*E_dot) / R**2  #g/cm^2    
    T0 = P0/(rho0*nc.knorm)  #K
    
    y0 = np.asarray([v0,rho0,T0]) 
    
    # integrating the equations
    
    r_bound = (R, 100*R)
    
    r_eval = np.linspace(r_bound[0],r_bound[1],resol)

    sol = si.solve_ivp(diff_system, r_bound, y0, t_eval=r_eval, args=(params,))
    
    if sol.success == False:
        print('Error in integration procedure')
    elif sol.success == True:
        print('Integration completed successfully')
    
    r = sol.t #cm
    v = sol.y[0] #cm/s
    rho = sol.y[1]
    n = rho / (nc.mus*nc.mp) #cm^-3
    T = sol.y[2] #K
    
    profiles = []
    
    profiles.append(v)
    profiles.append(n)
    profiles.append(T)
    
    return sol_profiles(radius=r, variables=profiles, params=params)
    


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
    
    h = np.linspace(min(profiles.r),min(max(profiles.r),10.*1000*nc.pc), h_resol)
   
    sigma_CII = np.zeros_like(h)
    
    sigma_CII[:] = 2* np.trapz(epsilon[profiles.r>h[:]] * profiles.r[profiles.r>h[:]] / np.sqrt((profiles.r[profiles.r>h[:]])**2 - h[:]**2),\
                 profiles.r[profiles.r>h[:]])        

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



def get_intensity_convolved(h, intensity_raw, params, params_obs):
    """
    computes the ionization state for hydrogen (x_HI) and Carbon (x_CII)
    
    Parameters
    ==========
    r: array
        
    intensity_raw: array
    
    params_obs: dict
        parameters from observations
    
    Returns
    =======
    intensity_convolved: lum_profile class element

    """    
    
    
    # params definition
    
    if "r_beam" in params_obs:
        if "r_beam" == None:
            raise ValueError("No beam given")
        else:
            r_beam = params_obs["r_beam"]
    else: 
        raise ValueError("No beam given")


    if "beam" in params_obs:
        if "beam" == None:
            raise ValueError("No beam given")
        else:
            beam = params_obs["beam"]
    else: 
        raise ValueError("No beam given")

    
    # creates the 2d profiles
    
    x, y, profile_2d = twod_making(intensity_raw, h)
    x, y, beam_2d = twod_making(beam, r_beam)

    #makes the 2dfft
    
    f_beam = np.fft.fft2(beam_2d)
    f_imag = np.fft.fft2(profile_2d)
                        
    #convolution
    f_cimag = f_beam * f_imag
    cimage = np.real(np.fft.ifftshift(np.fft.ifft2(f_cimag)))
    cprofile_2d = cimage[:, cimage.shape[1]//2]
    
    intensity_convolved = np.interp(h, x, cprofile_2d)
    
    #normalizes the convolved intensity
    
    norm_intensity = np.trapz(2*np.pi * h * intensity_raw, h) / np.trapz(2*np.pi * h * intensity_convolved, h)
    
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
        
    
    #    for el_h in h:
#        
#        integral = np.trapz(epsilon[r>el_h] * r[r>el_h] / np.sqrt((r[r>el_h])**2 - el_h**2), r[r>el_h])        
#        sigma_CII.append(2.*integral)
#                
#    sigma_CII = np.asarray(sigma_CII)
    
    
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


        