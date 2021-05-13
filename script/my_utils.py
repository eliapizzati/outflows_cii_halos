# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 14:50:57 2021

@author: anna
"""


import numpy as np

import scipy.interpolate
import natconst as nc

from astropy.cosmology import Planck13 as cosmo



def sersic(r, r_e, sersic_index, central):
        
    b = 2 * sersic_index - 1/3
    
    I_e = central * np.exp(b*((1/r_e)**(1/sersic_index)-1))
    
    return I_e*np.exp(-b*((r/r_e)**(1/sersic_index)-1.))

def halo(r, r_n, central):
        
    C = central * np.exp(1/r_n)
    
    return C*np.exp(-r/r_n)


def twod_making(profile, x_axis, nimage=1000):

    x_ext = np.linspace(-x_axis.max(), x_axis.max(), nimage)
    
    xy, yx = np.meshgrid(x_ext,x_ext)
                        
    function_ext = scipy.interpolate.interp1d(x_axis, profile, \
                                              fill_value = (profile[0], 0.), bounds_error=False)
        
    z = function_ext(np.sqrt(xy**2 + yx**2))
    
    return x_ext, z


def to_habing_flux(I_UV):
    
    return  4*np.pi*I_UV*(13.6-6.)*nc.ev/nc.hh/nc.cc/5.29e-14



def get_concentration(M_vir, z):
    """
    concentration parameter from Dutton&Macciò2014
    
    Parameters
    ==========    
    M_vir: float
        halo mass in solar masses
    
    z: float
        redshift
    
    Returns
    =======
    float
        concentration parameter

    """    

    b = -0.097 + 0.024 * z
    a = 0.537 + (1.025 - 0.537) * np.exp(-0.718*z**1.08)
    
    logc = a + b * np.log10(M_vir*cosmo.h/1e12)
    
    return 10**logc

def get_virial_radius(M_vir, z, overdensity=200.):
    """
    gets the virial radius given the virial mass 
    
    Parameters
    ==========    
    M_vir: float
        halo mass in solar masses
    
    z: float
        redshift
    
    overdensity: float
        overdensity used for the halo definition (default definition is 200)

    
    Returns
    =======
    float
        virial radius in cm

    """    
    return np.cbrt(3*M_vir*nc.ms/(cosmo.critical_density(z).value * 4*np.pi*overdensity))


def get_mass_profile_NFW(r,M_vir,z):
    """
    gets the mass inside a radius r using a NFW profile
    
    Parameters
    ==========    
    r: float
        radius of interest in cm
        
    M_vir: float
        halo mass in solar masses
    
    z: float
        redshift
    
    Returns
    =======
    float
        mass inside r in solar masses

    """
    c = get_concentration(M_vir, z)
    
    A_NFW = np.log(1+c) - c/(1.+c)
    
    r_s = get_virial_radius(M_vir,z) / c 
    
    M = M_vir/A_NFW * (np.log(1.+r/r_s)+r_s/(r_s+r) - 1)
    
    return M
    

    

def get_circular_velocity_profile_NFW(r,M_vir,z):
    """
    gets the mass inside a radius r using a NFW profile
    
    Parameters
    ==========    
    r: float
        radius of interest in cm
        
    M_vir: float
        halo mass in solar masses
    
    z: float
        redshift
    
    Returns
    =======
    float
        circular velocity at r in cm/s

    """
    M_r = get_mass_profile_NFW(r,M_vir,z)
    
    return np.sqrt(nc.gg*M_r*nc.ms/r)
    

def get_virial_mass_from_vc(v_c, z, overdensity=200.):
    """
    gets the virial mass from the parameter v_c (maximun circular velocity)
    
    Parameters
    ==========    
    v_c: float
        maximum circular velocity in cm/s
    
    z: float
        redshift
    
    overdensity: float
        overdensity used for the halo definition (default definition is 200)
    
    Returns
    =======
    float
        virial mass in solar masses

    """
    return v_c**3 / nc.gg**1.5 * (4*np.pi*overdensity*cosmo.critical_density(z).value/3)**(-0.5) / nc.ms


def get_vc_from_virial_mass(M_vir, z):
    """
    gets the parameter v_c (maximun circular velocity) from the virial mass
    
    Parameters
    ==========    
    M_vir: float
        halo mass in solar masses
    
    z: float
        redshift
        
    Returns
    =======
    float
        v_c in cm/s

    """
    R_vir = get_virial_radius(M_vir, z)
    
    return np.sqrt(nc.gg*M_vir*nc.ms/R_vir)