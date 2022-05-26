# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 14:50:57 2021

@author: anna
"""


import numpy as np
import os
import scipy.interpolate
import natconst as nc

import mydir

from astropy.cosmology import Planck15 as cosmo

input_filename_behroozi = os.path.join(mydir.script_dir, "input_data", "behroozi_z_5.dat")


def mstar_behroozi_from_file(M_vir):

    # First, you can open the file and check the information in it

    halo_masses, halo_mass_ratio, halo_mass_ratio_err_up, halo_mass_ratio_err_down = np.loadtxt(input_filename_behroozi, \
                                                                                                unpack=True)

    stellar_masses = halo_masses + halo_mass_ratio

    stellar_masses_fine = np.linspace(9.,12.,1000)

    halo_masses_fine = np.interp(stellar_masses_fine, stellar_masses, halo_masses)

    index_masses = np.searchsorted(halo_masses_fine, np.log10(M_vir))

    log_M_halo = halo_masses_fine[index_masses]
    log_M_star = stellar_masses_fine[index_masses]

    return 10**log_M_star


def mstar_behroozi(M_vir, z=5.0):
    # Behroozi+13
    # https://arxiv.org/pdf/1207.6105.pdf

    aexp = 1.0 / (1.0 + z)

    nu = np.exp(-4.0 * aexp ** 2)
    log_epsilon = -1.777 + (-0.006 * (aexp - 1) + 0.0 * z) * nu - 0.119 * (aexp - 1)
    log_M1 = 11.514 + (-1.793 * (aexp - 1) - 0.251 * z) * nu
    alpha = -1.412 + (0.731 * (aexp - 1)) * nu
    delta = 3.508 + (2.608 * (aexp - 1) - 0.043 * z) * nu
    gamma = 0.316 + (1.319 * (aexp - 1) + 0.279 * z) * nu

    M1 = 10 ** log_M1
    epsilon = 10 ** log_epsilon

    # eq. 3
    def f_be(x, alpha=0.0, delta=0.0, gamma=0.0):
        return - np.log10(10.0 ** (alpha * x) + 1.0) + delta * ((np.log10(1.0 + np.exp(x))) ** gamma) / (
                    1.0 + np.exp(10 ** (-x)))

    #
    out = np.log10(epsilon * M1) + f_be(x=np.log10(M_vir / M1), alpha=alpha, gamma=gamma, delta=delta) - f_be(x=0.0,
                                                                                                           alpha=alpha,
                                                                                                           gamma=gamma,
                                                                                                           delta=delta)

    return 10**out


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


def from_habing_flux(I_g0):
    return I_g0 * nc.hh * nc.cc * 5.29e-14 / 4 / np.pi / (13.6 - 6.) / nc.ev

def plw_from_habing_flux(I_g0):
    I_UV = I_g0 * nc.hh * nc.cc * 5.29e-14 / 4 / np.pi / (13.6 - 6.) / nc.ev #erg/s/Hz/sr/cm2
    plw = 1.38e9 * I_UV #s^-1
    return plw

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

    r_vir = get_virial_radius(M_vir,z)
    r_s = r_vir / c

    return M_vir/A_NFW * (np.log(1.+r/r_s)+r_s/(r_s+r) - 1)


    
def get_mass_profile_disk(r, M_star, model="exp", R_pure = 0.3):

    if model == "exp":

        R = R_pure * 1e3 * nc.pc

        x = r/R

        return (M_star/2.) * (2. - (x*(x+2)+2)*np.exp(-x))

    elif model == "uniform":
        return M_star

    else:
        raise ValueError("model not supported")



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


def get_circular_velocity_profile_NFW_and_disk(r, M_vir, z):
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
    M_r = get_mass_profile_NFW(r, M_vir, z)

    M_star = mstar_behroozi(M_vir)

    r_vir = get_virial_radius(M_vir, z)

    M_disk_r = get_mass_profile_disk(r, M_star, model = "exp", R_pure=r_vir/nc.pc/1e3/1e2)

    return np.sqrt(nc.gg * (M_r+M_disk_r) * nc.ms / r)


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


#print(mstar_behroozi(6.8e11, z=4.55)/1e10)
#print(get_vc_from_virial_mass(6.8e11, z=4.55)/1e5)