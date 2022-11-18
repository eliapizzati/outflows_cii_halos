"""
This script contains the following functions:
    mstar_behroozi_from_file: M_star from Behroozi+13, from file; it uses z=5 as default and it can't do any other redshift
    mstar_behroozi: M_star from Behroozi+13 direclty from fit
    mvir_behroozi: inverse of mstar_behroozi, directly from fit
    twod_making: makes a 2D image from a 1D profile; the image will be used to be transformed for the convolution
    to_habing_flux: converts from UV flux to Habing flux
    from_habing_flux: converts from Habing flux to UV flux
    plw_from_habing_flux:  converts from Habing flux to PLW flux
    get_concentration: returns the concentration of a halo
    get_virial_radius: returns the virial radius given the virial mass
    get_mass_profile_NFW: returns the mass profile of a NFW halo
    get_mass_profile_disk: returns the mass profile of a disk
    get_circular_velocity_profile_NFW: returns the circular velocity inside a radius r using a NFW profile
    get_circular_velocity_profile_NFW_and_disk: returns the circular vvelocity inside a radius r using a NFW profile + the central disk
    get_virial_mass_from_vc: gets the virial mass from the parameter v_c (maximun circular velocity)
    get_vc_from_virial_mass: gets the parameter v_c (maximun circular velocity) from the virial mass
    from_xfitted_to_eta: converts from the fitted x parameter to the eta parameter
    from_eta_to_xfitted: converts from the eta parameter to the fitted x parameter

"""


import numpy as np
import os
import scipy.interpolate
import natconst as nc

import my_dir

from scipy import optimize

from astropy.cosmology import Planck15 as cosmo



def mstar_behroozi_from_file(M_vir):
    """
    M_star from Behroozi+13, from file; it uses z=5 as default and it can't do any other redshift
    Parameters
    ----------
    M_vir:
        halo mass in solar masses

    Returns
    -------
    float
        galaxy stellar mass in solar masses
    """
    input_filename_behroozi = os.path.join(my_dir.script_dir, "input_data", "behroozi_z_5.dat")

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
    """
    M_star from Behroozi+13 direclty from fit
    # https://arxiv.org/pdf/1207.6105.pdf

    Parameters
    ----------
    M_vir: float
        halo mass in solar masses
    z: float
        redshift

    Returns
    -------
    float
        galaxy stellar mass in solar masses
    """

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


def mvir_behroozi(M_star, z=5.0):
    """
    inverse of mstar_behroozi, directly from fit
    Parameters
    ----------
    M_star:     float
        galaxy stellar mass in solar masses
    z:        float
        redshift

    Returns
    -------
    float
        halo mass in solar masses
    """
    def zero_me(x, M_star):
        return M_star - mstar_behroozi(x, z=5.)

    mvir = optimize.newton(zero_me, x0=15 * M_star, args=(M_star,), maxiter=500, tol=1.)

    return mvir



def twod_making(profile, x_axis, nimage=1000):
    """
    makes a 2D image from a 1D profile; the image will be used to be transformed for the convolution

    Parameters
    ----------
    profile:    array
        1D profile
    x_axis:   array
    nimage:     int
        number of pixels of the image

    Returns
    -------
    array, array
        extended x_axis, 2D image
    """
    x_ext = np.linspace(-x_axis.max(), x_axis.max(), nimage)
    
    xy, yx = np.meshgrid(x_ext,x_ext)
                        
    function_ext = scipy.interpolate.interp1d(x_axis, profile, \
                                              fill_value = (profile[0], 0.), bounds_error=False)
        
    z = function_ext(np.sqrt(xy**2 + yx**2))
    
    return x_ext, z


def to_habing_flux(I_UV):
    """
    converts from UV flux to Habing flux
    Parameters
    ----------
    I_UV: float
        UV flux in erg/s/cm^2/Hz

    Returns
    -------
    float
        Habing flux in MW value
    """
    return  4*np.pi*I_UV*(13.6-6.)*nc.ev/nc.hh/nc.cc/5.29e-14


def from_habing_flux(I_g0):
    """
    converts from Habing flux to UV flux
    Parameters
    ----------
    I_g0: float
        Habing flux in MW value

    Returns
    -------
    float
        UV flux in erg/s/cm^2/Hz
    """
    return I_g0 * nc.hh * nc.cc * 5.29e-14 / 4 / np.pi / (13.6 - 6.) / nc.ev


def plw_from_habing_flux(I_g0):
    """
    converts from Habing flux to UV flux
    Parameters
    ----------
    I_g0

    Returns
    -------

    """
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
    
    overdensity: float, optional
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
    """
    gets the mass inside a radius r using a disk profile (exponential or uniform)

    Parameters
    ----------
    r: float
        radius of interest in cm
    M_star: float
        stellar mass in solar masses
    model: str, optional
        disk model (default is "exp")
    R_pure: float, optional
        scale radius of the exponential disk in kpc unit (default is 0.3 )

    Returns
    -------
    float
        mass inside r in solar masses
    """
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
    returns the circular velocity inside a radius r using a NFW profile
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
    returns the circular velocity inside a radius r using a NFW profile + disk

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


def from_xfitted_to_eta(x, a_fit, b_fit):
    """
    gives the eta parameter from the fitted x parameter (which can be either SFR or Mstar), using the
     relation between eta and SFR/Mstar in the paper (log eta = a_fit + b_fit * log SFR/1msun/yr or
    log eta = a_fit + b_fit * log Mstar/1e10msun)

    Parameters
    ----------
    x: float
        SFR or Mstar in solar masses per year or solar masses
    a_fit : float
        a_fit parameter from the fit
    b_fit : float
        b_fit parameter from the fit

    Returns
    -------
    float
        eta parameter
    """
    return 10**(b_fit*np.log10(x) + np.log10(a_fit))


def from_eta_to_xfitted(eta, a_fit, b_fit):
    """
    gives the x parameter (which can be either SFR or Mstar) from the eta paramter, using the relation
    between eta and SFR/Mstar in the paper (log eta = a_fit + b_fit * log SFR/1msun/yr or
    log eta = a_fit + b_fit * log Mstar/1e10msun)

    Parameters
    ----------
    eta: float
        eta parameter
    a_fit: float
        a_fit parameter from the fit
    b_fit: float
        b_fit parameter from the fit

    Returns
    -------
    float
        SFR or Mstar in solar masses per year or solar masses
    """
    return 10 ** ((np.log10(eta) - np.log10(a_fit))/ b_fit)

