# -*- coding: utf-8 -*-
"""
Created on Tue May 11 16:49:20 2021

@author: anna
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy.integrate as si

#from numba import jit

import time

import emcee

import natconst as nc

import mydir
from load_data import obs_data_list, names, names_CII_halo, names_wo_CII_halo,\
                      names_other,  observational_data_fuji

from fast_solver import diff_system_fast, stopping_condition

from radiation_fields import UVB_rates
from my_utils import twod_making

import gnedincooling as gc

gc.frtinitcf(0, os.path.join(mydir.script_dir, "input_data", "cf_table.I2.dat"))        

# preliminar parameters



data = obs_data_list[1]


redshift = data.params_obs["redshift"]

Plw = UVB_rates(redshift, quantity="LW rate")
Ph1 = UVB_rates(redshift, quantity="H rate")
Pg1 = UVB_rates(redshift, quantity="He rate")  
gamma_H = Ph1
gamma_CI = UVB_rates(redshift, quantity="CI rate")
gamma_CII = UVB_rates(redshift, quantity="CII rate")
UV_intensity = UVB_rates(redshift, quantity="UV intensity")

I_CMB = (2*nc.hh*nc.line_CII_rest_frame**3) /  \
        (nc.cc**2 * (np.exp((nc.hh*nc.line_CII_rest_frame) / (nc.kk*(1.+redshift)*nc.CMB_temperature)) - 1.))
        
B_coeff = nc.A_coeff_ul * nc.cc**2 / (2*nc.hh*nc.line_CII_rest_frame**3)

A_tilde = nc.A_coeff_kl*nc.A_coeff_ku / (nc.A_coeff_kl + nc.A_coeff_ku)
P_UV_ul = A_tilde * nc.cc**2 / (2*nc.hh*nc.line_UV_mixing_rest_frame**3) 
    
T_UV = 3.61e4
    
rmax = 30
h_resol = 500
r_resol = 500

cut = 45.

integrator = "RK45"
    

other_params = dict([("integrator", integrator),
                   ("cut", cut),
                   ("rmax", rmax),
                   ("h_resol", h_resol),
                   ("r_resol", r_resol),
                   ("redshift", redshift),
                   ("Plw", Plw),
                   ("Ph1", Ph1),
                   ("Pg1", Pg1),
                   ("gamma_H", gamma_H),
                   ("gamma_CI", gamma_CI),
                   ("gamma_CII", gamma_CII),
                   ("UV_intensity", UV_intensity),
                   ("FWHM_vel", data.params_obs["line_FWHM"]),
                   ("I_CMB", I_CMB),
                   ("B_coeff", B_coeff),
                   ("P_UV_ul", P_UV_ul),            
                   ("T_UV", T_UV)])
                
    


# MCMC definitions

"""
theta : beta, SFR, v_c
"""

def get_emission_fast(theta, data, other_params, print_time = True):
    
    # setting the parameters
    
    beta, SFR_pure, v_c_pure = theta
    
    f_esc_ion = 0.
    f_esc_FUV = 0.
    Zeta = 1. 
    R_in_pure = 0.3
    alfa = 1.
    
    Plw = other_params["Plw"]
    Ph1 = other_params["Ph1"]
    Pg1 = other_params["Pg1"]
    
    
    redshift = other_params["redshift"]
    
    
    overdensity = 200.
    hubble2 =  2.1962761244736533e-18**2 * (0.30712*(1+redshift)**3 + 5.384308416949404e-05*(1+redshift)**4 +  0.6913912010962934)
    critical_density = 3*hubble2 / (8*np.pi*nc.gg) # in g/cm^3
    M_vir_pure = (v_c_pure*1e5)**3 / nc.gg**1.5 * (4*np.pi*overdensity*critical_density/3)**(-0.5) / nc.ms

    # gravity part
        
    cosmo_h = 0.6774
    logc =  0.537 + (1.025 - 0.537) * np.exp(-0.718*redshift**1.08) + (-0.097 + 0.024 * redshift) * np.log10(M_vir_pure*cosmo_h/1e12)
    c = 10**logc
    A_NFW = np.log(1+c) - c/(1.+c)
    r_s = np.cbrt(3*M_vir_pure*nc.ms/(critical_density * 4*np.pi*overdensity)) / c / 1e3 / nc.pc # in kpc

    
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
    
    # changing the dimensions for the integration part
    
    R_kpc = R_in_pure # in kpc
    
    v0_kms = v0 / 1e5 # in km/s
    n0_cm3 = rho0 / (nc.mp*nc.mus) # in cm-3
    T0_K = T0 # in K
    
    y0 = np.asarray([v0_kms,n0_cm3,T0_K]) 
    
    # integrating the equations
    
    r_bound = (R_kpc, 100*R_kpc)
    
    r_eval = np.linspace(r_bound[0],r_bound[1],other_params["r_resol"])
    
    if print_time:
          
        t_ivp = time.perf_counter()

    cut = other_params["cut"]
    sol = si.solve_ivp(diff_system_fast, r_bound, y0, t_eval=r_eval,\
                       args=( SFR_pure, redshift, M_vir_pure, f_esc_ion, f_esc_FUV, Plw, Ph1, Pg1, Zeta, A_NFW, r_s, cut),\
                       method = other_params["integrator"], events=stopping_condition) #,rtol=1.0e-3
        
    if print_time:
        time_ivp = (time.perf_counter() - t_ivp)
        print("total time ivp (s)=", time_ivp)

    r_kpc = sol.t # in kpc
    v_kms = sol.y[0] # in km/s
    n_cm3 = sol.y[1] # in cm-3
    T_K = abs(sol.y[2]) # in K
        
    # getting back to the old dimensions
    
    v = v_kms * 1e5    # in cm/s
    n = n_cm3 # in cm-3
    T = T_K # in K
    
    
    # ionization part
            
    
    gamma_CI = other_params["gamma_CI"] + nc.Gamma_CI_FUV_1000 * (1./r_kpc)**2 * SFR_pure*f_esc_FUV
    gamma_CI += nc.Gamma_CI_EUV_1000 * (1./r_kpc)**2 * SFR_pure*f_esc_ion
    gamma_H = other_params["gamma_H"] + nc.Gamma_H_1000 * (1./r_kpc)**2 * SFR_pure*f_esc_ion        
    gamma_CII = other_params["gamma_CII"] + nc.Gamma_CII_1000 * (1./r_kpc)**2 * SFR_pure*f_esc_ion
    intensity_tot = other_params["UV_intensity"] + nc.intensity_UV_1000 * (1./r_kpc)**2 * SFR_pure*f_esc_FUV

        
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
    n_H = (1. - x_e) * nc.A_H * n

    x_CII = 1. / (1. + (gamma_CII)/(beta_CIII * n_e) + (beta_CII * n_e)/(gamma_CI + kappa_CI * n_e) + kappa_CII/beta_CIII)

    # emission part
    
    epsilon = 7.9e-20 *  n**2 * (nc.A_C * Zeta) * nc.A_H * x_e * x_CII * np.exp(-92./T) /  92.**0.5        

    C_ul_e = 8.63e-6 / 2. / np.sqrt(T) * 1.60
    C_ul_H = 20. * 1e-10 / 1.3

    B_coeff = other_params["B_coeff"]
    P_UV_ul = other_params["P_UV_ul"]
    T_UV = other_params["T_UV"]
    
    num = B_coeff * I_CMB + nc.A_coeff_ul + P_UV_ul * intensity_tot + n_e * C_ul_e + n_H * C_ul_H 
    den =  B_coeff * I_CMB + P_UV_ul * intensity_tot * np.exp(-nc.T_star/T_UV) +\
            n_e * C_ul_e * np.exp(-nc.T_star/T) + n_H * C_ul_H * np.exp(-nc.T_star/T)

    T_spin = nc.T_star / np.log(num/den)

    eta = 1. - (np.exp(nc.hh*nc.line_CII_rest_frame / (nc.kk*(1.+redshift)*T_spin)) - 1.) /\
           (np.exp((nc.hh*nc.line_CII_rest_frame) / (nc.kk*(1.+redshift)**2*nc.CMB_temperature)) - 1.)
        
    epsilon *= eta

    # intergrating along the line of sight
    
    h = np.linspace(min(r_kpc),max(max(r_kpc),rmax), h_resol)
   
    sigma_CII = np.zeros_like(h)
    
    for el_h, i_h in zip(h, range(len(h))):
        
        integral = np.trapz(epsilon[r_kpc>el_h] * nc.pc * 1e3 * r_kpc[r_kpc>el_h] / np.sqrt((r_kpc[r_kpc>el_h])**2 - el_h**2), r_kpc[r_kpc>el_h])                
        sigma_CII[i_h] = 2.*integral
                
    
    # transforming sigma to the intensity
    
    FWHM_vel = other_params["FWHM_vel"]
    intensity_raw = sigma_CII / (nc.line_CII_rest_frame*4*np.pi*(1.+redshift)**3*FWHM_vel/nc.cc)
        
    # changing the units
            
    intensity_raw *= 1e26 #transformation to mJy 
    intensity_raw /= 4.2e10 #transforming sr to arcsec^2
    
        
    x, y, profile_2d = twod_making(intensity_raw, h, nimage=1000)
    
    beam_interp = np.interp(h, data.x_beam/1e3/nc.pc, data.beam, right=0.)
    
    beam_interp[beam_interp<0.] = 0.

    x_beam, y_beam, beam_2d = twod_making(beam_interp, h, nimage=1000)
    
    #makes the 2dfft
    
    f_beam = np.fft.fft2(beam_2d)
    f_imag = np.fft.fft2(profile_2d)
                        
    #convolution
    f_cimag = f_beam * f_imag
    cimage = np.real(np.fft.ifftshift(np.fft.ifft2(f_cimag)))
    cprofile_2d = cimage[:, cimage.shape[1]//2]
    
    intensity_convolved = np.interp(h, x, cprofile_2d, right=0.)
    
    #normalizes the convolved intensity
    
    norm_intensity = np.trapz(2*np.pi * h * intensity_raw, h)\
                                / np.trapz(2*np.pi * h * intensity_convolved, h)
    
    intensity_convolved *= norm_intensity

    return h, intensity_convolved
    
def log_likelihood(theta, data, other_params):
    
    
    h, intensity_convolved = get_emission_fast(theta, data, other_params)
    
    emission_profile = interp1d(h, intensity_convolved)

    res = emission_profile(data.x/1e3/nc.pc) - data.data
    
    residuals = 2*res / (data.err_down + data.err_up) 
    
    chi2 = np.sum(residuals**2)
    
    log_likelihood.counter += 1

    print("iteration {}: chi2 = {:.1f} fot beta = {:.1f}, SFR = {:.1f}, v_c = {:.1f}".format(\
          log_likelihood.counter, chi2, theta[0], theta[1], theta[2]))
    
    
    return -0.5 * chi2

def log_prior_uniform(theta, data):
    """
    defines the priors for the set of parameters theta
    
    Parameters
    ==========
    theta: array
            
    Returns
    =======
    priors value: float

    """    
    
    beta, SFR, v_c = theta
    
    
    if 1.0 < beta < 6.0 and 10. < SFR < 100.0 and 150. < v_c < 250.:
        return 0.0
    
    return -np.inf


def log_prior_gaussian(theta, data):
    """
    defines the priors for the set of parameters theta
    
    Parameters
    ==========
    theta: array
            
    Returns
    =======
    priors value: float

    """    
    
    beta, SFR, v_c = theta
    
    if 1.0 < beta < 8.0:
        prior =  0.0
    else:
        prior = -np.inf
        
    prior += - 2*(SFR-data.params_obs["SFR"])**2/(data.params_obs["SFR_err_up"]+data.params_obs["SFR_err_down"])**2
    prior += - 2*(SFR-data.params_obs["v_c"])**2/(data.params_obs["v_c_err_up"]+data.params_obs["v_c_err_down"])**2


def log_probability(theta, data, other_params):
    """
    defines the probability as a combination of priors and likelihood
    
    Parameters
    ==========
    theta: array
    
    obs_data: array
            
    Returns
    =======
    priors value: float

    """    
    lp = log_prior_gaussian(theta, data)
    
    if not np.isfinite(lp):
        return -np.inf
    
    return lp + log_likelihood(theta, data, other_params)


if __name__ == "__main__":
    
    theta_true = [4.0, 50., 200.]
    
    ndim = len(theta_true)
    nwalkers= 32
    
    pos = theta_true + np.asarray([1.0, 50., 50.]) * np.random.randn(nwalkers, ndim)
    
    pos[pos < 0] = 5.


    folder = "data_emcee"
    
    if not os.path.exists(os.path.join(mydir.data_dir, folder)):
        os.mkdir(os.path.join(mydir.data_dir, folder))

    filename = os.path.join(mydir.data_dir, folder, "trial_run.h5")
    
    
    backend = emcee.backends.HDFBackend(filename)
    backend.reset(nwalkers, ndim)
        
    sampler = emcee.EnsembleSampler(nwalkers=32, ndim=ndim, log_prob_fn=log_probability,\
                                    args=(data,other_params), backend = backend)
    sampler.run_mcmc(pos, 100, progress=True);
    
    
    fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    labels = ["beta", "SFR", "v_c"]
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)
        
    axes[-1].set_xlabel("step number")
    
    folder = "plot_emcee"
    
    if not os.path.exists(os.path.join(mydir.plot_dir, folder)):
        os.mkdir(os.path.join(mydir.plot_dir, folder))

    plt.savefig(os.path.join(mydir.plot_dir, folder, "emission.png"))
    
    print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
    
    print(
    "Mean autocorrelation time: {0:.3f} steps".format(
    np.mean(sampler.get_autocorr_time())
    )
    )
    
    
    
    
    
    
