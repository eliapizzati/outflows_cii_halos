"""
PARAMETER DESCRIPTION

 - PARAMS

     - DM_model: dark matter profile (choose between "NFW" and "iso_shpere"; if None, "iso_sphere" is considered as default)

     - beta: beta parameter

     - SFR: star formation rate in Msun/yr

     - f_esc_ion: escape fraction of ionizing photons

     - f_esc_FUV: escape fraction of non-ionizing photons

     - v_c: circular velocity in km/s (if the DM profile is NFW then this represents the maximum velocity,
                                       since the circular velocity profile becomes radius-dependent)

     - M_vir: virial mass in solar masses (if the DM profile is "iso_sphere" this is not used; if it is "NFW", then at least
                                           one parameter between v_c and M_vir must be given)

     - redshift: redshift

     - Zeta: metallicity (solar units)

     - alfa: alpha parameter

     - R_in: inner radius in kpc



"""

import numpy as np

import natconst as nc

from sol_modules import get_profiles

from post_sol_modules import get_ionization_states, get_surface_density, get_intensity_raw, \
    get_intensity_convolved, get_chi2

from load_data import observational_data_fuji

import time

post_profiles = True  # switcher for the steps after the profiles integration (i.e. ionization states, sigma, emission)

# Creating a dictionary for the parameters


params = dict([("DM_model", "NFW"),
               ("beta", 3.0),
               ("SFR", 1000.),
               ("f_esc_ion", 0.0),
               ("f_esc_FUV", 0.0),
               ("v_c", 250.),
               ("redshift", 5.0),
               ("Zeta", 1.0),
               ("alfa", 1.0),
               ("R_in", 0.3)])

time_start = time.perf_counter()

print("##################################")
print("Run with the following parameters:")
print("##################################")
print(params)
print("##################################")

# getting the profiles using the integrator in sol_modules.py (this is the step that needs to be optimized)

import matplotlib.pyplot as plt

plotting = False

time_profile = time.perf_counter()
profiles = get_profiles(params, resol=1000, print_time=False, integrator="RK45")
time_profile = (time.perf_counter() - time_profile)


if profiles.check_nans() == True:
    string_nans = "Integration error: presence of nans"
else:
    string_nans = "Integration successful"

print(string_nans)
print("total profile time (s)=", time_profile)

# getting the ionization states, the sigma CII and the convolved intensity for the CII

if post_profiles:
    ionization_state = get_ionization_states(profiles, params)

    sigma_CII = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500, add_CMB_suppression=True)

    # here we need to specify some data; we import data from Fujimoto+19

    data = observational_data_fuji

    intensity_raw = get_intensity_raw(sigma_CII, params, data.params_obs)

    intensity_conv = get_intensity_convolved(intensity_raw, params, data.params_obs, data,
                                             add_central_contribution=False)

    chi2 = get_chi2(intensity_conv, data)


if plotting:
    profiles.plot()
    ionization_state.plot()
    intensity_conv.plot()

time_elapsed = (time.perf_counter() - time_start)

print("total time elapsed (s)=", time_elapsed)

intensity_conv.to_file()


