"""
This script generates some profiles (and post processes too for emission etc)
and stores data in files (no plots so that it can be run in machines too)
"""

import os
import numpy as np

import mydir

from sol_modules import get_profiles

from post_sol_modules import get_ionization_states, get_surface_density, get_intensity_raw, \
    get_intensity_convolved, get_chi2

from load_data import obs_data_list, observational_data_fuji

import time

post_profiles = False  # switcher for the steps after the profiles integration (i.e. ionization states, sigma, emission)
plotting = True

# Creating a dictionary for the parameters


params = dict([("DM_model", "NFW+disk"),
               ("beta", 4.4),
               ("SFR", 88.),
               ("f_esc_ion", 0.0),
               ("f_esc_FUV", 0.0),
               ("v_c", 217.),
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

resol = 10000

time_profile = time.perf_counter()
profiles = get_profiles(params, resol=resol, print_time=False, integrator="RK45")
time_profile = (time.perf_counter() - time_profile)
print("total profile time (s)=", time_profile)


if profiles.check_nans() == True:
    string_nans = "Integration error: presence of nans"
else:
    string_nans = "Integration successful"
print(string_nans)

profiles.to_file(attributes_in_name=f"DC881725_{resol}_per_laura")

if post_profiles:
    # getting the ionization states, the sigma CII and the convolved intensity for the CII

    ionization_state = get_ionization_states(profiles, params)

    sigma_CII = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500, add_CMB_suppression=True)

    # here we need to specify some data; we import data from Fujimoto+19

    data = observational_data_fuji

    intensity_raw = get_intensity_raw(sigma_CII, params, data.params_obs)

    intensity_conv = get_intensity_convolved(intensity_raw, params, data.params_obs, data,
                                             add_central_contribution=False)

    chi2 = get_chi2(intensity_conv, data)
    intensity_conv.to_file()
    intensity_raw.to_file()


if plotting:
    profiles.plot(savefig=True)
    if post_profiles:
        ionization_state.plot(savefig=True)
        sigma_CII.plot(savefig=True)
        intensity_conv.plot(savefig=True)
        intensity_raw.plot(savefig=True)

time_elapsed = (time.perf_counter() - time_start)
print("total time elapsed (s)=", time_elapsed)



