# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 15:16:16 2021

@author: anna
"""

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

import matplotlib.pyplot as plt

from sol_modules import get_profiles

from post_sol_modules import get_ionization_states, get_surface_density, get_intensity_raw, \
    get_intensity_convolved, get_chi2

from load_data import observational_data_fuji

import plot_config as pltc

import time

post_profiles = True  # switcher for the steps after the profiles integration (i.e. ionization states, sigma, emission)

# Creating a dictionary for the parameters

# N.B. The integration time is bigger for higher betas, higher SFR, and higher v_c;
#      The solver also slows down significantly for non-zero values of the escape fractions f_esc_ion and f_esc_FUV,
#      and for a NFW DM model;


betas = np.asarray([1.0,2.0,3.0])

params = dict([("DM_model", "NFW"),
               ("beta", 1.0),
               ("SFR", 1000.),
               ("f_esc_ion", 0.0),
               ("f_esc_FUV", 0.0),
               ("v_c", 250.),
               ("redshift", 5.0),
               ("Zeta", 1.0),
               ("alfa", 1.0),
               ("R_in", 0.3)])

fig_sol, axs_sol = pltc.plot_configurator(plot_type="sol")
fig_ion, axs_ion = pltc.plot_configurator(plot_type="ion")

fig_sol.suptitle("v_c = {:.1f} km/s, f_escFUV = 0.0, f_escION = 0.0; NFW model" \
                 .format(params["v_c"], params["SFR"]))
fig_ion.suptitle("v_c = {:.1f} km/s, f_escFUV = 0.0, f_escION = 0.0; NFW model" \
                 .format(params["v_c"], params["SFR"]))

fig_int_raw, ax_int_raw = pltc.plot_configurator(plot_type="int", xlim=10)
fig_int_conv, ax_int_conv = pltc.plot_configurator(plot_type="int", xlim=10)

ax_int_raw.set_ylim((1e-3, 1e2))
ax_int_conv.set_ylim((1e-3, 1e2))

fig_int_raw.suptitle("v_c = {:.1f} km/s, f_escFUV = 0.0, f_escION = 0.0; NFW model" \
                     .format(params["v_c"], params["SFR"]))
fig_int_conv.suptitle("v_c = {:.1f} km/s, f_escFUV = 0.0, f_escION = 0.0; NFW model" \
                      .format(params["v_c"], params["SFR"]))

time_start = time.perf_counter()

print("##################################")
print("Run with the following parameters:")
print("##################################")
print(params)
print("##################################")

# getting the profiles using the integrator in sol_modules.py (this is the step that needs to be optimized)

show_profile = True
alphas = ["-", "--"]

data = observational_data_fuji

for SFR, alpha_counter in zip([500,1000], range(len(alphas))):
    params.update(SFR=SFR)

    for beta_el, beta_counter in zip(betas, range(len(betas))):
        params.update(beta=beta_el)

        profiles = get_profiles(params, resol=1000, print_time=True, integrator="RK45")
        ionization_state = get_ionization_states(profiles, params)
        sigma_CII = get_surface_density(profiles, ionization_state, params, rmax=30, h_resol=500, add_CMB_suppression=True)
        intensity_raw = get_intensity_raw(sigma_CII, params, data.params_obs)
        intensity_conv = get_intensity_convolved(intensity_raw, params, data.params_obs, data,
                                                 add_central_contribution=False)
        if beta_counter == 0 and SFR == 500:
            int_500, = ax_int_conv.plot(intensity_conv.h/(1000*nc.pc),intensity_conv.var, color="gray".format(beta_counter),
                                linestyle=alphas[alpha_counter])

        elif beta_counter == 0 and SFR == 1000:
            int_1000, = ax_int_conv.plot(intensity_conv.h/(1000*nc.pc),intensity_conv.var, color="gray".format(beta_counter),
                                linestyle=alphas[alpha_counter])

        if SFR == 500:
            profiles.plot(ax=axs_sol, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter), linestyle=alphas[alpha_counter])
            ionization_state.plot(ax=axs_ion, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter), linestyle=alphas[alpha_counter])

            intensity_raw.plot(ax=ax_int_raw, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter), linestyle=alphas[alpha_counter])
            intensity_conv.plot(ax=ax_int_conv, label=r"$\beta$={:.1f}".format(beta_el), color="C{}".format(beta_counter), linestyle=alphas[alpha_counter])

        elif SFR == 1000:
            profiles.plot(ax=axs_sol, color="C{}".format(beta_counter), linestyle=alphas[alpha_counter])
            ionization_state.plot(ax=axs_ion, color="C{}".format(beta_counter), linestyle=alphas[alpha_counter])

            intensity_raw.plot(ax=ax_int_raw, color="C{}".format(beta_counter), linestyle=alphas[alpha_counter])
            intensity_conv.plot(ax=ax_int_conv, color="C{}".format(beta_counter), linestyle=alphas[alpha_counter])

fuji = ax_int_conv.errorbar(data.x / (1000 * nc.pc), data.data, yerr=data.err, \
                        markerfacecolor='maroon', markeredgecolor='maroon', marker='o', \
                        linestyle='', ecolor='maroon')

ax_int_conv.legend([fuji, int_500, int_1000], ["Fujimoto+19", "SFR = 500", "SFR = 1000"])  ##

ncol = 5

fig_sol.legend(loc="lower center", fontsize="small", ncol=ncol)
fig_ion.legend(loc="lower center", fontsize="small", ncol=ncol)

fig_int_raw.legend(loc="lower center", fontsize="small", ncol=ncol)
fig_int_conv.legend(loc="lower center", fontsize="small", ncol=ncol)

plt.show()



