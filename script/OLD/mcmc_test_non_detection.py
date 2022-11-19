
"""
This script is used to test the comparison between detections and non detections
NOT SURE THIS IS STILL UPDATED AND WORKS
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import itertools
from scipy.interpolate import interp1d
from collections import namedtuple

import plot_config as pltc
import my_dir
import natconst as nc

from load_data import obs_data_list_non_det, obs_data_list, names, names_wo_CII_halo, names_wo_CII_halo_short

data0 = obs_data_list_non_det[0]
data1 = obs_data_list[1]


datas = obs_data_list_non_det# [data0, data1]
colors = ["maroon", "navy", "olive", "magenta", "purple"]
print(datas)

"""
PLOT1
"""


from mcmc_main import get_emission_fast, get_other_params
from mcmc_main import h, grid

fig_int_conv, ax_int_conv = pltc.plot_configurator(plot_type="int", xlim=15)

ax_int_conv.set_ylim((1e-3, 1e2))


fig_sigma, ax_sigma = pltc.plot_configurator(plot_type="sigma", xlim=15)


for data, color in zip(datas, colors):

    other_params = get_other_params(data.params_obs["redshift"], data.params_obs["line_FWHM"])

    beam_interp = np.interp(h, data.x_beam / 1e3 / nc.pc, data.beam, right=0.)

    beam_interp[beam_interp < 0.] = 0.

    beam_func = interp1d(h, beam_interp, \
                         fill_value=(beam_interp[0], 0.), bounds_error=False)

    beam_2d = beam_func(np.sqrt(grid[0] ** 2 + grid[1] ** 2))
    f_beam = np.fft.fft2(beam_2d)


    log_mstar = np.log10(data.params_obs["M_star"])
    log_sfr = data.params_obs["log_SFR"]
    log_vc = data.params_obs["log_v_c"]

    log_beta_fit = -0.36 *(log_mstar - 10) + np.log10(4.9)


    for log_beta in [log_beta_fit]:
        theta = [log_beta, log_sfr, log_vc]

        intensity, sigma = get_emission_fast(theta, data, other_params, h, grid, f_beam, return_quantities="emission")

        ax_int_conv.plot(h, intensity, alpha=1, label="eta = {:.1f}".format(10**log_beta), color=color)
        ax_sigma.plot(h, sigma, alpha=1, label="eta = {:.1f}".format(10**log_beta), color=color)


    alpine = ax_int_conv.errorbar(data.x / (1000 * nc.pc), data.data, yerr=data.err, \
                                  markerfacecolor=color, markeredgecolor=color, marker='o', \
                                  linestyle='', ecolor=color)
    factor_data = data.data[0] / data.beam[0]

    ax_int_conv.plot(data.x_beam / 1e3 / nc.pc, data.beam * factor_data, linestyle="--", color="gray")


fig_int_conv.legend(loc="lower center", ncol=8, fontsize="small")

fig_sigma.legend(loc="lower center", ncol=8, fontsize="small")

plt.show()


