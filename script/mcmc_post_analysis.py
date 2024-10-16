# -*- coding: utf-8 -*-
"""
Created on Thu May 13 14:21:16 2021

@author: anna
"""

import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import corner
import emcee
import itertools
from collections import namedtuple


import plot_config as pltc
import my_dir
import natconst as nc

from load_data import obs_data_list, names, names_CII_halo, observational_data_fuji


stored_data_loc = "mac" # can be either quasar (for the machine is quasar in leiden) or mac (for mac laptop) or github (for github)
                        # or linux (for linux laptop)plot_logprob = False

model = "old"

print_autocorr_time = False
plot_logprob = False
plot_saving = True

# which plots (analysis) to perform
plot1 = True # chains of the walkers
plot2 = True # corner plot of the parameters
plot3 = False # emission from the sampled posterior (requires gnedincooling) + luminosity of the CII line

thin = 1
discard = 0

nwalkers = 48
nsteps = 100000#int(input("number of steps?"))


int_data = int(input("data number?"))
data = obs_data_list[int_data]

filename = "{}_{:.0f}_updated_{}".format(data.params_obs["name_short"], nsteps, model)

print("###################################################################")

print("postprocessing an MCMC with the following params:")
print("n steps = {}".format(nsteps))
print("n walkers = {}".format(nwalkers))
print("data object = {}".format(data.params_obs["name_short"]))
print("filename = {}".format(filename))

print("###################################################################")

if stored_data_loc == "quasar":
    path = os.path.join("/data2/pizzati/projects/outflows/data_mcmc", "{}.h5".format(filename))
elif stored_data_loc == "mac" or stored_data_loc == "local": #default local is the mac laptop
    path = os.path.join("/Users/eliapizzati/projects/outflows/data_mcmc", "{}.h5".format(filename))
elif stored_data_loc == "github":
    folder = "data_emcee"
    if not os.path.exists(os.path.join(my_dir.data_dir, folder)):
        os.mkdir(os.path.join(my_dir.data_dir, folder))
    path = os.path.join(my_dir.data_dir, folder, "{}.h5".format(filename))
else:
    raise ValueError("stored_data_loc not recognized")


folder_plot = "plot_emcee"
if not os.path.exists(os.path.join(my_dir.plot_dir, folder_plot)):
    os.mkdir(os.path.join(my_dir.plot_dir, folder_plot))


reader = emcee.backends.HDFBackend(path)
samples = reader.get_chain(thin=thin, discard=discard)
samples_flat = reader.get_chain(flat=True, thin=thin, discard=discard)

log_prob_samples = reader.get_log_prob(flat=True, thin=thin, discard=discard)
log_prior_samples = reader.get_blobs(flat=True, thin=thin, discard=discard)

print("flat log prob shape: {0}".format(log_prob_samples.shape))

if print_autocorr_time:
    print("tau = ", reader.get_autocorr_time(thin=1, discard=1000))
    print("average tau = ", np.mean(reader.get_autocorr_time(thin=1, discard=1000)))
    print("N/100 = ", reader.get_chain(thin=1, discard=1000).shape[0]/100)


all_samples = np.concatenate(
    (samples_flat, log_prob_samples[:, None]),
    axis=1)

ndim = 3
labels = ["log beta", "log SFR", "log v_c"]

if plot1:
    """
    PLOT1 - chains
    """

    fig, axes = plt.subplots(ndim, figsize=(8, 7), sharex=True)

    linewidth = 0.8
    alpha = 0.8

    color_tuple = ('black', 'green', 'red', 'cyan', 'magenta', 'blue', 'darkorange', 'yellow', 'dodgerblue', 'purple',
                       'lightgreen', 'cornflowerblue')
    colors = itertools.cycle(color_tuple)
    color_vec = [next(colors) for iwalk in range(nwalkers)]

    for i in range(ndim):
        ax = axes[i]

        for iwalk in range(nwalkers):

            ax.plot(samples[:, iwalk, i], alpha=alpha, color=color_vec[iwalk], linewidth=linewidth)


        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i], labelpad=2)
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number")


    plt.subplots_adjust(left=0.14,  # the left side of the subplots of the figure
                        right=0.97,  # the right side of the subplots of the figure
                        bottom=0.13,  # the bottom of the subplots of the figure
                        top=0.9,  # the top of the subplots of the figure
                        wspace=0.1,  # the amount of width reserved for space between subplots,
                        # expressed as a fraction of the average axis width
                        hspace=0.1)  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height
    fig.suptitle("{0:}, v_c = {1:.1f} km/s, SFR = {2:.1f}".format(data.params_obs["name"], data.params_obs["v_c"],
                                                                  data.params_obs["SFR"]))

    if plot_saving:
        plt.savefig(os.path.join(my_dir.plot_dir, folder_plot, "chain_{}.png".format(filename)))


if plot2:
    """
    PLOT2 - corner
    """

    if plot_logprob:
        labels += ["log prob"]#, "log prior"]


    kwargs = dict(
        bins=50, smooth=0.9, labels=labels,  label_kwargs=dict(labelpad=2),
        # title_kwargs=dict(fontsize=14),
        color='#0072C1',
        truth_color='C3', quantiles=[0.16, 0.5, 0.84],
        levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
        plot_density=False, plot_datapoints=True, fill_contours=True,
        max_n_ticks=3)

    if plot_logprob:
        naxes = ndim + 1
    else:
        naxes = ndim

    fig, ax = plt.subplots(naxes, naxes, figsize=(10, 10))

    if plot_logprob:
        corner.corner(np.nan_to_num(all_samples, neginf = -200), **kwargs, fig=fig)
    else:
        corner.corner(np.nan_to_num(samples_flat, neginf=-200), **kwargs, fig=fig)

    axes = np.array(fig.axes).reshape((naxes, naxes))

    # Loop over the diagonal

    if plot_logprob:
        true_values = [None, data.params_obs["log_SFR"], data.params_obs["log_v_c"], None]

        err_ups = [None, data.params_obs["log_SFR"] + data.params_obs["log_SFR_err_up"],
                   data.params_obs["log_v_c"] + data.params_obs["log_v_c_err_up"], None]
        err_downs = [None, data.params_obs["log_SFR"] - data.params_obs["log_SFR_err_down"],
                     data.params_obs["log_v_c"] - data.params_obs["log_v_c_err_down"], None]

    else:
        true_values = [None, data.params_obs["log_SFR"], data.params_obs["log_v_c"]]

        err_ups = [None, data.params_obs["log_SFR"] + data.params_obs["log_SFR_err_up"],
                   data.params_obs["log_v_c"] + data.params_obs["log_v_c_err_up"]]
        err_downs = [None, data.params_obs["log_SFR"] - data.params_obs["log_SFR_err_down"],
                     data.params_obs["log_v_c"] - data.params_obs["log_v_c_err_down"]]

    # Loop over the diagonal
    for i in range(naxes):

        ax = axes[i, i]
        if i == 1 or i == 2:
            ax.axvline(true_values[i], color=kwargs["truth_color"], linestyle='--', linewidth=1.3)
            ax.axvspan(err_downs[i], err_ups[i], alpha=0.2, color='red')

        if i == 3:
            ax.set_xlim(-500, 100)

        summary = namedtuple('summary', ['median', 'lower', 'upper', 'string'])
        quantiles = kwargs['quantiles']
        quants_to_compute = np.array([quantiles[0], 0.5, quantiles[2]])
        quants = np.percentile(samples_flat.T[i], quants_to_compute * 100)
        summary.median = quants[1]
        summary.plus = quants[2]
        summary.minus = quants[0]

        if i == 0:
            ax.set_title(r"$ %.2f_{-%.2f}^{+%.2f} $ " % (10**summary.median, 10**summary.median-10**summary.minus, 10**summary.plus-10**summary.median))

        if i == 1:
            ax.set_title(r"$ %.2f_{-%.2f}^{+%.2f}" %  (10**summary.median, 10**summary.median-10**summary.minus, 10**summary.plus-10**summary.median))


        if i == 2:
            ax.set_title(
                r"$ %.2f_{-%.2f}^{+%.2f} " % (10**summary.median, 10**summary.median-10**summary.minus, 10**summary.plus-10**summary.median))

    # Loop over the histograms
    for yi in range(naxes):
        for xi in range(yi):

            ax = axes[yi, xi]

            if yi == 3:
                ax.set_ylim(-500, 200)

            if xi == 1 or xi == 2:
                ax.axvline(true_values[xi], color=kwargs["truth_color"], linestyle='--', linewidth=1.3)
                ax.axvspan(err_downs[xi], err_ups[xi], alpha=0.2, color='red')
                if xi == 1:
                    ax.set_xlim(0., 2.4)
                if xi == 2:
                    ax.set_xlim(2., 2.6)

            if yi == 1 or yi == 2:
                ax.axhline(true_values[yi], color=kwargs["truth_color"], linestyle='--', linewidth=1.3)
                ax.axhspan(err_downs[yi], err_ups[yi], alpha=0.2, color='red')
                if yi == 1:
                    ax.set_ylim(0., 2.4)
                if yi == 2:
                    ax.set_ylim(2., 2.6)

    #fig.suptitle("{0:}, v_c = {1:.1f} km/s, SFR = {2:.1f}".format(data.params_obs["name"], data.params_obs["v_c"],
    #                                                              data.params_obs["SFR"]))

    if plot_saving:
        plt.savefig(os.path.join(my_dir.plot_dir, folder_plot, "corner_{}.png".format(filename)))



if plot3:
    """
    PLOT3 - emission from samples
    """

    from mcmc_fast_emission import get_emission_fast, get_other_params
    from mcmc_main import h, grid

    number_of_lines_to_plot = 100

    other_params = get_other_params(data.params_obs["redshift"], data.params_obs["line_FWHM"])

    beam_interp = np.interp(h, data.x_beam / 1e3 / nc.pc, data.beam, right=0.)

    beam_interp[beam_interp < 0.] = 0.

    beam_func = interp1d(h, beam_interp, \
                         fill_value=(beam_interp[0], 0.), bounds_error=False)

    beam_2d = beam_func(np.sqrt(grid[0] ** 2 + grid[1] ** 2))
    f_beam = np.fft.fft2(beam_2d)

    fig_int_conv, ax_int_conv = pltc.plot_configurator(plot_type="int", xlim=15)

    ax_int_conv.set_ylim((1e-3, 1e2))

    counter = 0
    intx = np.random.randint(len(samples_flat), size=number_of_lines_to_plot)

    luminosities = []

    for theta in samples_flat[intx]:
        counter += 1
        print("computing emission for theta =", theta, "\t", "iteration number {}/{}" \
              .format(counter, number_of_lines_to_plot))

        intensity, sigma = get_emission_fast(theta, data, other_params, h, grid, f_beam, return_quantities="emission")

        h_ext = np.linspace(0, np.max(h), 100)
        sigma_ext = np.interp(h_ext, h, sigma)
        luminosity_cii = np.trapz(sigma_ext * 2 * np.pi * h_ext * nc.pc * 1e3, h_ext * nc.pc * 1e3)

        ax_int_conv.plot(h, intensity, alpha=0.1, color="gray")

        luminosities.append(luminosity_cii / nc.ls)

    luminosities = np.asarray(luminosities)
    # print(luminosities)
    print("model CII luminosities",np.percentile(luminosities/1e9, [18,50,84])[1], np.percentile(luminosities/1e9, [18,50,84])[1]-np.percentile(luminosities/1e9, [18,50,84])[0],
          np.percentile(luminosities/1e9, [18,50,84])[2]-np.percentile(luminosities/1e9, [18,50,84])[1])
    print("data CII luminosities", data.params_obs["luminosity_CII"]/1e9,  data.params_obs["luminosity_CII_err_up"]/1e9, data.params_obs["luminosity_CII_err_down"]/1e9)


    alpine = ax_int_conv.errorbar(data.x / (1000 * nc.pc), data.data, yerr=data.err, \
                                  markerfacecolor='maroon', markeredgecolor='maroon', marker='o', \
                                  linestyle='', ecolor='maroon')

    fig_int_conv.suptitle(
        "{0:}, v_c = {1:.1f} km/s, SFR = {2:.1f}".format(data.params_obs["name"], data.params_obs["v_c"],
                                                         data.params_obs["SFR"]))

    fig_int_conv.legend(loc="lower center", ncol=8, fontsize="small")

    if plot_saving:
        plt.savefig(os.path.join(my_dir.plot_dir, folder_plot, "emission_{}.png".format(filename)))





plt.show()


