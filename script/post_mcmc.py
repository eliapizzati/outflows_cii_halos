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

import plot_config as pltc
import mydir
import natconst as nc

from load_data import obs_data_list, names, names_CII_halo, names_wo_CII_halo,\
                      names_other,  observational_data_fuji

"""
1:  DC396844
8:  DC683613
13: DC881725
17: vc...875
3:  DC488399
7:  DC630594
12: DC880016
14: vc...582
"""            

log = True

emission = True

nwalkers= 96
nsteps = 1e4

sample_step = int(40 * (nsteps/1e4))
walker_step = int(12 * (nwalkers/96))

int_data = int(input("data number?"))
data = obs_data_list[int_data]

filename = "{}_{:.0f}".format(data.params_obs["name_short"], nsteps)

if log:
    filename += "_new_priors"

print("###################################################################")

print("postprocessing an MCMC with the following params:")
print("n steps = {}".format( nsteps))
print("n walkers = {}".format( nwalkers))
print("data object = {}".format( data.params_obs["name_short"]))
print("filename = {}".format( filename))
print("beta log =", log)
             
print("###################################################################")



folder_data = "data_emcee"
    
if not os.path.exists(os.path.join(mydir.data_dir, folder_data)):
    os.mkdir(os.path.join(mydir.data_dir, folder_data))


folder_plot = "plot_emcee"

if not os.path.exists(os.path.join(mydir.plot_dir, folder_plot)):
    os.mkdir(os.path.join(mydir.plot_dir, folder_plot))


path = os.path.join(mydir.data_dir, folder_data, "{}.h5".format(filename))

reader = emcee.backends.HDFBackend(path)

samples = reader.get_chain()
samples_flat = reader.get_chain(flat=True)

ndim = 3

fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)

if log:
    labels = ["log beta", "SFR", "v_c"]
else:
    labels = ["beta", "SFR", "v_c"]

for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)
    
axes[-1].set_xlabel("step number")

fig.suptitle("{0:}, v_c = {1:.1f} km/s, SFR = {2:.1f}".format(data.params_obs["name"], data.params_obs["v_c"], data.params_obs["SFR"]))

plt.savefig(os.path.join(mydir.plot_dir, folder_plot, "chain_{}.png".format(filename)))


#tau = reader.get_autocorr_time()
#burnin = int(2 * np.max(tau))
#thin = int(0.5 * np.min(tau))
log_prob_samples = reader.get_log_prob(flat=True)
#log_prior_samples = reader.get_blobs(flat=True)

#print("burn-in: {0}".format(burnin))
#print("thin: {0}".format(thin))
#print("flat chain shape: {0}".format(samples.shape))
print("flat log prob shape: {0}".format(log_prob_samples.shape))
#print("flat log prior shape: {0}".format(log_prior_samples.shape))


all_samples = np.concatenate(
    (samples_flat, log_prob_samples[:, None]),\
     #log_prior_samples[:, None]), \
     axis=1)

ndim = 3

if log:
    labels = ["log beta", "SFR", "v_c"]
else:
    labels = ["beta", "SFR", "v_c"]
#labels += ["log prob"]# "log prior"]


kwargs = dict(
            bins=50, smooth=0.9, labels=labels, #label_kwargs=dict(fontsize=14),
            #title_kwargs=dict(fontsize=14), 
            color='#0072C1',
            truth_color='C3', quantiles=[0.16, 0.5, 0.84],
            levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
            plot_density=False, plot_datapoints=True, fill_contours=True,
            max_n_ticks=3)

naxes = ndim 

fig,ax = plt.subplots(naxes,naxes,figsize=(10,10))

corner.corner(samples_flat, **kwargs, fig=fig)

axes = np.array(fig.axes).reshape((naxes, naxes))


true_values = [None, data.params_obs["SFR"], data.params_obs["v_c"], None]

err_ups = [None, data.params_obs["SFR"]+data.params_obs["SFR_err_up"], data.params_obs["v_c"]+data.params_obs["v_c_err_up"], None ]
err_downs = [None, data.params_obs["SFR"]-data.params_obs["SFR_err_down"], data.params_obs["v_c"]-data.params_obs["v_c_err_down"], None ]

# Loop over the diagonal
for i in range(naxes):
    
    ax = axes[i, i]
    if i == 1 or i == 2:
        ax.axvline(true_values[i], color=kwargs["truth_color"], linestyle='--', linewidth=1.3)
        ax.axvspan(err_downs[i], err_ups[i], alpha=0.2, color='red')

    if i == 3:
        ax.set_xlim(-500,100)

    
# Loop over the histograms
for yi in range(naxes):
    for xi in range(yi):
        
        
        ax = axes[yi, xi]
        
        if yi == 3:
            ax.set_ylim(-500,200)
        
        
        if xi == 1 or xi == 2:
            ax.axvline(true_values[xi], color=kwargs["truth_color"], linestyle='--', linewidth=1.3)
            ax.axvspan(err_downs[xi], err_ups[xi], alpha=0.2, color='red')
            if xi == 1:
                ax.set_xlim(1.,250.)
            if xi == 2:
                ax.set_xlim(100.,450.)


        if yi == 1 or yi == 2:
            ax.axhline(true_values[yi], color=kwargs["truth_color"], linestyle='--', linewidth=1.3)
            ax.axhspan(err_downs[yi], err_ups[yi], alpha=0.2, color='red')
            if yi == 1:
                ax.set_ylim(1.,250.)
            if yi == 2:
                ax.set_ylim(100.,450.)


        
fig.suptitle("{0:}, v_c = {1:.1f} km/s, SFR = {2:.1f}".format(data.params_obs["name"], data.params_obs["v_c"], data.params_obs["SFR"]))

plt.savefig(os.path.join(mydir.plot_dir, folder_plot, "corner_{}.png".format(filename)))

# emission

if emission:
    from mcmc import get_emission_fast, get_other_params
    from mcmc import h, grid
        
    
    other_params = get_other_params(data.params_obs["redshift"], data.params_obs["line_FWHM"])
    
    beam_interp = np.interp(h, data.x_beam/1e3/nc.pc, data.beam, right=0.)
    
    beam_interp[beam_interp<0.] = 0.
    
    beam_func = interp1d(h, beam_interp, \
                         fill_value = (beam_interp[0], 0.), bounds_error=False)
    
    beam_2d = beam_func(np.sqrt(grid[0]**2 + grid[1]**2))
    f_beam = np.fft.fft2(beam_2d)

        
    fig_int_conv, ax_int_conv = pltc.plot_configurator(plot_type="int", xlim=15)
    
    ax_int_conv.set_ylim((1e-3,1e2))
    
    counter = 0
    for  walker in samples[::sample_step]:
        for theta in walker[::walker_step]:
            counter += 1
            print("computing emission for theta =", theta, "\t", "iteration number {}/{}"\
                  .format(counter, int( nsteps*nwalkers/sample_step/walker_step)))

            #if theta[2]<150. and theta[1]>50.:

            intensity = get_emission_fast(theta, data, other_params, h, grid, f_beam)

            ax_int_conv.plot(h, intensity, alpha=0.1, color="gray")

        
    alpine = ax_int_conv.errorbar(data.x/(1000*nc.pc), data.data, yerr=data.err, \
            markerfacecolor='maroon',markeredgecolor='maroon', marker='o',\
            linestyle='', ecolor = 'maroon')
    
    fig_int_conv.suptitle("{0:}, v_c = {1:.1f} km/s, SFR = {2:.1f}".format(data.params_obs["name"], data.params_obs["v_c"], data.params_obs["SFR"]))

    fig_int_conv.legend(loc="lower center", ncol=8, fontsize="small")
    plt.savefig(os.path.join(mydir.plot_dir, folder_plot, "emission_{}.png".format(filename)))
        
            
    plt.show()
        
    
