# -*- coding: utf-8 -*-
"""
Created on Thu May 13 14:21:16 2021

@author: anna
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import corner
import emcee

import plot_config as pltc
import mydir
import natconst as nc

from load_data import obs_data_list, names, names_CII_halo, names_wo_CII_halo,\
                      names_other,  observational_data_fuji

             


emission = False 

nwalkers= 96
nsteps = 1e3

data = obs_data_list[8]



folder_data = "data_emcee"
    
if not os.path.exists(os.path.join(mydir.data_dir, folder_data)):
    os.mkdir(os.path.join(mydir.data_dir, folder_data))


folder_plot = "plot_emcee"

if not os.path.exists(os.path.join(mydir.plot_dir, folder_plot)):
    os.mkdir(os.path.join(mydir.plot_dir, folder_plot))

filename = "{}_{:.0f}".format(data.params_obs["name_short"], nsteps)

path = os.path.join(mydir.data_dir, folder_data, "{}.h5".format(filename))

reader = emcee.backends.HDFBackend(path)

samples = reader.get_chain()
samples_flat = reader.get_chain(flat=True)

ndim = 3

fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
labels = ["beta", "SFR", "v_c"]

for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)
    
axes[-1].set_xlabel("step number")

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

labels = ["beta", "SFR", "v_c"]
labels += ["log prob"]# "log prior"]

corner.corner(all_samples, labels=labels)
plt.savefig(os.path.join(mydir.plot_dir, folder_plot, "corner_{}.png".format(filename)))


# emission

if emission:
    from mcmc import get_emission_fast
    from mcmc import other_params, h, f_beam, grid

        
    fig_int_conv, ax_int_conv = pltc.plot_configurator(plot_type="int", xlim=15)
    
    ax_int_conv.set_ylim((1e-3,1e2))
    
        
    for  walker in samples[::20]:
        for theta in walker[::8]:
    
            print(theta)
            
            intensity = get_emission_fast(theta, data, other_params, h, grid, f_beam)    
            
            ax_int_conv.plot(h, intensity, alpha=0.1, color="gray")
        
        
    alpine = ax_int_conv.errorbar(data.x/(1000*nc.pc), data.data, yerr=data.err, \
            markerfacecolor='maroon',markeredgecolor='maroon', marker='o',\
            linestyle='', ecolor = 'maroon')
    
       
    fig_int_conv.legend(loc="lower center", ncol=8, fontsize="small")
    plt.savefig(os.path.join(mydir.plot_dir, folder_plot, "emission_{}.png".format(filename)))
        
            
            
        
    
