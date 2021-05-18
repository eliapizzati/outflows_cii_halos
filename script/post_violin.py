# -*- coding: utf-8 -*-
"""
Created on Mon May 17 14:49:56 2021

@author: anna
"""



import os
import matplotlib.pyplot as plt
import numpy as np
import emcee

import plot_config as pltc
import mydir
import natconst as nc

from load_data import obs_data_list, names, names_CII_halo, names_wo_CII_halo,\
                      names_other,  observational_data_fuji, names_CII_halo_short


names_list = ["DC_396844", "DC_683613", "DC_881725"]#, "DC_488399"]#, "DC_630594", "DC_880016"]

datas = []
names_plot = []

for data in obs_data_list:

    if data.params_obs["name_short"] in names_list:
        datas.append(data)
        names_plot.append(data.params_obs["name_short"])


folder_data = "data_emcee"
    
if not os.path.exists(os.path.join(mydir.data_dir, folder_data)):
    os.mkdir(os.path.join(mydir.data_dir, folder_data))


folder_plot = "plot_emcee"

if not os.path.exists(os.path.join(mydir.plot_dir, folder_plot)):
    os.mkdir(os.path.join(mydir.plot_dir, folder_plot))


nwalkers = 96
nsteps = 1e3

betas = []
logprobs = []

for data in datas:
    filename = "{}_{:.0f}".format(data.params_obs["name_short"], nsteps)

    print("###################################################################")
    
    print("postprocessing an MCMC with the following params:")
    print("n steps = {}".format( nsteps))
    print("n walkers = {}".format( nwalkers))
    print("data object = {}".format( data.params_obs["name_short"]))
    print("filename = {}".format( filename))
                 
    print("###################################################################")


    path = os.path.join(mydir.data_dir, folder_data, "{}.h5".format(filename))
    
    reader = emcee.backends.HDFBackend(path)
    
    samples_flat = reader.get_chain(flat=True)
    
    log_prob_samples = reader.get_log_prob(flat=True)

    all_samples = np.concatenate(
            (samples_flat, log_prob_samples[:, None]),\
            #log_prior_samples[:, None]), \
            axis=1)

    betas.append(all_samples[:,0])
    logprobs.append(all_samples[:,3][all_samples[:,3]>-20])

betas = np.asarray(betas)
logprobs = np.asarray(logprobs)

ndim = 3


fig, [ax_betas, ax_probs] = plt.subplots(1,2, figsize=(13,7), sharey=True)

parts = ax_betas.violinplot(betas, positions=None, vert=False, widths=0.8, showmeans=False, showextrema=False, showmedians=False,\
              #quantiles=[[0.25,0.75]]*len(betas), \
               points=100, bw_method=None)

i=0
for pc in parts['bodies']:
    pc.set_facecolor('C{}'.format(i))
    pc.set_edgecolor('gray')
    pc.set_alpha(1)
    i=i+1


parts = ax_probs.violinplot(logprobs, positions=None, vert=False, widths=0.8, showmeans=False, showextrema=False, showmedians=False,\
              #quantiles=[[0.25,0.75]]*len(betas), \
               points=100, bw_method=None)

i=0
for pc in parts['bodies']:
    pc.set_facecolor('C{}'.format(i))
    pc.set_edgecolor('gray')
    pc.set_alpha(1)
    i=i+1
    
ax_betas.set_xlabel("beta")

ax_probs.set_xlim(-10,1)
ax_probs.set_xlabel("log probability")

ax_betas.set_yticks(np.arange(1,len(names_plot)+1))
ax_betas.set_yticklabels(names_plot)
