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

import mydir


folder = "data_emcee"
    
if not os.path.exists(os.path.join(mydir.data_dir, folder)):
    os.mkdir(os.path.join(mydir.data_dir, folder))

filename = os.path.join(mydir.data_dir, folder, "trial_run.h5")

reader = emcee.backends.HDFBackend(filename)

samples = reader.get_chain(flat=True)

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
    (samples, log_prob_samples[:, None]),\
     #log_prior_samples[:, None]), \
     axis=1)

ndim = 3

labels = ["beta", "SFR", "v_c"]
labels += ["log prob"]# "log prior"]

corner.corner(all_samples, labels=labels);
