# -*- coding: utf-8 -*-
"""
Created on Mon May 17 14:49:56 2021

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


names_list = ["DC_396844, DC_]
datas = []

for data in obs_data_list:

    if data.params_obs["name_short"] in names_list:
        datas.append(data)


emission = True 

nwalkers= 96
nsteps = 1e3

sample_step = int(20 * (nsteps/1e3))
walker_step = int(12 * (nwalkers/96))

data = obs_data_list[17]

filename = "{}_{:.0f}".format(data.params_obs["name_short"], nsteps)

print("###################################################################")

print("postprocessing an MCMC with the following params:")
print("n steps = {}".format( nsteps))
print("n walkers = {}".format( nwalkers))
print("data object = {}".format( data.params_obs["name_short"]))
print("filename = {}".format( filename))
             
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
