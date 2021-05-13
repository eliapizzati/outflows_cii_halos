# -*- coding: utf-8 -*-
"""
Created on Thu May 13 12:32:57 2021

@author: anna
"""


import numpy as np

import matplotlib.pyplot as plt


from load_data import obs_data_list, names, names_CII_halo, names_wo_CII_halo,\
                      names_other,  observational_data_fuji

data = obs_data_list[1]


def log_prior_gaussian(theta, data):
    """
    defines the priors for the set of parameters theta
    
    Parameters
    ==========
    theta: array
            
    Returns
    =======
    priors value: float

    """    
    
    beta, SFR, v_c = theta
    
    if 1.0 < beta < 8.0 and SFR >= 1. and v_c >= 50.:
        prior =  0.0
        
        prior += - 2*(SFR-data.params_obs["SFR"])**2/(data.params_obs["SFR_err_up"]+data.params_obs["SFR_err_down"])**2
        prior += - 2*(v_c-data.params_obs["v_c"])**2/(data.params_obs["v_c_err_up"]+data.params_obs["v_c_err_down"])**2

        return prior
    else:
        return -np.inf

betas = np.linspace(0.,9., 1000)
SFRs = np.linspace(10.,300.,1000)
vcs = np.linspace(10.,300.,1000)

prior_beta = np.asarray([log_prior_gaussian([beta, 50., 250.], data)  for beta in betas])
prior_SFR = np.asarray([log_prior_gaussian([3.0, SFR, 200.], data)  for SFR in SFRs])
prior_vc = np.asarray([log_prior_gaussian([3.0, 50., vc], data)  for vc in vcs])

fig,ax = plt.subplots()

#ax.plot(betas,10**prior_beta)
#ax.plot(SFRs,10**prior_SFR)
ax.plot(vcs,10**prior_vc)

