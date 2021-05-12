# -*- coding: utf-8 -*-
"""
Created on Wed May 12 14:35:50 2021

@author: anna
"""

import matplotlib.pyplot as plt

import numpy as np
import os
import time
import itertools

import mydir
from trial_emcee import data, other_params, get_emission_fast
import plot_config as pltc
import natconst as nc


betas = [2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5]
SFRs = [10.,50.,100.]
v_cs = [150.,200.,250.,300.]

for  beta, SFR, v_c in itertools.product(betas, SFRs, v_cs):

    theta = [beta, SFR,v_c]
    print(theta)
    fig_int_conv, ax_int_conv = pltc.plot_configurator(plot_type="int", xlim=15)
    
    ax_int_conv.set_ylim((1e-3,1e2))
     
        
    folder = "plot_fast_emission"
        
    if not os.path.exists(os.path.join(mydir.plot_dir, folder)):
        os.mkdir(os.path.join(mydir.plot_dir, folder))
        
    time_profile = time.perf_counter()
        
    h, intensity = get_emission_fast(theta, data, other_params)    
            
    time_profile = (time.perf_counter() - time_profile)
        
    print("total emission time (s)=", time_profile)
        
    ax_int_conv.plot(h, intensity)

alpine = ax_int_conv.errorbar(data.x/(1000*nc.pc), data.data, yerr=data.err, \
        markerfacecolor='maroon',markeredgecolor='maroon', marker='o',\
        linestyle='', ecolor = 'maroon')
                        

   
fig_int_conv.legend(loc="lower center", ncol=8, fontsize="small")
plt.savefig(os.path.join(mydir.plot_dir, folder, "emission.png"))
    
        
        
    
    
    
