# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 17:16:03 2021

@author: anna
"""



import os
#import numpy as np
import matplotlib.pyplot as plt
import mydir

from solver_profiles import profile, solve_eq


params = dict([("type", "n")
               ("SFR", 20.),
               ("beta", 1.0), 
               ("f_esc", 0.), 
               ("v_c", 175.),
               ("Zeta", 1.0),
               ("alfa", 1.0),
               ("R_in", 0.3)])

r, profiles = solve_eq(beta=params["beta"], SFR_pure=params["SFR"], v_c_pure=params["v_c"])
    
sample = profile(radius=r, variable=profiles[0], params=params)

sample.to_file()

print(sample.check_nans())

fig, ax = plt.subplots()
sample.plot(ax=ax) 

fig.savefig(os.path.join(mydir.plot_dir, "trial.png"))
    




