# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 19:08:15 2021

@author: anna
"""

import os
import numpy as np
import matplotlib.pyplot as plt

import mydir

from load_data import names_wo_CII_halo_short


data_container_name = "wo_CII_halo.npy"


out_filename = os.path.join(mydir.data_dir, "data_chi2", data_container_name)

chi2_names = np.load(out_filename)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 4))


ax.imshow(chi2_names)

ax.set_yticks([y  for y in range(len(names_wo_CII_halo_short))])
ax.set_xlabel('Four separate samples')
ax.set_ylabel('betas')

# add x-tick labels
plt.setp(ax, yticks=[y for y in range(len(names_wo_CII_halo_short))],
         yticklabels=names_wo_CII_halo_short)

for tick in ax.get_xticklabels():
    tick.set_rotation(90)