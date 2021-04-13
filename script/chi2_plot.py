# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 19:08:15 2021

@author: anna
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors
 

import mydir

from load_data import names_wo_CII_halo_short, names_CII_halo_short, names_other_short, obs_data_list

names_plot = names_other_short

betas_plot = np.asarray([1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6])

data_container_name = "other"


out_filename = os.path.join(mydir.data_dir, "data_chi2", "{}.npy".format(data_container_name))

chi2_names = np.load(out_filename)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(1.8*8.27,1.5*4.))


im = ax.imshow(chi2_names, cmap=cm.viridis,norm=colors.PowerNorm(gamma=0.2), interpolation='nearest',\
               origin='lower', aspect='auto') #extent=extent)

# We want to show all ticks...
ax.set_xticks(np.arange(len(betas_plot)))
ax.set_yticks(np.arange(len(names_plot)))
# ... and label them with the respective list entries
ax.set_xticklabels(betas_plot)
ax.set_yticklabels(names_plot)


# Create colorbar
cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.set_ylabel("likelihood", rotation=-90, va="bottom")


# Turn spines off and create white grid.

for string in ["left", "right", "bottom","top"]:
    ax.spines[string].set_visible(False)
    
ax.set_xticks(np.arange(len(betas_plot)+1)-.5, minor=True)
ax.set_yticks(np.arange(len(names_plot)+1)-.5, minor=True)
ax.grid(which="minor", color="w", linestyle='-', linewidth=2,axis='y')
ax.tick_params(which="minor", bottom=False, left=False)


#
ax.set_xlabel('betas')
ax.set_ylabel('')




plt.subplots_adjust(left = 0.15,  # the left side of the subplots of the figure
    right = 0.98,   # the right side of the subplots of the figure
    bottom = 0.15,  # the bottom of the subplots of the figure
    top = 0.9,     # the top of the subplots of the figure
    wspace = 0.1,  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    hspace = 0.1)  # the amount of height reserved for space between subplots,
                  # expressed as a fraction of the average axis height

ax.set_title(data_container_name)




fig, ax = plt.subplots(figsize=(8.27,5.))

cmap_rend_col = cm.get_cmap('viridis_r')

v_c_min = 140.
v_c_max = 250.

norm = colors.Normalize(vmin=v_c_min, vmax=v_c_max)

    
names_short_fit = ["DC_881725", "DC_630594", "DC_488399","DC_683613","vc_5110377875","DC_880016",\
                   "DC_733857","DC_709575", "DC_539609",\
                   "DC_494057", "ve_530029038","vc_5100994794", "vc_5100969402", "DC_834764"]

betas_fit = [2.1,3.4,3.4,2.3,2.0,3.2,\
             2.9,3.1,3.6,\
             2.9,1.9,3.6,3.6,3.2]
  
for data in obs_data_list:
    
    for name, beta in zip(names_short_fit, betas_fit):
        if name == data.params_obs["name_short"]:
            print(name, beta)
            data.params_obs.update(beta_best_fit=beta)

    if data.params_obs["halo_class"] == "CII_halo":
        marker = "o"
        cii = ax.scatter(data.params_obs["SFR"], data.params_obs["beta_best_fit"],\
               marker=marker, color = cmap_rend_col((data.params_obs["v_c"]-v_c_min)/(v_c_max-v_c_min)),\
               s=50)

    elif data.params_obs["halo_class"] == "wo_CII_halo":
        marker = "D"
        wo_cii = ax.scatter(data.params_obs["SFR"], data.params_obs["beta_best_fit"],\
               marker=marker, color = cmap_rend_col((data.params_obs["v_c"]-v_c_min)/(v_c_max-v_c_min)),\
               s=50)

    elif data.params_obs["halo_class"] == "other":
        marker = "s"
        oth = ax.scatter(data.params_obs["SFR"], data.params_obs["beta_best_fit"],\
               marker=marker, color = cmap_rend_col((data.params_obs["v_c"]-v_c_min)/(v_c_max-v_c_min)),\
               s=70)
       


ax.set_ylim(1.3,3.7)
ax.set_xlabel("SFR")
ax.set_ylabel("beta")

plt.subplots_adjust(left = 0.15,  # the left side of the subplots of the figure
    right = 0.98,   # the right side of the subplots of the figure
    bottom = 0.15,  # the bottom of the subplots of the figure
    top = 0.9,     # the top of the subplots of the figure
    wspace = 0.1,  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    hspace = 0.1)  # the amount of height reserved for space between subplots,
                  # expressed as a fraction of the average axis height


cmap = cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
cmap.set_array([])

cb = fig.colorbar(cmap, orientation='vertical')
cb.set_label(r'$v_c$', rotation=0.)


plt.legend([cii, wo_cii, oth], ["CII", "wo CII", "other"])##

ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis

ax2.hist([data.params_obs["SFR"] for data in obs_data_list], alpha=0.2, color='g', bins=500, density=False)
ax2.set_yticks([])
ax2.set_ylim(0.,1.)



