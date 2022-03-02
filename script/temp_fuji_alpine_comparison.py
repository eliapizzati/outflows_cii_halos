


import os
import numpy as np
import mydir
import natconst as nc


import matplotlib.pyplot as plt
import matplotlib


from model_classes import load_from_file
from load_data import obs_data_list, names, names_CII_halo, names_wo_CII_halo,\
                      names_other,  observational_data_fuji



matplotlib.rcParams.update({
        "font.size": 16.0,
#        "font.family": 'sans-serif',
#        "font.sans-serif": ['Helvetica'],
        "axes.titlesize": 16.,
        "axes.labelsize": 16.,
        "xtick.labelsize": 16.,
        "ytick.labelsize": 16.,
        "xtick.major.size": 10.0,
        "ytick.major.size": 10.0,
        "xtick.minor.size": 6.0,
        "ytick.minor.size": 6.0,
 #       "legend.fontsize": 12.0,
 #       "figure.dpi": 200,
 #       "savefig.dpi": 200,
#        "text.usetex": True
})


size = 14


datas = obs_data_list

data_counter = 0


data_container_name = "CII_halo_NFW"

fig, ax = plt.subplots(figsize = (7,6))


ax.set_yscale("log")

ax.set_ylim(1e-3, 1e2)

ax.set_ylabel("flux [mJy/arcsec$^2$]", size=size)
ax.set_xlabel("b [kpc]", size=size)

#x_axis =
#data_tot = []

x_axis = np.linspace(0.,15.,100)

y = np.zeros_like(x_axis)

count_data = 0

for data in datas:

    #if data.params_obs["name"] not in names_CII_halo:
    #    pass
    #else:

        y += np.interp(x_axis,  data.x / (1000 * nc.pc), data.data)
        count_data += 1

y /=  count_data

alpine, = ax.plot(x_axis, y, color='maroon')

fuji = ax.errorbar(observational_data_fuji.x / (1000 * nc.pc), observational_data_fuji.data,
                            yerr=observational_data_fuji.err, \
                            markerfacecolor='navy', markeredgecolor='navy', marker='d', \
                            linestyle='', ecolor='navy')

ax.legend([alpine, fuji], ["ALPINE CII halo avg", "Fujimoto+19"])  ##

plt.subplots_adjust(left=0.15,  # the left side of the subplots of the figure
                    right=0.97,  # the right side of the subplots of the figure
                    bottom=0.15,  # the bottom of the subplots of the figure
                    top=0.95,  # the top of the subplots of the figure
                    wspace=0.,  # the amount of width reserved for space between subplots,
                    # expressed as a fraction of the average axis width
                    hspace=0.1)  # the amount of height reserved for space between subplots,
# expressed as a fraction of the average axis height


plt.show()

