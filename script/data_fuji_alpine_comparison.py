import numpy as np
import natconst as nc


import matplotlib.pyplot as plt
import matplotlib

from load_data import obs_data_list, observational_data_fuji



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



fig, ax = plt.subplots(figsize = (7,6))


ax.set_yscale("log")

ax.set_ylim(1e-3, 1e2)

ax.set_ylabel("flux [mJy/arcsec$^2$]", size=size)
ax.set_xlabel("b [kpc]", size=size)

data = datas[1]
alpine = ax.errorbar(data.x / (1000 * nc.pc), data.data,
                            yerr=data.err, \
                            markerfacecolor='navy', markeredgecolor='navy', marker='d', \
                            linestyle='', ecolor='navy', label="Data from Seiji")

print(data.params_obs["name_short"])


#x_mas = np.asarray([0.1471169994270768, 0.45771480203273307, 0.7640018153279446,1.058074158768218,1.3580970451490058,1.6636289628455216,1.961920729382407,2.2669263421365904])
#y =  np.asarray([0.9535182095209478, 0.7580882071309459, 0.43701675214435326, 0.21444930153684946,  0.07282819740825795 , 0.04594872925847865, 0.012911312372689132, 0.01895828537743512])

x_mas = np.asarray([0.16147798742138353, 0.46525157232704395, 0.7617924528301886,1.0583333333333331,1.362106918238994,1.6586477987421384,1.9588050314465408,2.262578616352202])
y =  np.asarray([1.0824133522788788, 0.8314527714489289, 0.48578475070320354, 0.22338456451964553,  0.077239758523539 , 0.00979145507086018, 0.017094666840436493, 0.021441311325220903])
y_up =  np.asarray([1.0824133522788788, 0.8314527714489289, 0.48578475070320354, 0.2635485430338834,  0.1099917981852589 , 0.04026437119529853, 0.041132301336074845, 0.03521062380461966])
y_down = np.asarray([1.0824133522788788, 0.8314527714489289, 0.48578475070320354, 0.19705095177516463,  0.050653339307195004 , 0.0010902784377834506, 0.001532717760955783, 0.007382730909593146])



from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

x =    x_mas / cosmo.arcsec_per_kpc_proper(data.params_obs["redshift"]).value

ax.errorbar(x, y,
                yerr=[y-y_down, y_up-y], \
                markerfacecolor='maroon', markeredgecolor='maroon', marker='o', \
             linestyle='', ecolor='maroon', label="Fujimoto+20 paper")

ax.legend()
fig.tight_layout()
plt.show()

