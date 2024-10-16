

"""
This script fenerates the plot for the paper for the non detection analysis
N.B. Fits of the eta-sfr and eta-mstar relations need to be updated!
"""

import numpy as np
import os
import matplotlib.pyplot as plt


from load_data import obs_data_list, names, names_CII_halo, observational_data_fuji
from scipy.interpolate import interp1d


import plot_config as pltc
import my_dir
import my_utils
import natconst as nc

import matplotlib
from matplotlib.patches import Rectangle



matplotlib.rcParams.update({
        "font.size": 14.0,
        "font.family": 'sans-serif',
#        "font.sans-serif": ['Helvetica'],
        "axes.titlesize": 14.,
        "axes.labelsize": 14.,
        "xtick.labelsize": 14.,
        "ytick.labelsize": 14.,
        "xtick.major.size": 4.0,
        "ytick.major.size": 4.0,
        "xtick.minor.size": 2.0,
        "ytick.minor.size": 2.0,
        "legend.fontsize": 13.0,
        "legend.frameon": False
 #       "figure.dpi": 200,
 #       "savefig.dpi": 200,
#        "text.usetex": True
})


data_saving = False
data_loading = True


plot1 = True  # dependencies

folder = "final_dependencies"

path = os.path.join(my_dir.data_dir, folder)

if not os.path.exists(path):
    os.mkdir(path)

int_data = 0#int(input("data number?"))
data = obs_data_list[int_data]

redshift = 4.54

if not data_loading:
    from mcmc_main import get_emission_fast, get_other_params
    from mcmc_main import h, grid

    other_params = get_other_params(data.params_obs["redshift"], data.params_obs["line_FWHM"])

    beam_interp = np.interp(h, data.x_beam / 1e3 / nc.pc, data.beam, right=0.)

    beam_interp[beam_interp < 0.] = 0.

    beam_func = interp1d(h, beam_interp, \
                         fill_value=(beam_interp[0], 0.), bounds_error=False)

    beam_2d = beam_func(np.sqrt(grid[0] ** 2 + grid[1] ** 2))
    f_beam = np.fft.fft2(beam_2d)



if plot1:
    """
    # dependencies
    """
    fig, [ax_sigma, ax_flux] = plt.subplots(2,1,figsize=(1.4*4.2,1.4*5.7), sharex=True)

    # ax_sigma.set_xscale('log')
    ax_sigma.set_yscale('log')
    ax_sigma.set_ylim((1e-6, 5e-2))
    ax_sigma.set_xlim((0.3, 12))
    ax_sigma.tick_params(length=4, axis="both", which="major")
    ax_sigma.tick_params(length=2, axis="both", which="minor")

    # ax_sigma.set_xlabel("b [kpc]")
    ax_sigma.set_ylabel(r"$\Sigma_{\rm CII}$ [erg cm$^{-2}$ s$^{-1}$]")
    # ax.xaxis.set_minor_locator(MultipleLocator(1))

    # ax_flux.set_xscale('log')
    ax_flux.set_yscale('log')
    ax_flux.set_ylim((1e-4, 1e1))
    ax_flux.set_xlim((0.3, 12))
    ax_flux.tick_params(length=4, axis="both", which="major")
    ax_flux.tick_params(length=2, axis="both", which="minor")

    ax_flux.set_xlabel("b [kpc]")
    ax_flux.set_ylabel('flux [mJy/arcsec$^2$]')
    # ax.xaxis.set_minor_locator(MultipleLocator(1))

    # INIT PHASE

    mstars = [1e9,1e10]
    colors = ["limegreen", "teal"]
    lws = [2.5,3]

    for mstar0, color, lw in zip(mstars, colors,lws):
        eta0 = my_utils.from_xfitted_to_eta(mstar0 / 1e10, 4.9, -0.36)
        sfr0 = my_utils.from_eta_to_xfitted(eta0, 14.2, -0.26)

        mvir0 = my_utils.mvir_behroozi(mstar0, z=redshift)

        print("mvir/1e11=",mvir0/1e11)
        print("mstar/1e9=",mstar0/1e9)

        vc0 = my_utils.get_vc_from_virial_mass(mvir0, z=redshift) / 1e5

        params0 = np.asarray([np.log10(eta0), np.log10(sfr0), np.log10(vc0)])

        print(eta0, sfr0, vc0)
        print("\n")

        if data_loading:
            h, sigma0, intensity0 = np.loadtxt(os.path.join(path, f"final_dep_m{np.log10(mstar0):.1f}.dat"))

        else:
            intensity0, sigma0, r0, n0, v0, T0 = get_emission_fast(params0, data, other_params, h, grid, f_beam,
                                                               return_quantities="all")

        qsigma, = ax_sigma.plot(h, sigma0, lw=lw, color=color, linestyle="-",
                                label=r'$\eta$={:.1f}, SFR={:.0f}M$_\odot$yr$^{{-1}}$, v$_c$={:.0f}kms$^{{-1}}$'.format(eta0,sfr0,vc0))

        qint, = ax_flux.plot(h, intensity0, lw=lw, color=color, linestyle="-")

        if data_saving:

            np.savetxt(os.path.join(path, f"final_dep_m{np.log10(mstar0):.1f}.dat"), (h, sigma0, intensity0))


    alpine = ax_flux.errorbar(data.x / (1000 * nc.pc), data.data, yerr=data.err, \
                          markerfacecolor='maroon', markeredgecolor='maroon', marker='o', \
                          linestyle='', ecolor='maroon')

    factor_data = data.data[0] / data.beam[0]

    ax_flux.plot(data.x_beam / 1e3 / nc.pc, data.beam * factor_data, linestyle="-.", color="gray",lw=2)

    #ax_sigma.legend()

    unc_down = np.full(h.shape, 1e-4)
    F0 = 0.110
    r0 = 3.
    unc_up = np.full(h.shape, F0)
    h_up_data = np.asarray([3.,6.,9.,12.])
    unc_up_data = np.asarray([0.110, 0.110*0.725, 0.110*0.43, 0.110*0.32])
    index_beam = np.searchsorted(h, r0)
    unc_up[index_beam::] = np.interp(h[index_beam::], h_up_data, unc_up_data)
    # unc_up[index_beam+1::] = F0 * (1./(np.pi*h[index_beam+1::]**2 - np.pi*h[index_beam]**2))**0.5
    # ax_flux.scatter(h_up_data, unc_up_data)
    ax_flux.fill_between(h, unc_down, unc_up,
                           alpha=0.15, color="gray")

    plt.subplots_adjust(left=0.17,  # the left side of the subplots of the figure
                        right=0.95,  # the right side of the subplots of the figure
                        bottom=0.11,  # the bottom of the subplots of the figure
                        top=0.96,  # the top of the subplots of the figure
                        wspace=0.05,  # the amount of width reserved for space between subplots,
                        # expressed as a fraction of the average axis width
                        hspace=0.1)  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height


    plt.show()

