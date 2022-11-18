import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.interpolate import interp1d


from mcmc_fast_emission import get_emission_fast, get_other_params
from mcmc_main import h, grid
from load_data import obs_data_list, names, names_CII_halo, observational_data_fuji


import plot_config as pltc
import my_dir
import my_utils
import natconst as nc

import matplotlib

matplotlib.rcParams.update({
        "font.size": 12.0,
        "font.family": 'sans-serif',
#        "font.sans-serif": ['Helvetica'],
        "axes.titlesize": 12.,
        "axes.labelsize": 12.,
        "xtick.labelsize": 12.,
        "ytick.labelsize": 12.,
        "xtick.major.size": 4.0,
        "ytick.major.size": 4.0,
        "xtick.minor.size": 2.0,
        "ytick.minor.size": 2.0,
        "legend.fontsize": 12.0,
        "legend.frameon": False
 #       "figure.dpi": 200,
 #       "savefig.dpi": 200,
#        "text.usetex": True
})


plot1 = True #single dependence
plot2 = False #multiple dependencies
plot3 = False #multiple dependencies w alpha


int_data = 0#int(input("data number?"))
data = obs_data_list[int_data]



other_params = get_other_params(data.params_obs["redshift"], data.params_obs["line_FWHM"])

beam_interp = np.interp(h, data.x_beam / 1e3 / nc.pc, data.beam, right=0.)

beam_interp[beam_interp < 0.] = 0.

beam_func = interp1d(h, beam_interp, \
                     fill_value=(beam_interp[0], 0.), bounds_error=False)

beam_2d = beam_func(np.sqrt(grid[0] ** 2 + grid[1] ** 2))
f_beam = np.fft.fft2(beam_2d)

axd = plt.figure(figsize=(15, 8.5), constrained_layout=False).subplot_mosaic(
    """
    AC
    AC
    AD
    BD
    BE
    BE
    """,
    gridspec_kw={
        "left": 0.1,
        "right": 0.97,
        "bottom": 0.13,
        "top": 0.96,
        "wspace": 0.2,
        "hspace": 0.35
    },
)

ax_flux = axd["A"]
ax_sigma = axd["B"]
ax_n = axd["C"]
ax_v = axd["D"]
ax_T = axd["E"]

# ax_sigma.set_xscale('log')
ax_sigma.set_yscale('log')
ax_sigma.set_ylim((1e-7, 1e0))
ax_sigma.set_xlim((0.3, 16))
ax_sigma.tick_params(length=4, axis="both", which="major")
ax_sigma.tick_params(length=2, axis="both", which="minor")

ax_sigma.set_xlabel("b [kpc]")
ax_sigma.set_ylabel('Sigma [erg cm-2 s-1]')
# ax.xaxis.set_minor_locator(MultipleLocator(1))

# ax_flux.set_xscale('log')
ax_flux.set_yscale('log')
ax_flux.set_ylim((1e-3, 1e1))
ax_flux.set_xlim((0.3, 16))
ax_flux.tick_params(length=4, axis="both", which="major")
ax_flux.tick_params(length=2, axis="both", which="minor")

ax_flux.set_xlabel("b [kpc]")
ax_flux.set_ylabel('flux [mJy/arcsec$^2$]')
# ax.xaxis.set_minor_locator(MultipleLocator(1))


ax_v.set_ylabel("v [1000 km/s] ")
ax_v.set_ylim((0, 0.8))

ax_v.set_xlim((0.3, 30))
ax_T.set_xlim((0.3, 30))
ax_n.set_xlim((0.3, 30))

ax_n.set_ylabel("n [cm$^{-3}$]) ")
ax_n.set_ylim((1e-5, 1e4))

ax_T.set_xlabel("r [kpc])")
ax_T.set_ylabel("T [K])")
ax_T.set_ylim((1e1, 1e7))

ax_n.set_xscale('log')
ax_n.set_yscale('log')
ax_v.set_xscale('log')
# ax_v.set_yscale('log')
ax_T.set_xscale('log')
ax_T.set_yscale('log')

# INIT PHASE
alfa0 = 1.

mstar0 = 1e10

eta0 = my_utils.from_xfitted_to_eta(mstar0 / 1e10, 4.9, -0.36)
sfr0 = my_utils.from_eta_to_xfitted(eta0, 14.2, -0.26)

mvir0 = my_utils.mvir_behroozi(mstar0, z=5.)

print(np.log10(mvir0))

vc0 = my_utils.get_vc_from_virial_mass(mvir0, z=5.) / 1e5

params0 = np.asarray([np.log10(eta0), np.log10(sfr0), np.log10(vc0)])

print(eta0, sfr0, vc0)

intensity0, sigma0, r0, n0, v0, T0 = get_emission_fast(params0, data, other_params, h, grid, f_beam,
                                                       return_quantities="all", alfa=alfa0)

qsigma, = ax_sigma.plot(h, sigma0, lw=2, color='C0', linestyle="-")
qsigma0 = ax_sigma.plot(h, sigma0, lw=2, color='gray', linestyle="--")

qint, = ax_flux.plot(h, intensity0, lw=2, color='C0', linestyle="-")
qint0 = ax_flux.plot(h, intensity0, lw=2, color='gray', linestyle="--")

r_axis = np.linspace(0.3, 20, 1000)

qn, = ax_n.plot(r_axis, np.interp(r_axis, r0, n0, left=0., right=0.), lw=2, color='C0', linestyle="-")
qn0 = ax_n.plot(r_axis, np.interp(r_axis, r0, n0, left=0., right=0.), lw=2, color='gray', linestyle="--")

qv, = ax_v.plot(r_axis, np.interp(r_axis, r0, v0 / 1e8, left=0., right=0.), lw=2, color='C0', linestyle="-")
qv0 = ax_v.plot(r_axis, np.interp(r_axis, r0, v0 / 1e8, left=0., right=0.), lw=2, color='gray', linestyle="--")

qT, = ax_T.plot(r_axis, np.interp(r_axis, r0, T0, left=0., right=0.), lw=2, color='C0', linestyle="-")
qT0 = ax_T.plot(r_axis, np.interp(r_axis, r0, T0, left=0., right=0.), lw=2, color='gray', linestyle="--")

alpine = ax_flux.errorbar(data.x / (1000 * nc.pc), data.data, yerr=data.err, \
                          markerfacecolor='maroon', markeredgecolor='maroon', marker='o', \
                          linestyle='', ecolor='maroon')

factor_data = data.data[0] / data.beam[0]

ax_flux.plot(data.x_beam / 1e3 / nc.pc, data.beam * factor_data, linestyle="-.", color="gray")

text = plt.figtext(0.63, 0.023,
                   s=f'eta ={eta0:.2f}   SFR = {sfr0:.0f}Msun/yr    vc = {vc0:.0f}kms', \
                   fontsize=14, bbox={'facecolor': 'lightgray', 'alpha': 0.8, \
                                      'boxstyle': "round", 'edgecolor': 'none', 'pad': 0.5})

if plot1:
    "single dep"

    # SLIDERS

    axcolor = 'lightgoldenrodyellow'

    axmstar = plt.axes([0.12, 0.02, 0.45, 0.02], facecolor=axcolor)

    s_mstar = Slider(axmstar, 'log mstar', 9.,11., valinit=np.log10(mstar0), valstep=0.01)


    def update(val):
        mstar = 10**(s_mstar.val)

        eta = my_utils.from_xfitted_to_eta(mstar/1e10, 4.9, -0.36)
        sfr = my_utils.from_eta_to_xfitted(eta, 14.2, -0.26)
        mvir = my_utils.mvir_behroozi(mstar, z=5.)
        vc = my_utils.get_vc_from_virial_mass(mvir, z=5.)/1e5

        params = np.asarray([np.log10(eta), np.log10(sfr), np.log10(vc)])

        intensity, sigma, r, n, v, T = get_emission_fast(params, data, other_params, h, grid, f_beam, return_quantities="all")

        qsigma.set_ydata(sigma)
        qint.set_ydata(intensity)
        qn.set_ydata(np.interp(r_axis, r, n, left=0., right=0.))
        qv.set_ydata(np.interp(r_axis, r, v/1e8, left=0., right=0.))
        qT.set_ydata(np.interp(r_axis, r, T, left=0., right=0.))

        text.set_text(f'eta ={eta:.2f}   SFR = {sfr:.0f}Msun/yr   vc = {vc:.0f}kms')

        print(eta, sfr, vc)

    s_mstar.on_changed(update)

    plt.show()



if plot2:
    "multiple deps"

    # SLIDERS

    axcolor = 'lightgoldenrodyellow'

    #
    #
    axeta = plt.axes([0.12, 0.05, 0.45, 0.013], facecolor=axcolor)
    axsfr = plt.axes([0.12, 0.03, 0.45, 0.013], facecolor=axcolor)
    #axvc = plt.axes([0.12, 0.01, 0.45, 0.01], facecolor=axcolor)
    axmstar = plt.axes([0.12, 0.01, 0.45, 0.013], facecolor=axcolor)


    bounds = [(0.1, 10.), (10., 150.), (100., 300.)]

    s_eta = Slider(axeta, r"eta", bounds[0][0], bounds[0][1], valinit=eta0, valstep=0.1)
    s_sfr = Slider(axsfr, 'SFR', bounds[1][0], bounds[1][1], valinit=sfr0, valstep=5.)
    # s_vc = Slider(axvc, 'v_c', bounds[2][0], bounds[2][1], valinit=vc0, valstep=5.)
    s_mstar = Slider(axmstar, 'log mstar', 8., 12., valinit=np.log10(mstar0), valstep=0.01)


    def update(val):
        mstar = 10 ** (s_mstar.val)
        eta = s_eta.val
        sfr = s_sfr.val
        mvir = my_utils.mvir_behroozi(mstar, z=5.)
        vc = my_utils.get_vc_from_virial_mass(mvir, z=5.) / 1e5

        params = np.asarray([np.log10(eta), np.log10(sfr), np.log10(vc)])

        intensity, sigma, r, n, v, T = get_emission_fast(params, data, other_params, h, grid, f_beam,
                                                         return_quantities="all")

        qsigma.set_ydata(sigma)
        qint.set_ydata(intensity)
        qn.set_ydata(np.interp(r_axis, r, n, left=0., right=0.))
        qv.set_ydata(np.interp(r_axis, r, v / 1e8, left=0., right=0.))
        qT.set_ydata(np.interp(r_axis, r, T, left=0., right=0.))

        text.set_text(f'eta ={eta:.2f}   SFR = {sfr:.0f}Msun/yr   vc = {vc:.0f}kms')

        print(eta, sfr, vc)


    s_eta.on_changed(update)
    s_sfr.on_changed(update)
    #s_vc.on_changed(update)
    s_mstar.on_changed(update)

    plt.show()


if plot3:
    "multiple deps w alpha"

    # SLIDERS

    axcolor = 'lightgoldenrodyellow'

    #
    #
    axalpha = plt.axes([0.12, 0.0625, 0.45, 0.01], facecolor=axcolor)
    axeta = plt.axes([0.12, 0.045, 0.45, 0.01], facecolor=axcolor)
    axsfr = plt.axes([0.12, 0.0275, 0.45, 0.01], facecolor=axcolor)
    axmstar = plt.axes([0.12, 0.01, 0.45, 0.01], facecolor=axcolor)


    bounds = [(0.1, 10.), (10., 150.), (100., 300.)]

    s_alpha = Slider(axalpha, 'alpha', 0.05, 2., valinit=alfa0, valstep=0.05)
    s_eta = Slider(axeta, r"eta", bounds[0][0], bounds[0][1], valinit=eta0, valstep=0.1)
    s_sfr = Slider(axsfr, 'SFR', bounds[1][0], bounds[1][1], valinit=sfr0, valstep=5.)
    s_mstar = Slider(axmstar, 'log mstar', 8., 12., valinit=np.log10(mstar0), valstep=0.01)


    def update(val):
        mstar = 10 ** (s_mstar.val)
        alfa = s_alpha.val
        eta = s_eta.val
        sfr = s_sfr.val
        mvir = my_utils.mvir_behroozi(mstar, z=5.)
        vc = my_utils.get_vc_from_virial_mass(mvir, z=5.) / 1e5

        params = np.asarray([np.log10(eta), np.log10(sfr), np.log10(vc)])

        intensity, sigma, r, n, v, T = get_emission_fast(params, data, other_params, h, grid, f_beam,
                                                         return_quantities="all", alfa=alfa)

        qsigma.set_ydata(sigma)
        qint.set_ydata(intensity)
        qn.set_ydata(np.interp(r_axis, r, n, left=0., right=0.))
        qv.set_ydata(np.interp(r_axis, r, v / 1e8, left=0., right=0.))
        qT.set_ydata(np.interp(r_axis, r, T, left=0., right=0.))

        text.set_text(f'eta ={eta:.2f}   SFR = {sfr:.0f}Msun/yr   vc = {vc:.0f}kms')

        print(eta, sfr, vc)

    s_alpha.on_changed(update)
    s_eta.on_changed(update)
    s_sfr.on_changed(update)
    s_mstar.on_changed(update)

    plt.show()