import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.interpolate import interp1d


from mcmc_new import get_emission_fast, get_other_params
from mcmc_new import h, grid
from load_data_updated import obs_data_list, names, names_CII_halo, observational_data_fuji


import plot_config as pltc
import mydir
import my_utils
import natconst as nc


int_data = int(input("data number?"))
data = obs_data_list[int_data]



other_params = get_other_params(data.params_obs["redshift"], data.params_obs["line_FWHM"])

beam_interp = np.interp(h, data.x_beam / 1e3 / nc.pc, data.beam, right=0.)

beam_interp[beam_interp < 0.] = 0.

beam_func = interp1d(h, beam_interp, \
                     fill_value=(beam_interp[0], 0.), bounds_error=False)

beam_2d = beam_func(np.sqrt(grid[0] ** 2 + grid[1] ** 2))
f_beam = np.fft.fft2(beam_2d)





axd = plt.figure(figsize=(9, 5.5), constrained_layout=False).subplot_mosaic(
    """
    A
    """,
    gridspec_kw={
        "left": 0.12,
        "right": 0.97,
        "bottom": 0.25,
        "top": 0.96,
        "wspace": 0.27,
        "hspace": 0.4
    },
)

ax = axd["A"]

ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_ylim((1e-7, 1e0))
#ax.set_xlim((0.3, 16))
ax.tick_params(length=4, axis="both", which="major")
ax.tick_params(length=2, axis="both", which="minor")

ax.set_xlabel("b [kpc]")
ax.set_ylabel('flux [mJy/arcsec$^2$]')
#ax.xaxis.set_minor_locator(MultipleLocator(1))

# INIT PHASE

#eta0 = 5.
#sfr0 = 50.
#vc0 = 200.

mstar0 =1e10
sfr0 = 5*(mstar0/1e9)
eta0 = 10**(-0.36*np.log10(mstar0/1e10) + np.log10(4.9))

from scipy import optimize
def zero_me(x):
  return mstar0 - my_utils.mstar_behroozi(x, z=5.)

mvir0 = optimize.newton(zero_me,1.e+11)

r_axis = np.linspace(0.3, 20,1000)
print(np.log10(mvir0))
vc0 = my_utils.get_vc_from_virial_mass(mvir0, z=5.)/1e5
params0 = np.asarray([np.log10(eta0), np.log10(sfr0), np.log10(vc0)])
print(eta0, sfr0, vc0)
intensity, sigma, n0, r0 = get_emission_fast(params0, data, other_params, h, grid, f_beam, return_sigma=True)

#q, = ax.plot(h, sigma,lw=2, color='C0', linestyle="-")

#q0 =  ax.plot(h, sigma,lw=2, color='gray', linestyle="--")

q, = ax.plot(r_axis, np.interp(r_axis, r0,n0, left=0., right=0.),lw=2, color='C0', linestyle="-")

q0 =  ax.plot(r_axis, np.interp(r_axis, r0,n0, left=0., right=0.),lw=2, color='gray', linestyle="--")

#alpine = ax.errorbar(data.x / (1000 * nc.pc), data.data, yerr=data.err, \
 #                              markerfacecolor='maroon', markeredgecolor='maroon', marker='o', \
 #                              linestyle='', ecolor='maroon')

factor_data = data.data[0] / data.beam[0]

ax.plot(data.x_beam / 1e3 / nc.pc, data.beam * factor_data, linestyle="-.", color="gray")

# SLIDERS

axcolor = 'lightgoldenrodyellow'

#axeta = plt.axes([0.22, 0.11, 0.65, 0.03], facecolor=axcolor)
#axsfr = plt.axes([0.22, 0.06, 0.65, 0.03], facecolor=axcolor)
#axvc = plt.axes([0.22, 0.01, 0.65, 0.03], facecolor=axcolor)
axmstar = plt.axes([0.22, 0.01, 0.65, 0.03], facecolor=axcolor)
#

bounds = [(0.1, 10.), (10., 150.), (100., 300.)]

#s_eta = Slider(axeta, r"eta", bounds[0][0], bounds[0][1], valinit=eta0, valstep=0.1)
#s_sfr = Slider(axsfr, 'SFR', bounds[1][0], bounds[1][1], valinit=sfr0, valstep=5.)
#s_vc = Slider(axvc, 'v_c', bounds[2][0], bounds[2][1], valinit=vc0, valstep=5.)
s_mstar = Slider(axmstar, 'log mstar', 8.,12., valinit=np.log10(mstar0), valstep=0.01)


def update(val):
    #eta = s_eta.val
    #sfr = s_sfr.val
    #vc = s_vc.val
    mstar = 10**(s_mstar.val)

    sfr = 5 * (mstar / 1e9)
    eta = 10**(-0.36 * np.log10(mstar / 1e10) + np.log10(4.9))

    def zero_me1(x):
        return mstar - my_utils.mstar_behroozi(x, z=5.)

    mvir = optimize.newton(zero_me1, 10*mstar)

    vc = my_utils.get_vc_from_virial_mass(mvir, z=5.)/1e5

    params = np.asarray([np.log10(eta), np.log10(sfr), np.log10(vc)])

    intensity, sigma, n, r = get_emission_fast(params, data, other_params, h, grid, f_beam, return_sigma=True)

    q.set_ydata(np.interp(r_axis, r, n, left=0., right=0.))

    print(eta, sfr, vc)

#s_eta.on_changed(update)
#s_sfr.on_changed(update)
#s_vc.on_changed(update)
s_mstar.on_changed(update)

plt.show()