
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

from cooling import cooling_vectorized, heating_vectorized, temperatures
import mydir
import plot_config


fig, ax = plt.subplots(figsize=(12,15))

ax.set_xlabel("T [K]")
ax.set_ylabel(r"$|\Lambda$(T)| $[\mathrm{erg}\,\mathrm{cm}^3\,\mathrm{s}^{-1}]$")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim((1e-28, 1e-21))

n0 = 1.
plw0 = 0.
ph10 = 0.
pg10 = 0.
pc60 = 0.

logn_init = np.log10(n0)
logplw_init = -13
logph1_init = -13
logpg1_init = -13
logpc6_init = -20

T_vec, cooling_vec0 = cooling_vectorized(temperatures, n0, plw0, ph10, pg10, pc60)
T_vec, heating_vec0 = heating_vectorized(temperatures, n0, plw0, ph10, pg10, pc60)

T_vec, cooling_vec_init = cooling_vectorized(temperatures, 10**logn_init, 10**logplw_init, 10**logph1_init, \
                                             10**logpg1_init, 10**logpc6_init)
T_vec, heating_vec_init = heating_vectorized(temperatures, 10**logn_init, 10**logplw_init, 10**logph1_init, \
                                             10**logpg1_init, 10**logpc6_init)

c0 = plt.plot(T_vec, cooling_vec0, lw=1, color='gray')
h0 = plt.plot(T_vec, heating_vec0, lw=1, color='gray', linestyle='--')

c, = plt.plot(T_vec, cooling_vec_init, lw=1, color='C0')
h, = plt.plot(T_vec, heating_vec_init, lw=1, color='C0', linestyle='--')

plt.subplots_adjust(left=0.1,  # the left side of the subplots of the figure
                    right=0.95,  # the right side of the subplots of the figure
                    bottom=0.35,  # the bottom of the subplots of the figure
                    top=0.98,  # the top of the subplots of the figure
                    wspace=0.1,  # the amount of width reserved for space between subplots,
                    # expressed as a fraction of the average axis width
                    hspace=0.1)  # the amount of height reserved for space between subplots,
# expressed as a fraction of the average axis height


axcolor = 'lightgoldenrodyellow'
axn   = plt.axes([0.2, 0.01, 0.65, 0.03], facecolor=axcolor)
axplw = plt.axes([0.2, 0.06, 0.65, 0.03], facecolor=axcolor)
axph1 = plt.axes([0.2, 0.11, 0.65, 0.03], facecolor=axcolor)
axpg1 = plt.axes([0.2, 0.16, 0.65, 0.03], facecolor=axcolor)
axpc6 = plt.axes([0.2, 0.21, 0.65, 0.03], facecolor=axcolor)

sn   = Slider(axn, 'log density', -5, 3, valinit=logn_init, valstep=0.5)
splw = Slider(axplw, 'log plw', -30, -6, valinit=logplw_init, valstep=1.)
sph1 = Slider(axph1, 'log ph1', -30, -6, valinit=logph1_init, valstep=1.)
spg1 = Slider(axpg1, 'log pg1', -30, -6, valinit=logpg1_init, valstep=1.)
spc6 = Slider(axpc6, 'log pc6', -30, -6, valinit=logpc6_init, valstep=1.)


def update(val):
    logn = sn.val
    logplw = splw.val
    logph1 = sph1.val
    logpg1 = spg1.val
    logpc6 = spc6.val

    c.set_ydata(cooling_vectorized(temperatures, 10**logn, 10**logplw, 10**logph1, 10**logpg1, 10**logpc6)[1])
    h.set_ydata(heating_vectorized(temperatures, 10**logn, 10**logplw, 10**logph1, 10**logpg1, 10**logpc6)[1])
    fig.canvas.draw_idle()

sn.on_changed(update)
splw.on_changed(update)
sph1.on_changed(update)
spg1.on_changed(update)
spc6.on_changed(update)


plt.show()