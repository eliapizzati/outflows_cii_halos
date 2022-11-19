"""
This script was used to plot some figures for the thesis (e.g. power spectrum, starburst99 emission, nonlinear collapse, ecc
"""



import numpy as np
import matplotlib.pyplot as plt
import natconst as nc

from colossus.cosmology import cosmology
from scipy.optimize import curve_fit
from colossus.lss import mass_function

import scipy.optimize as so

from functools import partial

import matplotlib

matplotlib.rcParams.update({
        "font.size": 16.0,
        "font.family": 'sans-serif',
#        "font.sans-serif": ['Helvetica'],
        "axes.titlesize": 16.,
        "axes.labelsize": 16.,
        "xtick.labelsize": 16.,
        "ytick.labelsize": 16.,
        "xtick.major.size": 6.0,
        "ytick.major.size": 6.0,
        "xtick.minor.size": 4.0,
        "ytick.minor.size": 4.0,
        "legend.fontsize": 13.0,
        "legend.frameon": False
 #       "figure.dpi": 200,
 #       "savefig.dpi": 200,
#        "text.usetex": True
})


cosmo = cosmology.setCosmology('planck15')

size = 14

plot1 = False   # sigma_m & ps
plot2 = False   # Mass function & density profile
plot3 = True    # starburst99


if plot1:
    """
    ## sigma_m & ps
    """
    
    redshifts = [0.,3.,10.,25.,50.,100.]
    
    
    fig, [ax_ps, ax_sigma] = plt.subplots(1,2,figsize=(1.3*8.27,1.3*3.2))
    
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    
    ax_ps.set_xlabel('log( k / h Mpc$^{-1}$)',labelpad=-0.8)
    ax_ps.set_ylabel('log( P(k) / h$^3$ Mpc$^3$)')
    
    
    k = 10**np.arange(-5,2,0.02)
    
    
    for z in redshifts:
        
        Pk = cosmo.matterPowerSpectrum(k, z=z)
    
        ax_ps.plot(np.log10(k), np.log10(Pk), label = 'z={:.0f}'.format(z))
    
    
    ax_ps.legend(frameon=False, ncol=2)# fontsize='large')
    
    #ax_sigma.set_xscale('log')
    #ax_sigma.set_yscale('log')
    
    ax_sigma.set_xlabel('log( M / M$_\odot$ )')
    ax_sigma.set_ylabel('log( $\sigma$ (M) )')
    
    R = 10**np.arange(-2.5,2,0.005)
    
    M = (R/1.85)**3 * 1e12
    
    for z in redshifts:
        
        sigma_tophat = cosmo.sigma(R=R, z=z)
    
        ax_sigma.plot(np.log10(M), np.log10(sigma_tophat), label = 'z={:.0f}'.format(z))
    
          
    plt.subplots_adjust(left = 0.08,  # the left side of the subplots of the figure
    right = 0.98,   # the right side of the subplots of the figure
    bottom = 0.16,  # the bottom of the subplots of the figure
    top = 0.95,     # the top of the subplots of the figure
    wspace = 0.2,  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    hspace = 0.3)  # the amount of height reserved for space between subplots,
                  # expressed as a fraction of the average axis height




if plot2:
    """
    #Mass function & density profile
    """

    fig, [ax_mass, ax_profile] = plt.subplots(1,2,figsize=(1.3*8.27,1.3*3.2))

    redshifts = [0.,3.,5.,10.,20.,30.]
    
    M = 10**np.arange(2.0, 16.0, 0.01)
    
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    
    ax_mass.set_xlabel('log( M / M$_\odot$)')
    ax_mass.set_ylabel('log( dn / dlog(M) [Mpc$^3$])', size=size)
    
    ax_mass.set_ylim(-9, 6)
    
    for z in redshifts:
        mfunc = mass_function.massFunction(M, z, q_out='dndlnM', mdef = 'fof', model = 'press74')
        
        ax_mass.plot(np.log10(M), np.log10(mfunc), '-', label = 'z = {:.0f}'.format(z))
        
    ax_mass.legend()
    
    
    radius = np.logspace(-2,2)
    
    density_nfw = 4/radius/(1+radius)**2
    density_iso = radius**(-2)
    
    ax_profile.plot(np.log10(radius), np.log10(density_nfw), label = "NFW")
    
    ax_profile.plot(np.log10(radius), np.log10(density_iso), label = "iso sphere")
    
    ax_profile.set_xlabel('log( r / r$_s$ )')
    ax_profile.set_ylabel(r'log( $\rho$ / $\rho(r_s)$ )')
  
    ax_profile.legend()

    
    
    plt.subplots_adjust(left = 0.1,  # the left side of the subplots of the figure
    right = 0.98,   # the right side of the subplots of the figure
    bottom = 0.16,  # the bottom of the subplots of the figure
    top = 0.95,     # the top of the subplots of the figure
    wspace = 0.2,  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    hspace = 0.3)  # the amount of height reserved for space between subplots,
                  # expressed as a fraction of the average axis height




if plot3:
    """
    # starburst99
    """
    
    fig, ax = plt.subplots(1,1,figsize=(8,5))

    
    data = np.loadtxt("starburst_data.dat", unpack=True)
    
    indices = [1,3,7,12,23,29,33,36] 
    times = [1,3,5,10,30,100,500,900]

    for indice, time in zip(indices,times): 
        ax.plot(np.log10(data[0]), data[indice], label="{} Myr".format(time))
    
        
    ax.set_xlabel(r'log( Wavelength [${\rm\AA}$] )')
    ax.set_ylabel(r'log( Luminosity [erg s$^{-1}$ ${\rm \AA}^{-1}$] )')
    
    ax.set_xlim(2,4)
    ax.set_ylim(37,41)
    
    ax.legend(ncol=3, fontsize=15)

    
       
    plt.subplots_adjust(left = 0.15,  # the left side of the subplots of the figure
    right = 0.95,   # the right side of the subplots of the figure
    bottom = 0.15,  # the bottom of the subplots of the figure
    top = 0.95,     # the top of the subplots of the figure
    wspace = 0.2,  # the amount of width reserved for space between subplots,
    # expressed as a fraction of the average axis width
    hspace = 0.3)  # the amount of height reserved for space between subplots,
                  # expressed as a fraction of the average axis height

 

"""
# non-linear collapse
"""
#t_i = 3e5 * nc.year
#
#delta_i = 1e-3
#
#Mass = 1e6 * nc.ms # grammi
#
#
#radius_max = (243.* (nc.gg*Mass*t_i**2)/250)**(1./3) / delta_i
#
#
#A = radius_max /2
#
#
#B = np.sqrt(A**3/(nc.gg*Mass))
#
#
#A = 1
#
#B = 1
#
#theta = np.linspace(0.,2*np.pi, 1000)
#
#theta = theta[:-1:]
#
#radius = A*(1-np.cos(theta))
#
#time = B*(theta-np.sin(theta)) # + t_i
#
#
#delta = (9./2) * ((theta-np.sin(theta))**2 / (1-np.cos(theta))**3) - 1.
#
#
#fig, ax = plt.subplots()
#
#
#
##ax.set_xscale('log')
##ax.set_yscale('log')
#
##ax.set_xlabel('time',size=size)
##ax.set_ylabel('radius', size=size)
#
#ax.tick_params(labelsize=size)
#
#k = 10**np.arange(-5,2,0.02)
#
#ax.set_xticks((0.,np.pi, 2*np.pi))
#ax.set_xticklabels(("0","$\pi$", "$2\pi$"))
#
#r_max = np.max(radius)
#r_vir = r_max/2.
#
#ax.set_yticks((0.,r_vir, r_max))
#ax.set_yticklabels(("r$_i$", "r$_{vir}$", "r$_m$"))
#
#
#
#ax.plot(time, [0]*len(time), linestyle='--', color='gray')
#ax.plot(time, [r_vir]*len(time), linestyle='--', color='gray')
#ax.plot(time, [r_max]*len(time), linestyle='--', color='gray')
#
#
#ax.plot(time, radius)
#
#
#
#fig, ax = plt.subplots()
#
#ax.plot(time, np.log10(delta), label='real growth')
##ax.set_yscale("log")
#ax.set_ylim(-2,2.9)
#
#
#def func(x, a):
#     return a * x**(2./3)
# 
#time = time[1::]
#delta = delta[1::]
#
#xdata = time[time<np.pi/16]
#ydata = delta[time<np.pi/16]
#
#popt, pcov = curve_fit(func, xdata, ydata)
#a = popt
#
#ax.plot(time, np.log10(func(time,a)), label='linear growth')
#
#ax.set_xticks((0.,np.pi, 2*np.pi))
#ax.set_xticklabels(("0","$\pi$", "$2\pi$"))
#
#ax.plot([0]*len(delta), np.log10(delta), linestyle='--', color='gray')
#ax.plot([np.pi]*len(delta), np.log10(delta), linestyle='--', color='gray')
#ax.plot([np.pi*2]*len(delta), np.log10(delta), linestyle='--', color='gray')
#
##ax.set_xscale('log')
##ax.set_yscale('log')
#
##ax.set_xlabel('time',size=size)
#ax.set_ylabel('log( $\delta$ )', size=size)
#
#ax.tick_params(labelsize=size)
#
#ax.legend(frameon=False, fontsize='large')
#





"""
# collapsing mass
"""

#
#delta_crit = 1.69
#
#
#redshifts = np.linspace(1.,30., 1000)
#
#
#fig, ax = plt.subplots()
#
##ax.set_xscale('log')
##ax.set_yscale('log')
#
#ax.set_xlabel('log( M / M$_\odot$)',size=size)
#ax.set_ylabel('$\sigma$ (M)', size=size)
#
#ax.tick_params(labelsize=size)
#
#
#def sol_tophat(R, z): 
#    cosmo.sigma(R=R, z=z) - delta_crit
#
#def sol(z):
#    return so.brentq(partial(sol_tophat,z),0.01,10)
#    
#def sol_vectorized(redshifts):
#    
#    vector = []
#    
#    for z in redshifts:
#        vector.append(sol(z))
#
#    return vector
#    
#ax.plot(redshifts, sol_vectorized(redshifts))
#
#
##ax.legend(frameon=False, fontsize='large')



"""
# sigma_M and delta_crit
"""



#
#redshifts = [0.,2.,5.,10.,20.,30.]
#
#
#fig, ax1 = plt.subplots()
#
##ax.set_xscale('log')
##ax.set_yscale('log')
#
#ax1.set_xlabel('log( M / M$_\odot$)',size=size)
#ax1.set_ylabel('log( $\sigma$ (M) )', size=size)
#
#ax1.tick_params(labelsize=size)
#
#
#R = 10**np.arange(-2.5,1.,0.005)
#
#M = (R/1.85)**3 * 1e12
#
#sigma_tophat_0 = cosmo.sigma(R=R, z=0.)
#
#ax1.plot(np.log10(M), np.log10(sigma_tophat_0), color='gray')
#
#
#
#for z in redshifts:
#    
#    sigma_tophat_z = cosmo.sigma(R=R, z=z)
#    
#    delta_crit = 1.69*sigma_tophat_0/sigma_tophat_z
#
#    ax1.plot(np.log10(M), np.log10(delta_crit), linestyle='--',alpha=0.8, label='z = {:.0f}'.format(z))
#
#
#ax1.legend(frameon=False, fontsize='large', ncol=2)
#
#
#fig.tight_layout()  # otherwise the right y-label is slightly clipped
#plt.show()
#
#

