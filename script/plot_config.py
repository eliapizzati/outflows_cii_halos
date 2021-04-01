# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 15:39:03 2021

@author: anna
"""

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

import natconst as nc


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
        "xtick.minor.size": 12.0,
        "ytick.minor.size": 12.0,
 #       "legend.fontsize": 12.0,
 #       "figure.dpi": 200,
 #       "savefig.dpi": 200,
#        "text.usetex": True
})


size = 14


def plot_configurator(plot_type = "sol", bottom=0.28, wspace=0.1):


    if plot_type == "sol":
        fig_sol, axs_sol = plt.subplots(1, 3, sharex=True, figsize=(1.5*12.27,1.5*5.))
        
        ax_v = axs_sol[0]
        ax_n = axs_sol[1]
        ax_T = axs_sol[2]
        
        ax_v.set_xlabel("log (r [kpc])", size=size)
        ax_v.set_ylabel("v [1000 km/s] ", size=size)
        ax_v.set_xlim((np.log10(0.3),np.log10(30)))
        
        ax_n.set_xlabel("log (r [kpc])", size=size)
        ax_n.set_ylabel("log (n [cm$^{-3}$]) ", size=size)
        
        ax_T.set_xlabel("log (r [kpc])", size=size)
        ax_T.set_ylabel("log (T [K])", size=size)
        ax_T.set_ylim((2,8))
    
        plt.subplots_adjust(left = 0.1,  # the left side of the subplots of the figure
            right = 0.95,   # the right side of the subplots of the figure
            bottom = bottom,  # the bottom of the subplots of the figure
            top = 0.9,     # the top of the subplots of the figure
            wspace = wspace,  # the amount of width reserved for space between subplots,
            # expressed as a fraction of the average axis width
            hspace = 0.1)  # the amount of height reserved for space between subplots,
                          # expressed as a fraction of the average axis height

        return fig_sol, axs_sol    


    elif plot_type == "ion":

        fig_ion, axs_ion = plt.subplots(1, 2, sharex=True, figsize=(1.5*8.27,1.5*4.))

        ax_xe = axs_ion[0]
        ax_xCII = axs_ion[1]
        
        ax_xe.set_xlabel("log (r [kpc])", size=size)
        ax_xe.set_ylabel("n$_\\mathrm{HI}$/n$_\\mathrm{H}$", size=size)
        ax_xe.set_ylim((0,1))
        ax_xe.set_xlim((np.log10(0.3),np.log10(30)))
        
        ax_xCII.set_xlabel("log (r [kpc])", size=size)
        ax_xCII.set_ylabel("n$_\\mathrm{CII}$/n$_\\mathrm{C}$", size=size)
        ax_xCII.ticklabel_format(axis='y', style='plain')
        ax_xCII.set_xlim((np.log10(0.3),np.log10(30)))
        ax_xCII.set_ylim((0,1))

        plt.subplots_adjust(left = 0.1,  # the left side of the subplots of the figure
            right = 0.95,   # the right side of the subplots of the figure
            bottom = bottom,  # the bottom of the subplots of the figure
            top = 0.9,     # the top of the subplots of the figure
            wspace = wspace,  # the amount of width reserved for space between subplots,
            # expressed as a fraction of the average axis width
            hspace = 0.1)  # the amount of height reserved for space between subplots,
                          # expressed as a fraction of the average axis height

        return fig_ion, axs_ion
    
    elif plot_type == "sigma": 
        
        fig_sigma, ax_sigma = plt.subplots(1, 1, sharex=True, figsize=(1.5*8.27,1.5*4.))
            
        ax_sigma.set_xlabel("b [kpc]", size=size)
        ax_sigma.set_ylabel(r"log ($\Sigma_{CII}$ [erg/cm s^2])", size=size)
        ax_sigma.set_yscale("log")
        ax_sigma.set_xlim((0.,16))

        plt.subplots_adjust(left = 0.1,  # the left side of the subplots of the figure
            right = 0.95,   # the right side of the subplots of the figure
            bottom = bottom,  # the bottom of the subplots of the figure
            top = 0.9,     # the top of the subplots of the figure
            wspace = wspace,  # the amount of width reserved for space between subplots,
            # expressed as a fraction of the average axis width
            hspace = 0.1)  # the amount of height reserved for space between subplots,
                          # expressed as a fraction of the average axis height

        return fig_sigma, ax_sigma
    
    elif plot_type == "int":
    
        fig_int, ax_int = plt.subplots(1, 1, sharex=True, figsize=(1.5*8.27,1.5*4.))
            
        ax_int.set_xlabel("b [kpc]", size=size)
        ax_int.set_ylabel("flux [mJy/arcsec$^2$]", size=size)
        ax_int.set_yscale("log")
        ax_int.set_xlim((0.,16))

        plt.subplots_adjust(left = 0.1,  # the left side of the subplots of the figure
            right = 0.95,   # the right side of the subplots of the figure
            bottom = bottom,  # the bottom of the subplots of the figure
            top = 0.9,     # the top of the subplots of the figure
            wspace = wspace,  # the amount of width reserved for space between subplots,
            # expressed as a fraction of the average axis width
            hspace = 0.1)  # the amount of height reserved for space between subplots,
                          # expressed as a fraction of the average axis height
        
        return fig_int, ax_int
    
    else:
        raise ValueError("No correct type")
        

