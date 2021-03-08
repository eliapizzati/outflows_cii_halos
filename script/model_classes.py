# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 15:20:36 2021

@author: anna
"""



import os
import numpy as np
import mydir

import matplotlib.pyplot as plt

import natconst as nc

# CLASSES

class sol_profiles():
    
    def __init__(self, radius, variables, params):
        
        self.r = radius
        
        self.var = variables
        
        self.v = variables[0]
        self.n = variables[1]
        self.T = variables[2]
        
        self.params = params

    
    def to_file(self):
        
        folder = 'data_profiles'
        
        if not os.path.exists(os.path.join(mydir.data_dir, folder)):
            os.mkdir(os.path.join(mydir.data_dir, folder))

        np.savetxt(os.path.join(mydir.data_dir, folder, "profiles_beta{:.2f}_SFR{}_vc{:.1f}.dat".\
                                format(self.params["beta"], self.params["SFR"], self.params["v_c"])), \
                                (self.r,self.v,self.n,self.T))
        
        return
    
        
    def plot(self, size=14):
        
        folder = 'data_profiles'
        
        if not os.path.exists(os.path.join(mydir.plot_dir, folder)):
            os.mkdir(os.path.join(mydir.plot_dir, folder))


        fig, axs = plt.subplots(1, 3, sharex=True, figsize=(1.5*12.27,1.5*5.))

        ax_v = axs[0]
        ax_n = axs[1]
        ax_T = axs[2]

        ax_v.set_xlabel("log (r [kpc])", size=size)
        ax_v.set_ylabel("v [1000 km/s] ", size=size)
        ax_v.tick_params(labelsize=size)
        ax_v.set_xlim((np.log10(0.3),np.log10(30)))

        ax_n.set_xlabel("log (r [kpc])", size=size)
        ax_n.set_ylabel("log (n [cm$^{-3}$]) ", size=size)
        ax_n.tick_params(labelsize=size)

        ax_T.set_xlabel("log (r [kpc])", size=size)
        ax_T.set_ylabel("log (T [K])", size=size)
        ax_T.tick_params(labelsize=size)
        ax_T.set_ylim((2,8))
    
        ax_v.plot(np.log10(self.r/(1000*nc.pc)),self.v/10**8)       
        ax_n.plot(np.log10(self.r/(1000*nc.pc)),np.log10(self.n))
        ax_T.plot(np.log10(self.r/(1000*nc.pc)),np.log10(self.T))
        
        plt.savefig(os.path.join(mydir.data_dir, folder, "profiles_beta{:.2f}_SFR{}_vc{:.1f}.dat".\
                                format(self.params["beta"], self.params["SFR"], self.params["v_c"])))

        return
    
    
    def check_nans(self):
        
        var = np.asarray(self.var)
        
        return np.isnan(np.sum(var))
    
   
class ion_profiles():
    
    def __init__(self, radius, variables, params):
        
        self.r = radius
        
        self.var = variables
        
        self.x_H = variables[0]
        self.x_CII = variables[1]
        
        self.params = params


    def to_file(self):
        
        folder = 'data_ionization'
        
        if not os.path.exists(os.path.join(mydir.data_dir, folder)):
            os.mkdir(os.path.join(mydir.data_dir, folder))

        np.savetxt(os.path.join(mydir.plot_dir, folder, "ionization_beta{:.2f}_SFR{}_vc{:.1f}.dat".\
                                format(self.params["beta"], self.params["SFR"], self.params["v_c"])), \
                                (self.r,self.x_H,self.x_CII))
        
        return
    
        
    def plot(self, size=14):
        
        folder = 'data_ionization'
        
        if not os.path.exists(os.path.join(mydir.plot_dir, folder)):
            os.mkdir(os.path.join(mydir.plot_dir, folder))


        fig, axs = plt.subplots(1, 2, sharex=True, figsize=(1.5*8.27,1.5*4.))

        ax_xe = axs[0]
        ax_xCII = axs[1]
        
        ax_xe.set_xlabel("log (r [kpc])", size=size)
        ax_xe.set_ylabel("n$_\\mathrm{HI}$/n$_\\mathrm{H}$", size=size)
        ax_xe.tick_params(labelsize=size)
        ax_xe.set_ylim((0,1))
        ax_xe.set_xlim((np.log10(0.3),np.log10(30)))

        ax_xCII.set_xlabel("log (r [kpc])", size=size)
        ax_xCII.set_ylabel("n$_\\mathrm{CII}$/n$_\\mathrm{C}$", size=size)
        ax_xCII.tick_params(labelsize=size)
        ax_xCII.ticklabel_format(axis='y', style='plain')
        ax_xCII.set_ylim((0,1))

        ax_xe.plot(np.log10(self.r/(1000*nc.pc)),1.-self.x_e)
        ax_xCII.plot(np.log10(self.r/(1000*nc.pc)),self.x_CII)
        
        plt.savefig(os.path.join(mydir.plot_dir, folder, "ionization_beta{:.2f}_SFR{}_vc{:.1f}.dat".\
                                format(self.params["beta"], self.params["SFR"], self.params["v_c"])))

        return    
    
    
    def check_nans(self):
        
        var = np.asarray(self.var)
        
        return np.isnan(np.sum(var))
    

class lum_profile():
    
    def __init__(self, radius, variable, params):
        
        self.h = radius
        
        self.var = variable
        
        self.params = params

    
    def to_file(self):
        
        folder = 'data_profiles_SFR{}_vc{:.1f}'.format(self.params["SFR"], self.params["v_c"])
        
        if not os.path.exists(os.path.join(mydir.data_dir, folder)):
            os.mkdir(os.path.join(mydir.data_dir, folder))

        np.savetxt(os.path.join(mydir.data_dir, folder, "{}_beta{:.2f}.dat".format(self.params["type"], self.params["beta"])), \
                   (self.r,self.var))
        
        return
    
        
    def plot(self, ax):
        
        ax.plot(self.r, self.var)
        return
    
    
    def check_nans(self):
        
        return np.isnan(np.sum(self.var))
    

    