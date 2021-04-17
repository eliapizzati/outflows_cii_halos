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

# general functions (to generate filenames&folders, and to load files)

def get_name_core(params):
    
    name_base = "_beta{:.2f}_SFR{:.1f}".format(params["beta"], params["SFR"])
    
    if params["NFW"] == True:
        
        pass
    
    elif params["NFW"] == False:

        name_specs = "_vc{:.1f}".format(params["v_c"])
            
        if params["f_esc_FUV"] != 0.0:
                
            name_specs += "fescfuv{:.2f}".format(params["f_esc_FUV"])
            
        if params["f_esc_ion"] != 0.0:
                
            name_specs += "fescion{:.2f}".format(params["f_esc_ion"])
    
    else:
            raise ValueError("No NFW switcher; you need to select a gravity model (isothermal sphere vs NFW)")
     
    return str(name_base + name_specs)

def get_folder(params, class_type):
    """
    class_type can be profiles, ionization, emission
    """
    
    if params["NFW"] == True:
            
            folder = 'data_{}_NFW'.format(class_type)
            
            if not os.path.exists(os.path.join(mydir.data_dir, folder)):
                os.mkdir(os.path.join(mydir.data_dir, folder))
            
    elif params["NFW"] == False:
            
            folder = 'data_{}'.format(class_type)

            if not os.path.exists(os.path.join(mydir.data_dir, folder)):
                os.mkdir(os.path.join(mydir.data_dir, folder))
            
    else:
            raise ValueError("No NFW switcher; you need to select a gravity model (isothermal sphere vs NFW)")

    return folder

def load_from_file(params, filename = None, class_type = "profiles", extension = ".dat", category = None):
    """
    class_type can be profiles, ionization, emission
    """
        
    if class_type == "profiles":
        
        name_prefix = "profiles"

        if filename is None:
            
            name_core = get_name_core(params)

            folder = get_folder(params, class_type)
            
            data = np.loadtxt(os.path.join(mydir.data_dir, folder, name_prefix + name_core + extension))     
                
        elif filename is not None:
            
            data = np.loadtxt(filename)
            
        try:
            r, v, n, T = data
        except:
            raise ValueError("No correct format for input data from sol_profiles")
                
        profile = sol_profiles(radius=r, variables=[v,n,T], params=params)

        
    elif class_type == "ionization":

        name_prefix = "ionization"
        
        if filename is None:
            
            name_core = get_name_core(params)

            folder = get_folder(params, class_type)
            
            data = np.loadtxt(os.path.join(mydir.data_dir, folder, name_prefix + name_core + extension))     
                
        elif filename is not None:
            
            data = np.loadtxt(filename)            

        try:
            r, x_e, x_CII = data
        except:
            raise ValueError("No correct format for input data from sol_profiles")
                
        profile = ion_profiles(radius=r, variables=[x_e, x_CII], params=params)

        
    elif class_type == "lum_profile":

        if category == "sigma":
            name_prefix = category
            
        elif category == "int_raw":
            name_prefix = "intensity_raw"
            
        elif category == "int_conv":
            name_prefix = "intensity_conv"
        
        else:
            raise ValueError("No correct category given")
    
        if filename is None:
            
            name_core = get_name_core(params)

            folder = get_folder(params, class_type)
            
            data = np.loadtxt(os.path.join(mydir.data_dir, folder, name_prefix + name_core + extension))     
                
        elif filename is not None:
            
            data = np.loadtxt(filename)            

        try:
            h, var = data
        except:
            raise ValueError("No correct format for input data from sol_profiles")
                
        profile = lum_profile(radius=h, variable=var, params=params, category=category)
    
    else:
        raise ValueError("No correct class selected, which kind of file are you searching for?")
                
    return profile
    

# CLASSES

class sol_profiles():
    
    def __init__(self, radius, variables, params):
        
        self.r = radius
        
        self.var = variables
        
        self.v = variables[0]
        self.n = variables[1]
        self.T = variables[2]
        
        self.params = params
        
        self.name_prefix = "profiles"
        

    def to_file(self, filename = None, extension = ".dat"):
        
        if filename is None:
                            
            name_core = get_name_core(self.params)
            
            folder = get_folder(self.params, self.name_prefix)
            
            np.savetxt(os.path.join(mydir.data_dir, folder,self.name_prefix + name_core + extension), \
                                (self.r,self.v,self.n,self.T))
            
          
        elif filename is not None:
            
            np.savetxt(filename, (self.r,self.v,self.n,self.T))
      
        return


    
        
    def plot(self, ax=None, savefig=False, size=14, label = None, extension = ".png", **kwargs):
        
        name_core = get_name_core(self.params)
            
        folder = get_folder(self.params, self.name_prefix)
        
        if ax is None:
            
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
        
            ax_v.plot(np.log10(self.r/(1000*nc.pc)),self.v/10**8, label=label, **kwargs)       
            ax_n.plot(np.log10(self.r/(1000*nc.pc)),np.log10(self.n), **kwargs)
            ax_T.plot(np.log10(self.r/(1000*nc.pc)),np.log10(self.T), **kwargs)

            if label is not None:
                
                fig.legend(loc="lower center")
            
            if savefig:
                
                plt.savefig(os.path.join(mydir.plot_dir, folder,self.name_prefix + name_core + extension))
                
                
        elif ax is not None:
            
            ax[0].plot(np.log10(self.r/(1000*nc.pc)),self.v/10**8, label=label, **kwargs)       
            ax[1].plot(np.log10(self.r/(1000*nc.pc)),np.log10(self.n), **kwargs)
            ax[2].plot(np.log10(self.r/(1000*nc.pc)),np.log10(self.T), **kwargs)
                
        return
    
    
    def check_nans(self):
        
        var = np.asarray(self.var)
        
        return np.isnan(np.sum(np.log10(var)))
    
   
    
class ion_profiles():
    
    def __init__(self, radius, variables, params):
        
        self.r = radius
        
        self.var = variables
        
        self.x_e = variables[0]
        self.x_CII = variables[1]
        
        self.params = params

        self.name_prefix = "ionization"
        
               
    def to_file(self, filename = None, extension = ".dat"):
        
        if filename is None:
                            
            name_core = get_name_core(self.params)
            
            folder = get_folder(self.params, self.name_prefix)
            
            np.savetxt(os.path.join(mydir.data_dir, folder, self.name_prefix + name_core + extension), \
                                (self.r,self.x_e,self.x_CII))
            
          
        elif filename is not None:
            
            np.savetxt(filename, (self.r,self.x_e,self.x_CII))
      
        return
    
            
    def plot(self, ax=None, savefig=False, size=14, label = None, extension = ".png", **kwargs):
                
        name_core = get_name_core(self.params)
        
        folder = get_folder(self.params, self.name_prefix)
        
        if ax is None:
            
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
    
            ax_xe.plot(np.log10(self.r/(1000*nc.pc)),1.-self.x_e, label=label, **kwargs)
            ax_xCII.plot(np.log10(self.r/(1000*nc.pc)),self.x_CII, **kwargs)

            if label is not None:
                
                fig.legend(loc="lower center")
            
            if savefig:
                
                plt.savefig(os.path.join(mydir.plot_dir, folder,self.name_prefix + name_core + extension))
                
                
        elif ax is not None:
            
            ax[0].plot(np.log10(self.r/(1000*nc.pc)),1.-self.x_e, label=label, **kwargs)
            ax[1].plot(np.log10(self.r/(1000*nc.pc)),self.x_CII, **kwargs)
                
        return
    
    
    def check_nans(self):
        
        var = np.asarray(self.var)
        
        return np.isnan(np.sum(var))
    

class lum_profile():
    
    def __init__(self, radius, variable, params, category, eta=None):
        """
        category can be sigma, int_raw, int_conv
        """
        
        self.h = radius
        
        self.var = variable
        
        self.params = params
        
        self.category = category
        
        self.eta = eta
        
        if self.category == "sigma":
            self.name_prefix = self.category
            
        elif self.category == "int_raw":
            self.name_prefix = "intensity_raw"
            
        elif self.category == "int_conv":
            self.name_prefix = "intensity_conv"
        
        else:
            raise ValueError("No correct category given")
    
    def to_file(self, filename = None, extension = ".dat"):
        
        if filename is None:
                            
            name_core = get_name_core(self.params)
            
            folder = get_folder(self.params, self.name_prefix)
            
            np.savetxt(os.path.join(mydir.data_dir, folder, self.name_prefix + name_core + extension), \
                                (self.h,self.var))
            
          
        elif filename is not None:
            
            np.savetxt(filename, (self.h,self.var))
      
        return
    
    
    def plot(self, ax=None, data = None, savefig=False, size=14, label = None, extension = ".png", **kwargs):
        
        name_core = get_name_core(self.params)
        
        folder = get_folder(self.params, self.name_prefix)
        
        if ax is None:
            
            fig, ax = plt.subplots(1, 1, sharex=True, figsize=(1.5*8.27,1.5*4.))
    
            ax.set_xlabel("b [kpc]", size=size)
            
            if self.category == "sigma":
                ax.set_ylabel(r"log ($\Sigma_{CII}$ [erg/cm s$^2$])", size=size)
            elif self.category == "int_raw":
                ax.set_ylabel("flux [mJy/arcsec$^2$]", size=size)
            elif self.category == "int_conv":
                ax.set_ylabel("flux [mJy/arcsec$^2$]", size=size)
            else:
                raise ValueError("No correct category")
    
            ax.tick_params(labelsize=size)
            ax.set_xlim((0.,10))
            
            ax.set_yscale("log")
    
            ax.plot(self.h/(1000*nc.pc),self.var, label=label, **kwargs)


            if data is not None:
                ax.errorbar(data.x, data.data, yerr=data.err, \
                            markerfacecolor='C3',markeredgecolor='C3', marker='o',\
                            linestyle='', ecolor = 'C3')
            
            if label is not None:
                
                fig.legend(loc="lower center")
            
            if savefig:
                
                plt.savefig(os.path.join(mydir.plot_dir, folder,self.name_prefix + name_core + extension))
                
        elif ax is not None:
            
            ax.plot(self.h/(1000*nc.pc),self.var, label=label, **kwargs)
                
        return
    
    def check_nans(self):
        
        var = np.asarray(self.var)
        
        return np.isnan(np.sum(var))

    
    
