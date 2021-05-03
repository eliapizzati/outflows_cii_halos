# -*- coding: utf-8 -*-
"""
Created on Mon May  3 10:50:09 2021

@author: anna
"""



import os
import numpy as np
import matplotlib.pyplot as plt
import itertools
from matplotlib import ticker, cm, colors


import natconst as nc
import mydir
import plot_config


save_to_file = True
load_from_file = False
plotting = False

temperatures = np.logspace(3,8,1000)

densities = np.logspace(-5,3, 5)

plws = np.logspace(-14,-7, 4)
ph1s = np.logspace(-14,-7, 4)
pg1s = np.logspace(-14,-7, 4)

plws = np.insert(plws, 0, 0.)
ph1s = np.insert(ph1s, 0, 0.)
pg1s = np.insert(pg1s, 0, 0.)



def cooling(T, n, Plw, Ph1, Pg1):
        """
        Cooling function (calling gnedincooling from Fortran functions)
        
        Parameters
        ==========
        T: float
            temperature
        n: float
            density
        Plw: float
            Lyman-Wernner PD rate
        Ph1: float
            H photoionization rate
        Pg1: float
            He photoionization rate
        """    
    
        Pc6 = 5e-17.
        Zeta = 1.
    
        return gc.frtgetcf_cool(T, n, Zeta, Plw, Ph1, Pg1, Pc6)


def heating(T, n, Plw, Ph1, Pg1):
        """
        Heating function (calling gnedincooling from Fortran functions)
        
        Parameters
        ==========
        T: float
            temperature
        n: float
            density
        Plw: float
            Lyman-Wernner PD rate
        Ph1: float
            H photoionization rate
        Pg1: float
            He photoionization rate    
        """    
        
        Pc6 = 5e-17.
        Zeta = 1.
    
        return gc.frtgetcf_heat(T, n, Zeta, Plw, Ph1, Pg1, Pc6)


def cooling_vectorized(T_vec, n, Plw, Ph1, Pg1): 
        
        cooling_values = [cooling(T, n, Plw, Ph1, Pg1) for T in T_vec]
        
        return np.asarray(cooling_values)
    
def heating_vectorized(T_vec, n, Plw, Ph1, Pg1): 
        
        heating_values = [heating(T, n, Plw, Ph1, Pg1) for T in T_vec]
        
        return np.asarray(heating_values)
    

    
def to_file(name_base, var_vec, T_vec, n, Plw, Ph1, Pg1):
        
        folder = 'data_cooling'
        
        if not os.path.exists(os.path.join(mydir.data_dir, folder)):
            os.mkdir(os.path.join(mydir.data_dir, folder))

        name = "{}_n{:.1e}_plw{:.1e}_ph{:.1e}_pg{:.1e}"\
                .format(name_base, n, Plw, Ph1, Pg1)
                
        
        path = os.path.join(mydir.data_dir, folder, name)

        np.savetxt(path, (T_vec,var_vec))
        
        print("INFO: saving data to {}".format(path))
        return

def from_file(name_base, n, Plw, Ph1, Pg1):
        
        folder = 'data_cooling'
        
        name = "{}_n{:.1e}_plw{:.1e}_ph{:.1e}_pg{:.1e}"\
                .format(name_base, n, Plw, Ph1, Pg1)
                
        path = os.path.join(mydir.data_dir, folder, name)

        T_vec, var_vec = np.loadtxt(path)
        
        print("INFO: loading data from {}".format(path))
        
        return T_vec, var_vec
    
    
    


if save_to_file == True:
    
    import gnedincooling as gc

    gc.frtinitcf(0, os.path.join(mydir.script_dir, "input_data", "cf_table.I2.dat"))
    
    print("###############################################")
    print("STARTING ITERATIONS")      
    print("###############################################")
          
    i=0              
    for n, plw, ph1, pg1 in itertools.product(densities, plws, ph1s, pg1s):
        
        cooling_values = cooling_vectorized(temperatures, n, plw, ph1, pg1)
        heating_values = heating_vectorized(temperatures, n, plw, ph1, pg1)
        
        to_file("cooling",cooling_values, temperatures, n, plw, ph1, pg1)
        to_file("heating",heating_values, temperatures, n, plw, ph1, pg1)
        
        i+=1
        print("ITERATION {} OF {}".format(i, densities.size*plws.size*ph1s.size*pg1s.size))
        

if plotting:
        
    # FIRST PLOT (density)
            
            
    fig, ax = plt.subplots(figsize=(8.27,5.))
    
    cmap_rend_col = cm.get_cmap('viridis')
    
    plw = plws[-1]
    ph1 = ph1s[0]
    pg1 = pg1s[0]
    
    n_min = min(densities)
    n_max = max(densities)
    
    norm = colors.LogNorm(vmin=n_min/10**0.5, vmax=n_max*10**0.5)
       
    
    for n in densities:
    
        if load_from_file == True:        
            T_vec, cooling_vec = from_file("cooling", n, plw,ph1,pg1)
            T_vec, heating_vec = from_file("heating", n, plw,ph1,pg1)
        else:
            T_vec, cooling_vec = cooling_vectorized(temperatures,n, plw,ph1,pg1)
            T_vec, heating_vec = heating_vectorized(temperatures,n, plw,ph1,pg1)
           
        ax.plot(T_vec, cooling_vec,  color = cmap_rend_col((np.log10(n)-np.log10(n_min))/(np.log10(n_max)-np.log10(n_min))))
        ax.plot(T_vec, heating_vec,   color = cmap_rend_col((np.log10(n)-np.log10(n_min))/(np.log10(n_max)-np.log10(n_min))),\
                linestyle=':')
    
    
    ax.set_xlabel("T [K]")
    ax.set_ylabel(r"$|\Lambda$(T)| $[\mathrm{erg}\,\mathrm{cm}^3\,\mathrm{s}^{-1}]$")
    
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim((1e-25,1e-21))
    
    #ax.axhline(1., linestyle='--', color='gray')
    
    ax.set_title("Plw={:.1e}, Ph1={:.1e}, Pg1={:.1e}".format(plw,ph1,pg1))    
    
    plt.subplots_adjust(left = 0.15,  # the left side of the subplots of the figure
        right = 0.98,   # the right side of the subplots of the figure
        bottom = 0.15,  # the bottom of the subplots of the figure
        top = 0.9,     # the top of the subplots of the figure
        wspace = 0.1,  # the amount of width reserved for space between subplots,
        # expressed as a fraction of the average axis width
        hspace = 0.1)  # the amount of height reserved for space between subplots,
                      # expressed as a fraction of the average axis height
    
    
    cmap = cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, orientation='vertical')
    cb.set_label(r'density [cm^-3]', rotation=90.)
    
    fig.tight_layout()
    
    
    # SECOND PLOT (plw)
    
    
    fig, ax = plt.subplots(figsize=(8.27,5.))
    
    cmap_rend_col = cm.get_cmap('plasma')
    
    n = densities[3]
    ph1 = ph1s[0]
    pg1 = pg1s[0]
    
    plw_min = min(plws[1::])
    plw_max = max(plws[1::])
    
    norm = colors.LogNorm(vmin=plw_min/10**0.5, vmax=plw_max*10**0.5)
       
    
    for plw in plws[1::]:
        print(plw)
    
        if load_from_file == True:        
            T_vec, cooling_vec = from_file("cooling", n, plw,ph1,pg1)
            T_vec, heating_vec = from_file("heating", n, plw,ph1,pg1)
        else:
            T_vec, cooling_vec = cooling_vectorized(temperatures,n, plw,ph1,pg1)
            T_vec, heating_vec = heating_vectorized(temperatures,n, plw,ph1,pg1)
           
        ax.plot(T_vec, cooling_vec,  color = cmap_rend_col((np.log10(plw)-np.log10(plw_min))/(np.log10(plw_max)-np.log10(plw_min))))
        ax.plot(T_vec, heating_vec,   color = cmap_rend_col((np.log10(plw)-np.log10(plw_min))/(np.log10(plw_max)-np.log10(plw_min))),\
                linestyle=':')
    
    
    if load_from_file == True:        
            T_vec, cooling_vec = from_file("cooling", n, plws[0],ph1,pg1)
            T_vec, heating_vec = from_file("heating", n, plws[0],ph1,pg1)
    else:
            T_vec, cooling_vec = cooling_vectorized(temperatures,n, plws[0],ph1,pg1)
            T_vec, heating_vec = heating_vectorized(temperatures,n, plws[0],ph1,pg1)
      
    
    ax.plot(T_vec, cooling_vec, linestyle='-', color='gray')
    ax.plot(T_vec, heating_vec, linestyle=':', color="gray")
    
    ax.set_xlabel("T [K]")
    ax.set_ylabel(r"$|\Lambda$(T)| $[\mathrm{erg}\,\mathrm{cm}^3\,\mathrm{s}^{-1}]$")
    
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim((1e-25,1e-21))
    
    
    ax.set_title("n={:.1e}, Ph1={:.1e}, Pg1={:.1e}".format(n,ph1,pg1))    
    
    plt.subplots_adjust(left = 0.15,  # the left side of the subplots of the figure
        right = 0.98,   # the right side of the subplots of the figure
        bottom = 0.15,  # the bottom of the subplots of the figure
        top = 0.9,     # the top of the subplots of the figure
        wspace = 0.1,  # the amount of width reserved for space between subplots,
        # expressed as a fraction of the average axis width
        hspace = 0.1)  # the amount of height reserved for space between subplots,
                      # expressed as a fraction of the average axis height
    
    
    cmap = cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, orientation='vertical')
    cb.set_label(r'Plw [s^-1]', rotation=90.)
    
    fig.tight_layout()
    
    
    
    # THIRD PLOT (ph1)
    
    
    fig, ax = plt.subplots(figsize=(8.27,5.))
    
    cmap_rend_col = cm.get_cmap('plasma')
    
    n = densities[0]
    plw = plws[2]
    pg1 = pg1s[0]
    
    ph1_min = min(ph1s[1::])
    ph1_max = max(ph1s[1::])
    
    norm = colors.LogNorm(vmin=ph1_min/10**0.5, vmax=ph1_max*10**0.5)
       
    
    for ph1 in ph1s[1::]:
    
        if load_from_file == True:        
            T_vec, cooling_vec = from_file("cooling", n, plw,ph1,pg1)
            T_vec, heating_vec = from_file("heating", n, plw,ph1,pg1)
        else:
            T_vec, cooling_vec = cooling_vectorized(temperatures,n, plw,ph1,pg1)
            T_vec, heating_vec = heating_vectorized(temperatures,n, plw,ph1,pg1)
           
        ax.plot(T_vec, cooling_vec,  color = cmap_rend_col((np.log10(ph1)-np.log10(ph1_min))/(np.log10(ph1_max)-np.log10(ph1_min))))
        ax.plot(T_vec, heating_vec,   color = cmap_rend_col((np.log10(ph1)-np.log10(ph1_min))/(np.log10(ph1_max)-np.log10(ph1_min))),\
                linestyle=':')
    
    
    if load_from_file == True:        
            T_vec, cooling_vec = from_file("cooling", n, plw,ph1s[0],pg1)
            T_vec, heating_vec = from_file("heating", n, plw,ph1s[0],pg1)
    else:
            T_vec, cooling_vec = cooling_vectorized(temperatures,n, plw,ph1s[0],pg1)
            T_vec, heating_vec = heating_vectorized(temperatures,n, plw,ph1s[0],pg1)
      
    
    #ax.plot(T_vec, cooling_vec, linestyle='-', color='gray')
    #ax.plot(T_vec, heating_vec, linestyle=':', color="gray")
    
    ax.set_xlabel("T [K]")
    ax.set_ylabel(r"$|\Lambda$(T)| $[\mathrm{erg}\,\mathrm{cm}^3\,\mathrm{s}^{-1}]$")
    
    ax.set_xscale("log")
    ax.set_yscale("log")
    #ax.set_ylim((1e-25,1e-21))
    
    
    ax.set_title("n={:.1e}, Plw={:.1e}, Pg1={:.1e}".format(n,plw,pg1))    
    
    plt.subplots_adjust(left = 0.15,  # the left side of the subplots of the figure
        right = 0.98,   # the right side of the subplots of the figure
        bottom = 0.15,  # the bottom of the subplots of the figure
        top = 0.9,     # the top of the subplots of the figure
        wspace = 0.1,  # the amount of width reserved for space between subplots,
        # expressed as a fraction of the average axis width
        hspace = 0.1)  # the amount of height reserved for space between subplots,
                      # expressed as a fraction of the average axis height
    
    
    cmap = cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, orientation='vertical')
    cb.set_label(r'Ph1 [s^-1]', rotation=90.)
    
    fig.tight_layout()
    
    
    # SECOND PLOT (pg1)
    
    
    fig, ax = plt.subplots(figsize=(8.27,5.))
    
    cmap_rend_col = cm.get_cmap('plasma')
    
    n = densities[3]
    ph1 = ph1s[0]
    pg1 = pg1s[0]
    
    plw_min = min(plws[1::])
    plw_max = max(plws[1::])
    
    norm = colors.LogNorm(vmin=plw_min/10**0.5, vmax=plw_max*10**0.5)
       
    
    for plw in plws[1::]:
        print(plw)
    
        if load_from_file == True:        
            T_vec, cooling_vec = from_file("cooling", n, plw,ph1,pg1)
            T_vec, heating_vec = from_file("heating", n, plw,ph1,pg1)
        else:
            T_vec, cooling_vec = cooling_vectorized(temperatures,n, plw,ph1,pg1)
            T_vec, heating_vec = heating_vectorized(temperatures,n, plw,ph1,pg1)
           
        ax.plot(T_vec, cooling_vec,  color = cmap_rend_col((np.log10(plw)-np.log10(plw_min))/(np.log10(plw_max)-np.log10(plw_min))))
        ax.plot(T_vec, heating_vec,   color = cmap_rend_col((np.log10(plw)-np.log10(plw_min))/(np.log10(plw_max)-np.log10(plw_min))),\
                linestyle=':')
    
    
    if load_from_file == True:        
            T_vec, cooling_vec = from_file("cooling", n, plws[0],ph1,pg1)
            T_vec, heating_vec = from_file("heating", n, plws[0],ph1,pg1)
    else:
            T_vec, cooling_vec = cooling_vectorized(temperatures,n, plws[0],ph1,pg1)
            T_vec, heating_vec = heating_vectorized(temperatures,n, plws[0],ph1,pg1)
      
    
    ax.plot(T_vec, cooling_vec, linestyle='-', color='gray')
    ax.plot(T_vec, heating_vec, linestyle=':', color="gray")
    
    ax.set_xlabel("T [K]")
    ax.set_ylabel(r"$|\Lambda$(T)| $[\mathrm{erg}\,\mathrm{cm}^3\,\mathrm{s}^{-1}]$")
    
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim((1e-25,1e-21))
    
    
    ax.set_title("n={:.1e}, Ph1={:.1e}, Pg1={:.1e}".format(n,ph1,pg1))    
    
    plt.subplots_adjust(left = 0.15,  # the left side of the subplots of the figure
        right = 0.98,   # the right side of the subplots of the figure
        bottom = 0.15,  # the bottom of the subplots of the figure
        top = 0.9,     # the top of the subplots of the figure
        wspace = 0.1,  # the amount of width reserved for space between subplots,
        # expressed as a fraction of the average axis width
        hspace = 0.1)  # the amount of height reserved for space between subplots,
                      # expressed as a fraction of the average axis height
    
    
    cmap = cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, orientation='vertical')
    cb.set_label(r'Plw [s^-1]', rotation=90.)
    
    fig.tight_layout()
    
    
    
    
