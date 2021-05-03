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




load_from_file = False


temperatures = np.logspace(3,8,1000)

densities = np.logspace(-5,3, 5)

plws = np.logspace(-14,-7, 4)
ph1s = np.logspace(-14,-7, 4)
pg1s = np.logspace(-14,-7, 4)

plws = np.insert(plws, 0, 0.)
ph1s = np.insert(ph1s, 0, 0.)
pg1s = np.insert(pg1s, 0, 0.)


if load_from_file == False:
    
    import gnedincooling as gc

    gc.frtinitcf(0, '/home/elia.pizzati/projects/university/thesis/cf_table.I2.dat')


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
    
        Pc6 = 0.
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
        
        Pc6 = 0.
        Zeta = 1.
    
        return gc.frtgetcf_heat(T, n, Zeta, Plw, Ph1, Pg1, Pc6)


    def cooling_vec(T_vec, n, Plw, Ph1, Pg1): 
        
        cooling_values = [cooling(T, n, Plw, Ph1, Pg1) for T in T_vec]
        
        return np.asarray(cooling_values)
    
    def heating_vec(T_vec, n, Plw, Ph1, Pg1): 
        
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

    
    
    print("###############################################")
    print("STARTING ITERATIONS")      
    print("###############################################")
          
    i=0              
    for n, plw, ph1, pg1 in itertools.product(densities, plws, ph1s, pg1s):
        
        cooling_values = cooling_vec(temperatures, n, plw, ph1, pg1)
        heating_values = heating_vec(temperatures, n, plw, ph1, pg1)
        
        to_file("cooling",cooling_values, temperatures, n, plw, ph1, pg1)
        to_file("heating",heating_values, temperatures, n, plw, ph1, pg1)
        
        i+=1
        print("ITERATION {} OF {}".format(i, densities.size*plws.size*ph1s.size*pg1s.size))
        

        

if load_from_file == True:
    
    
    def from_file(name_base, n, Plw, Ph1, Pg1):
        
        folder = 'data_cooling'
        
        name = "{}_n{:.1e}_plw{:.1e}_ph{:.1e}_pg{:.1e}"\
                .format(name_base, n, Plw, Ph1, Pg1)
                
        path = os.path.join(mydir.data_dir, folder, name)

        T_vec, var_vec = np.loadtxt(path)
        
        print("INFO: loading data from {}".format(path))
        
        return T_vec, var_vec

    
    fig, ax = plt.subplots(figsize=(8.27,5.))

    cmap_rend_col = cm.get_cmap('viridis')
    
    
    
    n_min = min(densities)
    n_max = max(densities)
    
    norm = colors.Normalize(vmin=n_min, vmax=n_max)
   
    
    for n in densities:
        
        T_vec, cooling_vec = from_file("cooling", n, 0.,0.,0.)

        T_vec, heating_vec = from_file("heating", n, 0.,0.,0.)
    
        ax.plot(T_vec, cooling_vec,  color = cmap_rend_col((n-n_min)/(n_max-n_min)))
        ax.plot(T_vec, cooling_vec,  color = cmap_rend_col((n-n_min)/(n_max-n_min)), linestyle=':')


    ax.set_xlabel("T [K]")
    ax.set_ylabel("cooling")
    ax.set_xscale("log")
    ax.set_yscale("log")
    
    #ax.axhline(1., linestyle='--', color='gray')

    #ax.set_title("Eta dependence on n, T; I_UV=0, x_e=0.5, z=5.")    
    
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
    cb.set_label(r'temperature', rotation=90.)

    fig.tight_layout()
