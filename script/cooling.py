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


save_to_file = False
save_plot = True
load_from_file = False
plotting = True

values = np.asarray([1e1, 1e-13, 1e-13, 1e-13, 0.])

temperatures = np.logspace(3,8,1000)

densities = np.logspace(-5,3, 9)

plws = np.logspace(-14,-7, 8)
ph1s = np.logspace(-14,-7, 8)
pg1s = np.logspace(-14,-7, 8)

pc6s = np.array([0.])

plws = np.insert(plws, 0, 0.)
ph1s = np.insert(ph1s, 0, 0.)
pg1s = np.insert(pg1s, 0, 0.)



def cooling(T, n, Plw, Ph1, Pg1, Pc6=0.):
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
    
        Zeta = 1.
    
        return gc.frtgetcf_cool(T, n, Zeta, Plw, Ph1, Pg1, Pc6)


def heating(T, n, Plw, Ph1, Pg1, Pc6=0.):
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
        
        Zeta = 1.
    
        return gc.frtgetcf_heat(T, n, Zeta, Plw, Ph1, Pg1, Pc6)


def cooling_vectorized(T_vec, n, Plw, Ph1, Pg1, Pc6): 
        
        cooling_values = [cooling(T, n, Plw, Ph1, Pg1, Pc6) for T in T_vec]
        
        return np.asarray(cooling_values)
    
def heating_vectorized(T_vec, n, Plw, Ph1, Pg1, Pc6): 
        
        heating_values = [heating(T, n, Plw, Ph1, Pg1, Pc6) for T in T_vec]
        
        return np.asarray(heating_values)
    

    
def to_file(name_base, var_vec, T_vec, n, Plw, Ph1, Pg1, Pc6):
        
        folder = 'data_cooling'
        
        if not os.path.exists(os.path.join(mydir.data_dir, folder)):
            os.mkdir(os.path.join(mydir.data_dir, folder))

        name = "{}_n{:.1e}_plw{:.1e}_ph{:.1e}_pg{:.1e}_pc{:.1e}.dat"\
                .format(name_base, n, Plw, Ph1, Pg1, Pc6)
                
        
        path = os.path.join(mydir.data_dir, folder, name)

        np.savetxt(path, (T_vec,var_vec))
        
        print("INFO: saving data to {}".format(path))
        return

def from_file(name_base, n, Plw, Ph1, Pg1, Pc6):
        
        folder = 'data_cooling'
        
        name = "{}_n{:.1e}_plw{:.1e}_ph{:.1e}_pg{:.1e}_pc{:.1e}.dat"\
                .format(name_base, n, Plw, Ph1, Pg1, Pc6)
                
        path = os.path.join(mydir.data_dir, folder, name)

        T_vec, var_vec = np.loadtxt(path)
        
        print("INFO: loading data from {}".format(path))
        
        return T_vec, var_vec
    
    
    

def plot_quantity(quantity, quantity_array, values):
    
    if quantity == "n":
        n, plw, ph1, pg1, pc6 = values
        title = "Plw={:.1e}, Ph1={:.1e}, Pg1={:.1e}, Pc6={:.1e}".format(plw,ph1,pg1, pc6)
        label_cmap = r'density [cm^-3]'
        name = "plot_lamda_vs_T_vs_{}_plw{:.1e}_ph{:.1e}_pg{:.1e}_pc{:.1e}.png"\
                .format(quantity, plw, ph1, pg1, pc6)
                
    elif quantity == "plw":
        n, plw, ph1, pg1, pc6 = values
        title = "n={:.1e}, Ph1={:.1e}, Pg1={:.1e}, Pc6={:.1e}".format(n,ph1,pg1, pc6)
        label_cmap = r'plw [s^-1]'
        name = "plot_lamda_vs_T_vs_{}_n{:.1e}_ph{:.1e}_pg{:.1e}_pc{:.1e}.png"\
                .format(quantity, n, ph1, pg1, pc6)

    elif quantity == "ph1":
        n, plw, ph1, pg1, pc6 = values
        title = "n={:.1e}, Plw={:.1e}, Pg1={:.1e}, Pc6={:.1e}".format(n, plw, pg1, pc6)
        label_cmap = r'ph1 [s^-1]'
        name = "plot_lamda_vs_T_vs_{}_n{:.1e}_plw{:.1e}_pg{:.1e}_pc{:.1e}.png"\
                .format(quantity, n, plw, pg1, pc6)

    elif quantity == "pg1":
        n, plw, ph1, pg1, pc6 = values
        title = "n={:.1e}, Plw={:.1e}, Ph1={:.1e}, Pc6={:.1e}".format(n, plw, ph1, pc6 )
        label_cmap = r'pg1 [s^-1]'
        name = "plot_lamda_vs_T_vs_{}_n{:.1e}_plw{:.1e}_ph{:.1e}_pc{:.1e}.png"\
                .format(quantity, n, plw, ph1, pc6 )

    elif quantity == "pc6":
        n, plw, ph1, pg1, pc6 = values
        title = "n={:.1e}, Plw={:.1e}, Ph1={:.1e}, Pg1={:.1e}".format(n, plw, ph1, pg1 )
        label_cmap = r'pc6 [s^-1]'
        name = "plot_lamda_vs_T_vs_{}_n{:.1e}_plw{:.1e}_ph{:.1e}_pg{:.1e}.png"\
                .format(quantity, n, plw, ph1, pg1 )


    fig, ax = plt.subplots(figsize=(8.27,5.))
    
    cmap_rend_col = cm.get_cmap('viridis')
    
        
    v_min = min(quantity_array)
    v_max = max(quantity_array)
    
    norm = colors.LogNorm(vmin=v_min/10**0.5, vmax=v_max*10**0.5)
       
    for q in quantity_array:

        if quantity == "n":
            n = q
        elif quantity == "plw":
            plw = q
        elif quantity == "ph1":
            ph1 = q
        elif quantity == "pg1":
            pg1 = q
        elif quantity == "pg1":
            pc6 = q

            
        if load_from_file == True:        
            T_vec, cooling_vec = from_file("cooling", n, plw,ph1,pg1, pc6)
            T_vec, heating_vec = from_file("heating", n, plw,ph1,pg1, pc6)
        else:
            
            T_vec, cooling_vec = cooling_vectorized(temperatures,n, plw,ph1,pg1,pc6)
            T_vec, heating_vec = heating_vectorized(temperatures,n, plw,ph1,pg1,pc6)
           
        ax.plot(T_vec, cooling_vec,  color = cmap_rend_col((np.log10(q)-np.log10(v_min))/(np.log10(v_max)-np.log10(v_min))))
        ax.plot(T_vec, heating_vec,   color = cmap_rend_col((np.log10(q)-np.log10(v_min))/(np.log10(v_max)-np.log10(v_min))),\
                linestyle=':')
    
    
    ax.set_xlabel("T [K]")
    ax.set_ylabel(r"$|\Lambda$(T)| $[\mathrm{erg}\,\mathrm{cm}^3\,\mathrm{s}^{-1}]$")
    
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim((1e-25,1e-21))
        
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
    
    ax.set_title(title)    
    cb.set_label(label_cmap, rotation=90.)
    
    fig.tight_layout()
    
    if save_plot:
        folder = 'plot_cooling'
        
        if not os.path.exists(os.path.join(mydir.plot_dir, folder)):
            os.mkdir(os.path.join(mydir.plot_dir, folder))

        path = os.path.join(mydir.plot_dir, folder,  name)

        plt.savefig(name)
    

    
if save_to_file == True:

    import gnedincooling as gc
    
    gc.frtinitcf(0, os.path.join(mydir.script_dir, "input_data", "cf_table.I2.dat"))
    
    print("###############################################")
    print("STARTING ITERATIONS")      
    print("###############################################")
          
    i=0           
    
    for n, plw, ph1, pg1, pc6 in itertools.product(densities, plws, ph1s, pg1s, pc6s):
        
        cooling_values = cooling_vectorized(temperatures, n, plw, ph1, pg1, pc6)
        heating_values = heating_vectorized(temperatures, n, plw, ph1, pg1, pc6)
        
        to_file("cooling",cooling_values, temperatures, n, plw, ph1, pg1, pc6)
        to_file("heating",heating_values, temperatures, n, plw, ph1, pg1, pc6)
        
        i+=1
        print("ITERATION {} OF {}".format(i, densities.size*plws.size*ph1s.size*pg1s.size))
    




if plotting:
        
    if load_from_file == False:
        
        import gnedincooling as gc

        gc.frtinitcf(0, os.path.join(mydir.script_dir, "input_data", "cf_table.I2.dat"))            
            
    
    quantity = input("which quantity (n,plw,ph1,pg1,pc6)?")   
    
        
    if quantity == "n":
        quantity_array = densities
        other_values = np.delete(values, 0)
        
    elif quantity == "plw":
        quantity_array = plws
        other_values = np.delete(values, 1)

    elif quantity == "ph1":
        quantity_array = ph1s
        other_values = np.delete(values, 2)

    elif quantity == "pg1":
        quantity_array = pg1s
        other_values = np.delete(values, 3)

    elif quantity == "pc6":
        quantity_array = pc6s
        other_values = np.delete(values, 4)


    
    plot_quantity(quantity, quantity_array, values = other_values)
    

    
    
    
    
