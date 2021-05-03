# -*- coding: utf-8 -*-
"""
Created on Mon May  3 10:50:09 2021

@author: anna
"""



import os
import numpy as np
import matplotlib.pyplot as plt
import itertools

import natconst as nc
import mydir




load_from_file = False




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

        name = "{}_n{:.1e}_plw{:.1e}_ph1{:.1e}_pg1{:.1e}"\
                .format(name_base, n, Plw, Ph1, Pg1)
                
        
        path = os.path.join(mydir.data_dir, folder, name)

        np.savetxt(path, (T_vec,var_vec))
        
        print("INFO: saving data to {}".format(path))
        return

    
    temperatures = np.logspace(3,8,1000)
    
    densities = np.logspace(-5,3, 17)

    plws = np.logspace(-14,-7, 8)
    ph1s = np.logspace(-14,-7, 8)
    pg1s = np.logspace(-14,-7, 8)
    
    plws = np.insert(plws, 0, 0.)
    ph1s = np.insert(ph1s, 0, 0.)
    pg1s = np.insert(pg1s, 0, 0.)
    
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
        print("ITERATION {} OF {}".format(i, len(n)*len(plws)*len(ph1s)*len(pg1s)))
        

        

if load_from_file == True:
    
    
    pass
