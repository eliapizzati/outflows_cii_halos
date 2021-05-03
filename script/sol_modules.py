"""
Created on Fri Mar  5 16:12:54 2021

@author: anna
"""

import os
import natconst as nc
import mydir 
import time


import scipy.integrate as si

import numpy as np

from model_classes import sol_profiles
from radiation_fields import UVB_rates

from my_utils import get_circular_velocity_profile_NFW, get_virial_mass_from_vc

# FUNCTIONS
        
# cooling function


import gnedincooling as gc

gc.frtinitcf(0, os.path.join(mydir.script_dir, "input_data", "cf_table.I2.dat"))

def lamda(T, n, r, params, Plw, Ph1, Pg1):
    """
    Cooling function (calling gnedincooling from Fortran functions)
    
    Parameters
    ==========
    T: float
        temperature
    n: float
        density
    r: float
        radius
    params: dict
        parameters to be passed to all functions

    """    
    
    # params definition
     
    if "SFR" in params:
        SFR_pure = params["SFR"]
    else: 
        raise ValueError("No SFR given")
    
    if "f_esc_ion" in params:
        f_esc_ion = params["f_esc_ion"]
    else:
        f_esc_ion = 0.
        
    if "f_esc_FUV" in params:
        f_esc_FUV = params["f_esc_FUV"]
    else:
        f_esc_FUV = 0.
        
    if "Zeta" in params:
        Zeta = params["Zeta"]
    else:
        Zeta = 1.
    
    # loading photoionization rates
    #
    Pc6 = 0.

    # UV background (values around 1e-13 s^-1)
        
    Plw = Plw
    Ph1 = Ph1
    Pg1 = Pg1            
        
    # ADDING THE GALACTIC FLUX CONTRIBUTION 
        
    if f_esc_ion != 0.0:
        Ph1 += nc.Gamma_H_1000 * (1000*nc.pc/r)**2 * SFR_pure*f_esc_ion
        Pg1 += nc.Gamma_He_1000 * (1000*nc.pc/r)**2 * SFR_pure*f_esc_ion

    if f_esc_FUV != 0.0:
        Plw += nc.Gamma_LW_1000 * (1000*nc.pc/r)**2  * SFR_pure*f_esc_FUV
    
    return gc.frtgetcf_cool(T, n, Zeta, Plw, Ph1, Pg1, Pc6) - gc.frtgetcf_heat(T, n, Zeta, Plw, Ph1, Pg1, Pc6)



# SYSTEM OF EQUATIONS

def diff_system(r, y, params, Plw, Ph1, Pg1):
    """
    differential system for the Euler equations
    
    Parameters
    ==========
    r: float
        temperature
    y: array
        array of the variables (v, n, T)
    params: dict
        parameters to be passed to all functions

    Returns
    =======
    output from the differential system: array

    """    
    
    # params definition
     
    if params["DM_model"] == "NFW":    

        if "redshift" in params:
            redshift = params["redshift"]
        else: 
            raise ValueError("No redshift given")
        
        if "M_vir" in params:
            M_vir_pure = params["M_vir"]
        else: 
            if "v_c" in params:
                v_c_pure = params["v_c"]
            else: 
                raise ValueError("No v_c given")
            
            M_vir_pure = get_virial_mass_from_vc(v_c_pure*1e5, redshift)
            
    elif params["DM_model"] == "iso_sphere" or params["DM_model"] is None:

            if "v_c" in params:
                v_c_pure = params["v_c"]
            else: 
                raise ValueError("No v_c given")
                
    else:
        raise ValueError("No correct DM model; DM models are: NFW, iso_sphere (=None)")
        
    # defining the equations 
    
    assert len(y.shape) == 1, 'Dimensionality of the input array is wrong'

    v = y[0] # velocity
    rho = y[1] # density
    T = y[2] # temperature
        
    c_S2 = nc.gamma*nc.knorm * T
    c_T2 = nc.knorm * T
    
    n = rho / (nc.mus * nc.mp) #in cgs unit # NB CHECK CONSISTENCY WITH METALLICITY
    
    q = rho*lamda(T, n, r, params, Plw, Ph1, Pg1) / (nc.mus * nc.mp)**2
    
    if params["DM_model"] == "NFW":    
        
        v_c = get_circular_velocity_profile_NFW(r, M_vir_pure, redshift)
                
    elif params["DM_model"] == "iso_sphere" or params["DM_model"] is None:

        v_c = v_c_pure * 1e5 #cm/s 
    
    v_e = v_c * np.sqrt(2) #cm/s 

    
    output_a = ((2*v/r)*(c_S2 - v_e**2/4) + (nc.gamma-1)*q) / (v**2-c_S2)
    
    output_b = ((2*rho/r)*(v_e**2/4-v**2) - (nc.gamma-1)*q*rho/v) / (v**2-c_S2)

    output_c = ((nc.gamma-1)*(2*T/r)*(v_e**2/4-v**2) - (nc.gamma-1)*q*(v**2-c_T2)/(nc.knorm*v)) / (v**2-c_S2)        
    
    output = np.asarray([output_a,output_b,output_c])
    
    assert len(output.shape) == 1  , 'Dimensionality of the output array is wrong'
    
    return output



def get_profiles(params, resol=1000, print_time=False, integrator="RK45"):
    """
    computes the profiles for v, n, T as a function of r
    
    Parameters
    ==========
    params: dict
        parameters to be passed to all functions
    resol: int, optional
        number of r-steps
    
    Returns
    =======
    profiles: sol_profiles class element

    """    
    
    # params definition
    
    if "beta" in params:
        beta = params["beta"]
    else: 
        raise ValueError("No beta given")
    
    if "SFR" in params:
        SFR_pure = params["SFR"]
    else: 
        raise ValueError("No SFR given")
        
    if "redshift" in params:
        redshift = params["redshift"]
    else: 
        raise ValueError("No redshift given")

    if "Zeta" in params:
        Zeta = params["Zeta"]
    else:
        Zeta = 1.
    
    if "alfa" in params:
        alfa = params["alfa"]
    else:
        alfa = 1.
        
    if "R_in" in params:
        R_in_pure = params["R_in"]
    else:
        alfa = 0.3
    
    # getting the BC
    
    Plw = UVB_rates(redshift, quantity="LW rate")
    Ph1 = UVB_rates(redshift, quantity="H rate")
    Pg1 = UVB_rates(redshift, quantity="He rate")  
    
    
    SFR = SFR_pure/nc.year #1/s
    
    E_SN = 1e51*(SFR)/100 #erg (energy from supernovae-- energy of a single SN \
                        #per 100 M_sun of Star Formation)

    M_dot = beta*SFR*nc.ms  #mass from SN
    
    E_dot = alfa*E_SN #erg/s
  
    #M0 = 1. 
    v0 = np.sqrt(E_dot/M_dot)/np.sqrt(2)  #m/s
    #c0 = v0/M0

    R = R_in_pure*1000*nc.pc #cm

    rho0 = 0.1125395*np.sqrt(M_dot**3/E_dot) / R**2 #g/cm^3
    P0 = 0.0337618*np.sqrt(M_dot*E_dot) / R**2  #g/cm^2    
    T0 = P0/(rho0*nc.knorm)  #K
    
    y0 = np.asarray([v0,rho0,T0]) 
    
    # integrating the equations
    
    r_bound = (R, 100*R)
    
    r_eval = np.linspace(r_bound[0],r_bound[1],resol)

    if print_time:
      t_ivp = time.perf_counter()

    sol = si.solve_ivp(diff_system, r_bound, y0, t_eval=r_eval, args=(params,Plw,Ph1,Pg1),method = integrator) #,rtol=1.0e-3
    
    if print_time:
      time_ivp = (time.perf_counter() - t_ivp)
      print("total time ivp (s)=", time_ivp)

    if sol.success == False:
        print('Integration stopped before reaching the end of the array')
    elif sol.success == True:
        print('Integration completed successfully')
    
    r = sol.t #cm
    v = sol.y[0] #cm/s
    rho = sol.y[1]
    n = rho / (nc.mus*nc.mp) #cm^-3 # NB CHECK CONSISTENCY WITH METALLICITY
    T = sol.y[2] #K
    
    mask = v > 4.1e6

    profiles = []
    
    profiles.append(v[mask])
    profiles.append(n[mask])
    profiles.append(T[mask])
    
    return sol_profiles(radius=r[mask], variables=profiles, params=params)





        
