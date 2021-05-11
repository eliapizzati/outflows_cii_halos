# -*- coding: utf-8 -*-
"""
Created on Tue May 11 11:04:50 2021

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
from my_utils import get_virial_mass_from_vc, get_concentration, get_virial_radius
from sol_modules import get_profiles

import plot_config as pltc

# FUNCTIONS
        
# cooling function


import gnedincooling as gc

gc.frtinitcf(0, os.path.join(mydir.script_dir, "input_data", "cf_table.I2.dat"))        


# SYSTEM OF EQUATIONS

def diff_system_fast(r, y, SFR_pure, redshift, M_vir_pure, f_esc_ion, f_esc_FUV, Plw, Ph1, Pg1, Zeta, A_NFW, r_s):
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
        
    # loading constants
    kk = 1.3807e-16  # Bolzmann's constant [erg/K]
    mp = 1.6726e-24  # Mass of proton [g]
    gamma = 5./3
    mus = 0.61 # mean molecular weight of the sun
    pc = 3.08572e18   # Parsec [cm]
    ms = 1.99e33      # Solar mass [g]
    gg = 6.672e-8   # Gravitational constant

    Gamma_H_1000 = 5.48031935502901e-09 # s^-1
    Gamma_He_1000 = 1.7687762344020628e-09 # s^-1
    Gamma_LW_1000 = 1.4229125141877616e-08 # s^-1

    # defining the equations 
    
    v = y[0] # velocity in kms
    n = y[1] # density in cm-3
    T = y[2] # temperature in K
        
    knorm_kmsK = kk/(mus*mp) / 1e10 # in (km/s)**2 / K 
    print("##############################")

    c_S2 = gamma*knorm_kmsK * abs(T) # in (km/s)^2
    print("cs", np.sqrt(c_S2))
    c_T2 = knorm_kmsK * abs(T)    # in (km/s)^2
    
    # cooling part
        
    Plw = Plw + Gamma_LW_1000 * (1./r)**2  * SFR_pure*f_esc_FUV
    Ph1 = Ph1 + Gamma_H_1000 * (1./r)**2 * SFR_pure*f_esc_ion
    Pg1 = Pg1 + Gamma_He_1000 * (1./r)**2 * SFR_pure*f_esc_ion
          
    Pc6 = 0.
    
    lamda = gc.frtgetcf_cool(T, n, Zeta, Plw, Ph1, Pg1, Pc6) - gc.frtgetcf_heat(T, n, Zeta, Plw, Ph1, Pg1, Pc6)

    q = n * lamda / (mus * mp) / 1e10 # in (km/s)^2/s
    print("q", q)
    
    M_r = M_vir_pure/A_NFW * (np.log(1.+r/r_s)+r_s/(r_s+r) - 1)

    v_c = np.sqrt(gg*M_r*ms/(r*1e3*pc)) / 1e5 #in km/s
    print("v_c", v_c)
    v_e = v_c * np.sqrt(2)  # in km/s
    print("slow", )

    # final system
    print("v", v)
    print("T", T)
    print("n", n)
    print("##############################")

    derivative_a = (2*v/r)*(c_S2 - v_e**2/4) / (v**2-c_S2) + (gamma-1.) * q / ((v**2-c_S2)*1e2/pc)
    
    derivative_b = (2*n/r)*(v_e**2/4 - v**2) / (v**2-c_S2) - n / (1e2*v/pc) * (gamma-1.) * q / (v**2-c_S2)

    derivative_c = (gamma-1.) * (2*T/r) * (v_e**2/4-v**2) / (v**2-c_S2) - (gamma-1.)*q / (knorm_kmsK*(1e2*v/pc)) * (v**2-c_T2) / (v**2-c_S2)        
    
    output = np.asarray([derivative_a,derivative_b,derivative_c])
        
    return output


def stopping_condition(t, y, SFR_pure, redshift, M_vir_pure, f_esc_ion, f_esc_FUV, Plw, Ph1, Pg1, Zeta, A_NFW, r_s): 
    return y[0] - 50. # in km/s

stopping_condition.terminal = True
stopping_condition.direction = -1


def get_profiles_fast(params, resol=1000, print_time=False, integrator="RK45"):
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
        
    if "redshift" in params:
        redshift = params["redshift"]
    else: 
        raise ValueError("No redshift given")
    
    if "alfa" in params:
        alfa = params["alfa"]
    else:
        alfa = 1.
        
    if "R_in" in params:
        R_in_pure = params["R_in"]
    else:
        alfa = 0.3
        
    if params["DM_model"] == "NFW":    
        
        if "M_vir" in params:
            M_vir_pure = params["M_vir"]
        else: 
            if "v_c" in params:
                v_c_pure = params["v_c"]
            else: 
                raise ValueError("No v_c given")
            
            M_vir_pure = get_virial_mass_from_vc(v_c_pure*1e5, redshift)


    # getting some preliminar values 
    
    # gravity part
        
    c = get_concentration(M_vir_pure, redshift)

    A_NFW = np.log(1+c) - c/(1.+c)
    
    
    r_s = get_virial_radius(M_vir_pure, redshift) / c / 1e3 / nc.pc # in kpc

    
    # UVB part
    
    Plw = UVB_rates(redshift, quantity="LW rate")
    Ph1 = UVB_rates(redshift, quantity="H rate")
    Pg1 = UVB_rates(redshift, quantity="He rate")  
    
    # getting the BC
    
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
    
    # changing the dimensions for the integration part
    
    R_kpc = R_in_pure # in kpc
    
    v0_kms = v0 / 1e5 # in km/s
    n0_cm3 = rho0 / (nc.mp*nc.mus) # in cm-3
    T0_K = T0 # in K
    
    y0 = np.asarray([v0_kms,n0_cm3,T0_K]) 
    
    # integrating the equations
    
    r_bound = (R_kpc, 100*R_kpc)
    
    r_eval = np.linspace(r_bound[0],r_bound[1],resol)

    if print_time:
      t_ivp = time.perf_counter()

    sol = si.solve_ivp(diff_system_fast, r_bound, y0, t_eval=r_eval,\
                       args=( SFR_pure, redshift, M_vir_pure, f_esc_ion, f_esc_FUV, Plw, Ph1, Pg1, Zeta, A_NFW, r_s),\
                       method = integrator, events=stopping_condition) #,rtol=1.0e-3
    
    if print_time:
      time_ivp = (time.perf_counter() - t_ivp)
      print("total time ivp (s)=", time_ivp)

    if sol.success == False:
        print('Integration stopped before reaching the end of the array')
    elif sol.success == True:
        print('Integration completed successfully')
    
    r_kpc = sol.t # in kpc
    v_kms = sol.y[0] # in km/s
    n_cm3 = sol.y[1] # in cm-3
    T_K = sol.y[2] # in K
    
    print(sol.t_events)
    
    print(sol.y_events)
    
    # getting back to the old dimensions
    
    r = r_kpc * 1e3 * nc.pc # in cm
    v = v_kms * 1e5    # in cm/s
    n = n_cm3 # in cm-3
    T = T_K # in K
    
    mask = v > 0. 
    
    profiles = []
    
    profiles.append(v[mask])
    profiles.append(n[mask])
    profiles.append(T[mask])
    
    return sol_profiles(radius=r[mask], variables=profiles, params=params)



params = dict([("DM_model", "NFW"),
                   ("beta", 3.0), 
                   ("SFR", 50.),
                   ("f_esc_ion", 0.0), 
                   ("f_esc_FUV", 0.0), 
                   ("M_vir", 5e11),
                   ("redshift", 5.0),
                   ("Zeta", 1.0),
                   ("alfa", 1.0),
                   ("R_in", 0.3)])

time_start = time.perf_counter()

print("##################################")
print("Run with the following parameters:")
print("##################################")
print(params)
print("##################################")


import matplotlib.pyplot as plt


fig_sol, axs_sol = pltc.plot_configurator(plot_type="sol")    


integrator_list = ["RK45"]
#integrator_list = ["BDF","LSODA"]
#integrator_list = ["LSODA"]

show_profile    = True

folder = "plot_fast_solver"

if not os.path.exists(os.path.join(mydir.plot_dir, folder)):
    os.mkdir(os.path.join(mydir.plot_dir, folder))

for integrator in integrator_list:

    print(integrator)
    time_profile = time.perf_counter()
    profiles_new = get_profiles_fast(params, resol=1000,print_time=True,integrator=integrator)
    time_profile = (time.perf_counter() - time_profile)

    print("total profile time new (s)=", time_profile)
    
    
#    time_profile = time.perf_counter()
#    profiles_old = get_profiles(params, resol=1000,print_time=True,integrator=integrator)
#    time_profile = (time.perf_counter() - time_profile)
#
#    print("total profile time old (s)=", time_profile)

if show_profile:
    profiles_new.plot(ax=axs_sol, label=integrator)
#    profiles_old.plot(ax=axs_sol)

if show_profile:
    fig_sol.legend(loc="lower center", ncol=8, fontsize="small")
    plt.savefig(os.path.join(mydir.plot_dir, folder, "profiles.png"))

    
    



