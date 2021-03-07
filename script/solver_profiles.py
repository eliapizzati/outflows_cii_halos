"""
Created on Fri Mar  5 16:12:54 2021

@author: anna
"""

import os
import natconst as nc
import mydir 
import gnedincooling as gc
import radiation_fields

import scipy.integrate as si

import numpy as np


# cooling function


gc.frtinitcf(0, os.path.join(mydir.script_dir, "input_data", "cf_table.I2.dat"))




def lamda(T, n, r, SFR_pure, Zeta=1., f_esc=0.):
    
    Pc6 = 0.

    # UV background (values around 1e-13)
        
        
    Plw = nc.Gamma_LW_UVB
    Ph1 = nc.Gamma_H_UVB
    Pg1 = nc.Gamma_He_UVB         
        
    # ADDING THE GALACTIC FLUX CONTRIBUTION
        
    Plw += nc.Gamma_LW_1000 * (1000*nc.pc/r)**2  * SFR_pure*f_esc
    Ph1 += nc.Gamma_H_1000 * (1000*nc.pc/r)**2 * SFR_pure*f_esc
    Pg1 += nc.Gamma_He_1000 * (1000*nc.pc/r)**2 * SFR_pure*f_esc

    return gc.frtgetcf_cool(T, n, Zeta, Plw, Ph1, Pg1, Pc6) - gc.frtgetcf_heat(T, n, Zeta, Plw, Ph1, Pg1, Pc6)



# SYSTEM OF EQUATIONS

def diff_system(r, y, v_c_pure, SFR_pure, f_esc=0., Zeta=1.):

    assert len(y.shape) == 1, 'Dimensionality of the input array is wrong'

    v = y[0] # velocity
    rho = y[1] # density
    T = y[2] # temperature
        
    c_S2 = nc.gamma*nc.knorm * T
    c_T2 = nc.knorm * T
    
    n = rho / (nc.mus * nc.mp) #in cgs unit
    
    q = rho*lamda(T, n, r, Zeta=Zeta, f_esc=f_esc, SFR_pure=SFR_pure) / (nc.mus * nc.mp)**2
    
    v_c = v_c_pure *1e5 #cm/s 
    v_e = v_c * np.sqrt(2) #cm/s 

    
    output_a = ((2*v/r)*(c_S2 - v_e**2/4) + (nc.gamma-1)*q) / (v**2-c_S2)
    
    output_b = ((2*rho/r)*(v_e**2/4-v**2) - (nc.gamma-1)*q*rho/v) / (v**2-c_S2)

    output_c = ((nc.gamma-1)*(2*T/r)*(v_e**2/4-v**2) - (nc.gamma-1)*q*(v**2-c_T2)/(nc.knorm*v)) / (v**2-c_S2)        
    
    output = np.asarray([output_a,output_b,output_c])
    
    assert len(output.shape) == 1  , 'Dimensionality of the output array is wrong'
    
    return output



def solve_eq(beta, SFR_pure, v_c_pure, alfa = 1., R_pure = 0.3, resol = 1000):
    
    # getting the BC
    
    SFR = SFR_pure/nc.year #1/s
    
    E_SN = 1e51*(SFR)/100 #erg (energy from supernovae-- energy of a single SN \
                        #per 100 M_sun of Star Formation)

    M_dot = beta*SFR*nc.ms  #mass from SN
    
    E_dot = alfa*E_SN #erg/s
  
    #M0 = 1. 
    v0 = np.sqrt(E_dot/M_dot)/np.sqrt(2)  #m/s
    #c0 = v0/M0

    R = R_pure*1000*nc.pc #cm

    rho0 = 0.1125395*np.sqrt(M_dot**3/E_dot) / R**2 #g/cm^3
    P0 = 0.0337618*np.sqrt(M_dot*E_dot) / R**2  #g/cm^2    
    T0 = P0/(rho0*nc.knorm)  #K
    
    y0 = np.asarray([v0,rho0,T0]) 
    
    # integrating the equations
    
    r_bound = (R, 100*R)
    
    r_eval = np.linspace(r_bound[0],r_bound[1],resol)

    sol = si.solve_ivp(diff_system, r_bound, y0, t_eval=r_eval)
    
    if sol.success == False:
        print('Error in integration procedure')
        
    
    r = sol.t #cm
    v = sol.y[0] #cm/s
    rho = sol.y[1]
    n = rho / (nc.mus*nc.mp) #cm^-3
    T = sol.y[2] #K
    
    profiles = []
    
    profiles.append(v)
    profiles.append(n)
    profiles.append(T)
    
    return r, profiles
    

    
class profile():
    
    def __init__(self, radius, variable, params):
        
        self.r = radius
        self.var = variable
        self.params = params

    
    def to_file(self):
        
        folder = 'data_profiles_SFR{}_vc{:.1f}'.format(self.params["SFR"], self.params["v_c"])
            
        np.savetxt(os.path.join(mydir.data_dir, folder, "{}_beta{:.2f}.dat".format(self.params["type"], self.params["beta"])), \
                   (self.r,self.var))
        
        return
    
        
    def plot(self, ax):
        
        ax.plot(self.r, self.var)
        return
    
    def check_nans(self):
        
        return np.isnan(np.sum(self.var))
    
    

        