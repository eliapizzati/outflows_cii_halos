# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 19:49:55 2021

@author: anna
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 00:49:01 2019

@author: anna
"""


from functools import partial
import scipy.constants as sc
import scipy.integrate as si
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
import natconst as nc



def f_in(r,M):
    return ((6./(5.+1/M**2))**(9/14))*((1.+3./M**2)/4.)**(1/7)-r

def inn(r):
    a=0.00001
    b=1.
    return so.brentq(partial(f_in,r),a,b)


def f_outer(r,M):
    return M**3*((1.+3./M**2)/4.)**2-r**2

def outer(r):
    a=1.0
    b=1000
    return so.brentq(partial(f_outer,r),a,b)



def cc85_profiles(beta, SFR):

    #set galactic variables

    #set physical variables

    r_val=np.linspace(0.0001,1,10000)
    ys = []
    xs = []

    for r in r_val:
        try:
            y = inn(r)
        except ValueError:
        # Should we not be able to find a solution in this window.
            assert 'problema'
        else:
            ys.append(y)
            xs.append(r)

    M_in = np.asarray(ys)
    r_in = np.asarray(xs)

        
    r_val=np.linspace(1.0,40.,10000)
    ys = []
    xs = []

    for r in r_val:
        try:
            y = outer(r)
        except ValueError:
        # Should we not be able to find a solution in this window.
            assert 'problema'
        else:
            ys.append(y)
            xs.append(r)
            
    M_out = np.asarray(ys)
    r_out = np.asarray(xs)

    r_tot = np.concatenate((r_in,r_out))
    M_tot = np.concatenate((M_in,M_out))


    alfa=1. #energy loading factor

    E_SN = 1e51 * SFR /100. /nc.year #erg (energy from supernovae-- energy of a single SN \
                        #per 100 M_sun of Star Formation)


    R = 0.3*1000*nc.pc #cm

    mu=0.61 #mean molecular weight of the sun

    k=nc.kk/(mu*nc.mp) 


    V=4*np.pi*R**3/3

    M_dot = beta*SFR*nc.ms/nc.year  #mass from SN
    
    E_dot = alfa*E_SN #erg/s
  
    Q = E_dot/V

    q = M_dot/V
 
    
    #inner zone
    
    c_in=np.sqrt((Q/q)/(M_in**2/2+3/2))
    
    v_in=M_in*c_in
    
    p_in=(q*r_in*R)/(3*v_in)
    
    P_in=p_in*c_in**2*3/5
    
    T_in=P_in/(p_in*k)
    
    M0=1.
    
    v0=v_in[len(r_in)-1]
    
    c0=v0/M0
    
    p0=p_in[len(r_in)-1]
    
    P0=P_in[len(r_in)-1]
    
    T0=P0/(p0*k)  #K
    
    
    c_out=np.sqrt(((M0**2/2+3/2)*c0**2)/(M_out**2/2+3/2))
    
    v_out=M_out*c_out
    
    p_out=(v0*p0)/(v_out*(r_out**2))
    
    
    P_out=p_out**(5/3)*(P0/(p0**(5/3)))
    
    T_out=P_out/(p_out*k)  
   
    
    #total
    
    v=np.concatenate((v_in,v_out))
    T=np.concatenate((T_in,T_out))
    p=np.concatenate((p_in,p_out))

    return r_tot, v, p, T
