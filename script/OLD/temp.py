# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 22:35:50 2021

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



#set galactic variables

#set physical variables

alfa=1. #energy loading factor

beta=0.2 #mass loading factor

SFR=10./nc.year #1/s

E_SN=1e51*(SFR)/100 #erg (energy from supernovae-- energy of a single SN \
                        #per 100 M_sun of Star Formation)

M_dot=beta*SFR*nc.ms #mass from SN

E_dot=alfa*E_SN #erg/s

R=0.3*1000*nc.pc #cm

gamma=5/3

v_e=400*1e5 #cm/s

mu=0.61 #mean molecular weight of the sun

k=nc.kk/(mu*nc.mp) 


V=4*np.pi*R**3/3

Q=E_dot/V

q=M_dot/V

#inner zone

a=0.00001
b=1.


def f_in(r,M):
    return ((6./(5.+1/M**2))**(9/14))*((1.+3./M**2)/4.)**(1/7)-r

def inn(r):
    return so.brentq(partial(f_in,r),a,b)

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
        

M_in=np.asarray(ys)
r_in=np.asarray(xs)

plt.figure()
plt.plot(r_in, M_in)
plt.grid(b=True, which='both', color='0.65', linestyle='-')
plt.xscale('log')
plt.xlabel("r/R", size=18)
plt.ylabel("M", size=18)

print(len(r_val), len(r_in))
  

#outer zone


a=1.0
b=1000


r_val=np.linspace(1.0,40.,10000)
ys = []
xs = []


def f_outer(r,M):
    return M**3*((1.+3./M**2)/4.)**2-r**2

def outer(r):
    return so.brentq(partial(f_outer,r),a,b)

for r in r_val:
    try:
        y = outer(r)
    except ValueError:
        # Should we not be able to find a solution in this window.
        assert 'problema'
    else:
        ys.append(y)
        xs.append(r)
        

M_out=np.asarray(ys)
r_out=np.asarray(xs)

plt.figure()
plt.plot(r_out, M_out)
plt.grid(b=True, which='both', color='0.65', linestyle='-')
plt.xscale('log')
plt.xlabel("r/R", size=18)
plt.ylabel("M", size=18)


print(len(r_val), len(r_out))


#total



r_tot=np.concatenate((r_in,r_out))
M_tot=np.concatenate((M_in,M_out))



plt.figure()
plt.xscale('log')
plt.grid(b=True, which='both', color='0.65', linestyle='-')
plt.plot(r_tot, M_tot)
plt.xlim((1e-2,np.max(r_tot)))
plt.xlabel("r/R", size=18)
plt.ylabel("M", size=18)



#physical variables

'''
I PUT THE LEGEND AND EVERYTHING JUST FOR THE TEMPERATURE ONE, SO DO IT ALSO WITH THE OTHERS IF YOU NEED TO'''


while beta<=3.5:
    
    M_dot=beta*SFR*nc.ms  #mass from SN
    
    E_dot=alfa*E_SN #J/s
  
    Q=E_dot/V

    q=M_dot/V
 
    
    #inner zone
    
    
    c_in=np.sqrt((Q/q)/(M_in**2/2+3/2))
    
    v_in=M_in*c_in
    
    p_in=(q*r_in*R)/(3*v_in)
    
    
    P_in=p_in*c_in**2*3/5
    
    T_in=P_in/(p_in*k)
    
    '''
    #speed
    plt.figure()
    plt.grid(b=True, which='both', color='0.65', linestyle='-')
    plt.xscale('log')
    plt.plot(r_in*R/(1000*nc.pc), v_in/10**8, color='blue')
    plt.xlabel("r [kpc]", size=18)
    plt.ylabel("v [1000 km/s] ", size=18)
    
    #pressure
    plt.figure()
    plt.grid(b=True, which='both', color='0.65', linestyle='-')
    plt.xscale('log')
    plt.plot(r_in*R/(1000*nc.pc), P_in, color='red')
    plt.xlabel("r [kpc]", size=18)
    plt.ylabel("P [g/cm^2] ", size=18)
    
    #temperature
    plt.figure()
    plt.grid(b=True, which='both', color='0.65', linestyle='-')
    plt.xscale('log')
    plt.plot(r_in*R/(1000*nc.pc), T_in, color='gray')
    plt.xlabel("r [kpc]", size=18)
    plt.ylabel("T [K] ", size=18)
    
    #density
    plt.figure()
    plt.grid(b=True, which='both', color='0.65', linestyle='-')
    plt.xscale('log')
    plt.plot(r_in*R/(1000*nc.pc), p_in/(nc.mp*mu), color='orange')
    plt.xlabel("r [kpc]", size=18)
    plt.ylabel("n [cm^(-3)] ", size=18)
    
    '''
    #outer zone
    
    '''
    #value from the paper
    M0=1. 
    
    v0=np.sqrt(E_dot/M_dot)/np.sqrt(2)  #m/s
    
    c0=v0/M0
    
    p0=0.113*np.sqrt(M_dot**3/E_dot)/R**2 #kg/m^3
    
    P0=0.0338*np.sqrt(M_dot*E_dot)/R**2  #Pa
    
    T0=P0/(p0*k)  #K
    '''
    #value from integration
    
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
    
    '''
    #speed
    plt.figure()
    plt.grid(b=True, which='both', color='0.65', linestyle='-')
    plt.xscale('log')
    plt.plot(r_out*R/(1000*nc.pc), v_out/10**8, color='blue')
    plt.xlabel("r [kpc]", size=18)
    plt.ylabel("v [1000 km/s] ", size=18)
    
    #pressure
    plt.figure()
    plt.xscale('log')
    plt.grid(b=True, which='both', color='0.65', linestyle='-')
    plt.plot(r_out*R/(1000*nc.pc), P_out, color='red')
    plt.xlabel("r [kpc]", size=18)
    plt.ylabel("P [g/cm^2] ", size=18)
    
    #temperature
    plt.figure()
    plt.xscale('log')
    plt.grid(b=True, which='both', color='0.65', linestyle='-')
    plt.plot(r_out*R/(1000*nc.pc), T_out, color='gray')
    plt.xlabel("r [kpc]", size=18)
    plt.ylabel("T [K] ", size=18)
    
    #density
    plt.figure()
    plt.grid(b=True, which='both', color='0.65', linestyle='-')
    plt.xscale('log')
    plt.plot(r_out*R/(1000*nc.pc), p_out/(mu*nc.mp), color='orange')
    plt.xlabel("r [kpc]", size=18)
    plt.ylabel("n [cm^(-3)] ", size=18)
    '''
    #total
    
    v=np.concatenate((v_in,v_out))
    P=np.concatenate((P_in,P_out))
    T=np.concatenate((T_in,T_out))
    p=np.concatenate((p_in,p_out))
    '''
    #speed
    plt.figure()
    plt.grid(b=True, which='both', color='0.65', linestyle='-')
    plt.xscale('log')
    plt.plot(r_tot*R/(1000*nc.pc), v/10**8, color='blue')
    plt.xlim((1e-2,np.max(r_tot*R/(1000*nc.pc))))
    plt.xlabel("r [kpc]", size=18)
    plt.ylabel("v [1000 km/s] ", size=18)
    
    #pressure
    plt.figure()
    plt.xscale('log')
    plt.grid(b=True, which='both', color='0.65', linestyle='-')
    plt.plot(r_tot*R/(1000*nc.pc), P, color='red')
    plt.xlim((1e-2,np.max(r_tot*R/(1000*nc.pc))))
    plt.xlabel("r [kpc]", size=18)
    plt.ylabel("P [g/cm^2] ", size=18)
    '''
    #temperature
    y_axis=np.linspace(4.5,8.0,1000)
    plt.figure(0)
    plt.ylim((4.5,8))
    #plt.xscale('log')
    #plt.grid()
    plt.plot(np.log10(r_tot), np.log10(T), label='$\\beta$={0:.1f}'.format(beta))
    #plt.plot(np.asarray([0]*len(y_axis)), y_axis, color='gray',alpha=0.05)
    plt.xlim((-1.0,np.log10(np.max(r_tot))))
    plt.xlabel("log (r / R)", size=18)
    plt.ylabel("log (T [K]) ", size =18)
    plt.tick_params(labelsize=18)
    plt.tight_layout()
    plt.legend(loc='best')

    '''
    #density
    plt.figure()
    plt.grid(b=True, which='both', color='0.65', linestyle='-')
    plt.plot(r_tot*R/(1000*nc.pc), p/(sc.m_p*10**8), color='orange')
    plt.xlim((1e-2,np.max(r_tot*R/(1000*nc.pc))))
    plt.xscale('log')
    plt.xlabel("r [kpc]", size=18)
    plt.ylabel("n [cm^(-3)] ", size=18)
    
    #check of integral
    plt.figure()
    plt.grid(b=True, which='both', color='0.65', linestyle='-')
    plt.xscale('log')
    plt.plot(r_tot*R/(1000*nc.pc), v**2/2+(5/2)*P/p)
    plt.plot(r_tot*R/(1000*nc.pc), np.asarray([Q/q]*len(r_tot)))
    plt.xlabel("r [kpc]", size=18)
    plt.ylabel("B", size=18)
    '''
    #np.savetxt(r'C:\Users\anna\Google Drive (elia.pizzati@sns.it)\Università\Tesi\Script\Data\Data_CC85\beta{:.2f}.dat'.format(beta), (r_tot, T))


    beta=beta+0.9

