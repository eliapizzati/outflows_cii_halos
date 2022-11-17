# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 17:42:04 2021

@author: anna
"""


import numpy as np
import os

import natconst as nc
import mydir
from script.OLD.load_data import obs_data_list

from my_utils import get_concentration, get_virial_radius, get_circular_velocity_profile_NFW

input_filename_SFR =  os.path.join(mydir.script_dir, "input_data", "sfr_release.dat")

input_filename_behroozi =  os.path.join(mydir.script_dir, "input_data", "behroozi_z_5.dat")


halo_masses, halo_mass_ratio = np.loadtxt(input_filename_behroozi, unpack=True, usecols=(0,1))  
stellar_masses = halo_masses + halo_mass_ratio

    
if __name__ == "__main__":
    print("######################################################")
    print("name", "\t", "SFR","\t", "Mstar(1e10)","\t", "redshift","\t",\
          "Mhalo(1e11)","\t", "R_vir(kpc)","\t",\
          "concentration","\t", "v_c(1e5)")
    print("######################################################")



r = np.linspace(0.1*nc.pc*1e3, 50*nc.pc*1e3, 1000)

import matplotlib.pyplot as plt
from matplotlib import cm, colors

fig_rel, ax_rel = plt.subplots(figsize=(8.27,5.))

fig_abs, ax_abs = plt.subplots(figsize=(8.27,5.))

ax_rel.set_xscale("log")
ax_rel.set_xlabel("r/r_vir")
ax_rel.set_ylabel("v_c(r)/v_c,max")

ax_abs.set_xlabel("r [kpc]")
ax_abs.set_ylabel("v_c(r) [km/s]")
ax_abs.set_xlim(0.1,20)

cmap_rend_col = cm.get_cmap('viridis')
    
M_vir_min = 1
M_vir_max = 10
    
norm = colors.Normalize(vmin=M_vir_min, vmax=M_vir_max)
    
for data in obs_data_list:
    
    name_short = data.params_obs["name_short"]
    SFR = data.params_obs["SFR"]
    M_star = data.params_obs["M_star"]
    redshift = data.params_obs["redshift"]
    
    index_masses = np.searchsorted(stellar_masses, np.log10(M_star))

    log_M_halo = halo_masses[index_masses]

    M_halo = 10**log_M_halo # in solar masses
    
    R_vir = get_virial_radius(M_halo, redshift)
    
    concentration = get_concentration(M_halo, redshift)

    #v_c = 23.4 * (M_halo/1e8)**(1./3) * ((1.+redshift)/10)**(1./2) * cosmo.h**(1./3) * 1e5
    v_c = np.sqrt(nc.gg*M_halo*nc.ms/R_vir)
    
    
    print(name_short, "\t", round(SFR),"\t", round(M_star/1e10, 2),"\t", \
          round(redshift,2), "\t", round(M_halo/1e11,2), "\t",\
          round(R_vir/1e3/nc.pc,1),"\t", round(concentration,2), "\t", round(v_c/1e5))
    
    v_r = get_circular_velocity_profile_NFW(r, M_halo, redshift)
    
    ax_rel.plot(r/R_vir, v_r/v_c)
    ax_abs.plot(r/1e3/nc.pc, v_r/1e5, color = cmap_rend_col((M_halo/1e11-M_vir_min)/(M_vir_max-M_vir_min)))
    
     
ax_rel.axhline(1., linestyle='--', color='gray')


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
    
cb = fig_abs.colorbar(cmap, orientation='vertical')
cb.set_label(r'M_vir [1e11 msun]', rotation=90.)

fig_abs.tight_layout()
fig_rel.tight_layout()

plt.show()
 
    


