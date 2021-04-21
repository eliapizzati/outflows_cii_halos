# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 16:39:49 2021

@author: anna
"""



import os
import numpy as np


from load_data import obs_data_list, names, names_CII_halo, names_wo_CII_halo, names_other,  observational_data_fuji
from sol_modules import get_profiles
            

import time


verbose = True

v_cs = []
redshifts = []
SFRs = [] 

halo_class = input("which halo class (CII_halo, wo_CII_halo, other)?")

for data in obs_data_list:
    if data.params_obs["halo_class"] != halo_class: #or data.params_obs["name"] != "DEIMOS_COSMOS_881725":
    #if data.params_obs["name"] in names_wo_CII_halo or data.params_obs["name"] in names_CII_halo:#names_wo_CII_halodata.params_obs["name"] != "DEIMOS_COSMOS_881725":
    #if data.params_obs["name"] != "vuds_cosmos_5110377875":
        pass

    v_cs.append(data.params_obs["v_c"])
    redshifts.append(data.params_obs["redshift"])
    SFRs.append(data.params_obs["SFR"])
    
f_esc_ion = 0.0

f_esc_FUVs = [0.0]#,0.05,0.1,0.2]

betas = np.linspace(1.0,4.0, 31)

#v_cs = np.linspace(175.,300., 26)

params = dict([("DM_model", "NFW"),
           ("beta", 1.0), 
           ("SFR", 50.),
           ("f_esc_ion", f_esc_ion), 
           ("f_esc_FUV", 0.0), 
           ("v_c", 175.),
           ("redshift", 5.),
           ("Zeta", 1.0),
           ("alfa", 1.0),
           ("R_in", 0.3)])    

    
for v_c, SFR, z, counter_data in zip(v_cs, SFRs, redshifts, range(len(obs_data_list))):
    
    if verbose:
          print("#################################################################################")
          print("updating  obs parameters (v_c, SFR, z), number {}/{}".format(counter_data, len(obs_data_list)))
          print("#################################################################################")
    
    for f_esc_FUV, counter_fesc in zip(f_esc_FUVs, range(len(f_esc_FUVs))):
    
        if verbose:
                print("#################################################################################")
                print("updating  f_esc_FUV, number {}/{}".format(counter_fesc, len(f_esc_FUVs)))
                print("#################################################################################")
           
        for beta in betas:
            
            time_start = time.perf_counter()

            params.update(beta = beta)
            params.update(f_esc_FUV = f_esc_FUV)
            params.update(v_c = v_c)
            params.update(SFR = SFR)
            params.update(redshift = z)
            
            if verbose:
                print("run with beta = {:.1f}, f_esc_FUV = {:.2f}, v_c = {:.1f}, SFR = {:.1f}, z = {:.1f}".format(\
                      beta, f_esc_FUV, v_c, SFR, z))
            
            profiles = get_profiles(params, resol=1000)
        
            profiles.to_file()
                
            if profiles.check_nans() == True:
                string_nans = "Integration error: presence of nans"
            else:
                string_nans = "Integration successful"
                    
            print("beta=", beta, "\t", string_nans)
            
            time_elapsed = (time.perf_counter() - time_start)
        
            print("total time elapsed (s)=", time_elapsed)




