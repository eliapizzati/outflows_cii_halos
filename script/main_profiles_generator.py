# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 16:39:49 2021

@author: anna
"""



import os
import numpy as np
import logging


from load_data import obs_data_list, names, names_CII_halo, names_wo_CII_halo, names_other,  observational_data_fuji
from sol_modules import get_profiles
from my_utils import get_vc_from_virial_mass
import mydir

import time


verbose = True
err_switcher = True


filename_log = str(input("filename for the logger:"))

logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', \
                    datefmt='%m/%d/%Y %I:%M:%S %p',\
                    level=logging.DEBUG,\
                    handlers=[logging.FileHandler(os.path.join(mydir.log_dir, "{}.log".format(filename_log))),\
                              logging.StreamHandler()])


redshifts = []
SFRs = [] 
M_virs = []
datas = []

halo_class = input("which halo class (CII_halo, wo_CII_halo, other)?")

if err_switcher: 
    err = input("which error (none, up, down)?")
    
for data in obs_data_list:
    if data.params_obs["halo_class"] != halo_class: #or data.params_obs["name"] != "DEIMOS_COSMOS_881725":
    #if data.params_obs["name"] in names_wo_CII_halo or data.params_obs["name"] in names_CII_halo:#names_wo_CII_halodata.params_obs["name"] != "DEIMOS_COSMOS_881725":
    #if data.params_obs["name"] != "vuds_cosmos_5110377875":
        pass

    else:
        redshifts.append(data.params_obs["redshift"])
        SFRs.append(data.params_obs["SFR"])
        datas.append(data)
        
        if err_switcher:
            if err == "none":
                M_virs.append(data.params_obs["M_vir"])
            elif err == "up":
                M_virs.append(data.params_obs["M_vir"]+data.params_obs["M_vir_err_up"])
            elif err == "down":
                M_virs.append(data.params_obs["M_vir"]-data.params_obs["M_vir_err_down"])
        else:
            M_virs.append(data.params_obs["M_vir"])
        
        
f_esc_ion = 0.0

f_esc_FUVs = [0.0]#,0.05,0.1,0.2]

#○betas = np.linspace(1.0,4.0, 31)
betas = np.linspace(1.0,4.0, 31)

#v_cs = np.linspace(175.,300., 26)

params = dict([("DM_model", "NFW"),
           ("beta", 1.0), 
           ("SFR", 50.),
           ("f_esc_ion", f_esc_ion), 
           ("f_esc_FUV", 0.0), 
           ("v_c", 175.),
           ("M_vir", 1e11),
           ("redshift", 5.),
           ("Zeta", 1.0),
           ("alfa", 1.0),
           ("R_in", 0.3)])    

    
for SFR, M_vir, z, counter_data in zip(SFRs, M_virs, redshifts, range(len(datas))):
    
    if verbose:
          logging.info("#################################################################################")
          logging.info("updating  obs parameters (v_c, M_vir, SFR, z), number {}/{}".format(counter_data+1, len(datas)))
          logging.info("#################################################################################")
    
    for f_esc_FUV, counter_fesc in zip(f_esc_FUVs, range(len(f_esc_FUVs))):
    
        if verbose:
                logging.info("#################################################################################")
                logging.info("updating  f_esc_FUV, number {}/{}".format(counter_fesc+1, len(f_esc_FUVs)))
                logging.info("#################################################################################")
           
        for beta in betas:
            
            time_start = time.perf_counter()

            v_c = get_vc_from_virial_mass(M_vir, z)/1e5
            
            params.update(beta = beta)
            params.update(f_esc_FUV = f_esc_FUV)
            params.update(v_c = v_c)
            params.update(M_vir = M_vir)
            params.update(SFR = SFR)
            params.update(redshift = z)
            
            if verbose:
                logging.info("run with beta = {:.1f}, f_esc_FUV = {:.2f} ({:d}/{:d}), M_vir = {:.2e}, v_c = {:.1f}, SFR = {:.1f}, z = {:.1f} ({:d}/{:d})"\
                      .format(beta, f_esc_FUV, counter_fesc+1, len(f_esc_FUVs), \
                              M_vir, v_c, SFR, z, counter_data+1, len(datas)))
            
            profiles = get_profiles(params, resol=1000)
        
            profiles.to_file()
                
            if profiles.check_nans() == True:
                string_nans = "Integration error: presence of nans"
            else:
                string_nans = "Integration successful"
                    
            logging.info("beta= %.2f \t %s", beta, string_nans)
            
            time_elapsed = (time.perf_counter() - time_start)
        
            logging.info("total time elapsed (s) = %.1f", time_elapsed)




