"""
This script is created to read the data stored with the script main_profile_saver.py
"""

import os
import numpy as np

import my_dir
import natconst as nc


from load_data import obs_data_list, observational_data_fuji
from model_classes import load_from_file, get_data_folder

import time


# Creating a dictionary for the parameters


params = dict([("DM_model", "NFW+disk"),
               ("beta", 6.8),
               ("SFR", 15.),
               ("f_esc_ion", 0.0),
               ("f_esc_FUV", 0.0),
               ("v_c", 198.),
               ("redshift", 5.0),
               ("Zeta", 1.0),
               ("alfa", 1.0),
               ("R_in", 0.3)])


# reading the profiles from file

folder = get_data_folder(params, "profiles")
filename = "profiles_beta6.80_SFR15.0_vc198.0_NFW+disk_VC5100537582_10000_per_laura_log_grid.dat"

profiles = load_from_file(params, filename=os.path.join(my_dir.data_dir, folder, filename),
                          class_type="profiles")

if profiles.check_nans() == True:
    string_nans = "Integration error: presence of nans"
else:
    string_nans = "Integration successful"
print(string_nans)

print(profiles.r[:100:]/nc.pc)
print(profiles.T[:100:])

