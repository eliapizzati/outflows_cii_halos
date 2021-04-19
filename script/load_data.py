# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 10:42:59 2021

@author: anna
"""

import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table


import mydir
import natconst as nc 

from my_utils import get_circular_velocity_profile_NFW, get_virial_radius
from data_classes import obs_data

from astropy.cosmology import Planck15 as cosmo


"""
ALPINE
"""

#define the input file


names = ["DEIMOS_COSMOS_351640", "DEIMOS_COSMOS_396844", "DEIMOS_COSMOS_416105", "DEIMOS_COSMOS_488399",\
         "DEIMOS_COSMOS_494057", "DEIMOS_COSMOS_494763", "DEIMOS_COSMOS_539609", "DEIMOS_COSMOS_630594",\
         "DEIMOS_COSMOS_683613", "DEIMOS_COSMOS_709575", "DEIMOS_COSMOS_733857", "DEIMOS_COSMOS_834764",#"DEIMOS_COSMOS_848185",\
         "DEIMOS_COSMOS_880016", "DEIMOS_COSMOS_881725", #"vuds_cosmos_510596653",\
         "vuds_cosmos_5100537582", "vuds_cosmos_5100969402", "vuds_cosmos_5100994794",\
         #"vuds_cosmos_5101218326", 
         "vuds_cosmos_5110377875", "vuds_efdcs_530029038"]


names_CII_halo = ["DEIMOS_COSMOS_396844", "DEIMOS_COSMOS_488399", "DEIMOS_COSMOS_630594", "DEIMOS_COSMOS_683613",\
                  "DEIMOS_COSMOS_880016", "DEIMOS_COSMOS_881725", \
                  "vuds_cosmos_5100537582", "vuds_cosmos_5110377875"]

names_wo_CII_halo = ["DEIMOS_COSMOS_351640", "DEIMOS_COSMOS_416105", "DEIMOS_COSMOS_539609", \
                     "DEIMOS_COSMOS_709575", "DEIMOS_COSMOS_733857"]#"vuds_cosmos_510596653",

names_other = ["DEIMOS_COSMOS_494057", "DEIMOS_COSMOS_494763", "DEIMOS_COSMOS_834764", \
               "vuds_cosmos_5100969402", "vuds_cosmos_5100994794", "vuds_efdcs_530029038"]


names_short = ["DC_351640", "DC_396844", "DC_416105", "DC_488399",\
         "DC_494057", "DC_494763", "DC_539609", "DC_630594",\
         "DC_683613", "DC_709575", "DC_733857", "DC_834764",#"DC_848185",\
         "DC_880016", "DC_881725", #"vc_510596653",\
         "vc_5100537582", "vc_5100969402", "vc_5100994794",\
         #"vc_5101218326", 
         "vc_5110377875", "ve_530029038"]


names_CII_halo_short = ["DC_396844", "DC_488399", "DC_630594", "DC_683613", "DC_880016",\
                        "DC_881725", "vc_5100537582", "vc_5110377875"]

names_wo_CII_halo_short = ["DC_351640", "DC_416105", "DC_539609", \
                     "DC_709575", "DC_733857"]#"vc_510596653",

names_other_short = ["DC_494057", "DC_494763", "DC_834764", "vc_5100969402", "vc_5100994794", "ve_530029038" ]


obs_data_list = []

redshift=6

size = 14


if __name__ == "__main__":
    fig_tot, ax_tot = plt.subplots(1, 1, sharex=True, figsize=(1.5*8.27,1.5*4.))
    ax_tot.set_xlabel("b [kpc]", size=size)
    ax_tot.set_ylabel("data [mJy/arcsec^2]", size=size)
    ax_tot.tick_params(labelsize=size)
    #ax_tot.set_xlim(0.,15)
    ax_tot.set_ylim(1e-7,1e1)
    ax_tot.set_yscale("log")
    
    plt_star, ax_star = plt.subplots()




input_filename_general = os.path.join(mydir.script_dir, "input_data", "ALPINE_merged_catalogs.fits")

input_filename_behroozi =  os.path.join(mydir.script_dir, "input_data", "behroozi_z_5.dat")

#First, you can open the file and check the information in it
fits_table_hdu = fits.open(input_filename_general)

table_header = fits_table_hdu[1].header
table_data = fits_table_hdu[1].data

evt_data = Table(table_data)

halo_masses, halo_mass_ratio = np.loadtxt(input_filename_behroozi, unpack=True, usecols=(0,1))  

stellar_masses = halo_masses + halo_mass_ratio


    
if __name__ == "__main__":
    print("######################################################")
    print("name", "\t", "SFR","\t", "Mstar(1e10)","\t", "Mhalo(1e11)","\t", "redshift", "\t", "v_c(1e5)")
    print("######################################################")
    

for name, name_short in zip(names, names_short):
    
    
    L_CII = evt_data[evt_data["name"] == name]["LCII"][0]
    redshift = evt_data[evt_data["name"] == name]["zCII"][0]
    CII_FWHM = evt_data[evt_data["name"] == name]["CII_FWHM"][0]
    log_M_star = evt_data[evt_data["name"] == name]["logMstar"][0]
    log_SFR = evt_data[evt_data["name"] == name]["logSFR_SED"][0]
    
    
    index_masses = np.searchsorted(stellar_masses, log_M_star)
    
    log_M_halo = halo_masses[index_masses]
    
    M_halo = 10**log_M_halo # in solar masses
        
    R_vir = get_virial_radius(M_halo, redshift)
    
    #v_c = 23.4 * (M_halo/1e8)**(1./3) * ((1.+redshift)/10)**(1./2) * cosmo.h**(1./3) * 1e5
    v_c = np.sqrt(nc.gg*M_halo*nc.ms/R_vir)
    
    
    
    input_filename_CII = os.path.join(mydir.script_dir, "input_data", "{}.dat".format(name))
    input_filename_psf = os.path.join(mydir.script_dir, "input_data", "{}_psf.dat".format(name))
    
    data = np.loadtxt(input_filename_CII, unpack=True)  
    
    fuji_x_arcsec, fuji_y, err_down, err_up = data
    
    fuji_x = fuji_x_arcsec / cosmo.arcsec_per_kpc_proper(redshift).value * 1e3*nc.pc  # in cm
    
    fuji_y *= 1e3 # from Jy to mJy
    err_up *= 1e3 # from Jy to mJy
    err_down *= 1e3 # from Jy to mJy
    
    
    data_beam = np.loadtxt(input_filename_psf, unpack=True)  
    
    beam_x_arcsec, beam_y = data_beam
    
    beam_x = beam_x_arcsec / cosmo.arcsec_per_kpc_proper(redshift).value * 1e3*nc.pc  # in cm
    
    #beam_y *= 1e3 # from Jy to mJy
    
    #beam_y /= np.trapz(2*np.pi*beam_y, beam_x)
    
    
    if name in names_CII_halo:
        halo_class = "CII_halo"
    elif name in names_wo_CII_halo:
        halo_class = "wo_CII_halo"
    elif name in names_other:
        halo_class = "other"
    else:
        print(name)
        raise ValueError("No correct halo class")
        
    params_obs = dict([("name", name),
                       ("name_short", name_short),
                       ("halo_class", halo_class),
                       ("redshift", redshift),
                       ("line_FWHM", CII_FWHM*nc.km),
                       ("M_star", 10**log_M_star),
                       ("M_vir", M_halo),
                       ("SFR", 10**log_SFR),
                       ("v_c", v_c/1e5),
                       ("sersic_effective_radius", 1.1),
                       ("sersic_index", 1.),
                       ("beta_best_fit", None),
                       ("likelihood_best_fit", 0.)])
                    
    
    err = [err_down,err_up]
    
    observational_data = obs_data(x_data=fuji_x, y_data=fuji_y, x_beam=beam_x, y_beam=beam_y, err=err, params_obs=params_obs)
    
    
    obs_data_list.append(observational_data)
    
    if __name__ == "__main__":
        #observational_data.plot()
        observational_data.plot(ax=ax_tot)
        #observational_data.print_values()

        ax_star.scatter(log_M_star, log_SFR)
        print(name_short, "\t", round(10**log_SFR),"\t", round(10**log_M_star/1e10, 2),"\t", round(M_halo/1e11,2),\
              "\t", round(redshift,2), "\t", round(v_c/1e5))
    

#
#

"""
OLD --> FUJIMOTO+19
"""

input_filename = os.path.join(mydir.script_dir, "input_data", "radial_cii_ALMA-ALL.dat")

data = np.loadtxt(input_filename, unpack=True)  

fuji_x_arcsec, fuji_y, err_down, err_up = data

fuji_x = fuji_x_arcsec / 0.17076006194998467 * 1e3*nc.pc  # in cm

fuji_y *= 1e3 # from Jy to mJy
err_up *= 1e3 # from Jy to mJy
err_down *= 1e3 # from Jy to mJy


input_filename_beam = os.path.join(mydir.script_dir, "input_data", "radial_psf_ALMA-ALL.dat")

data_beam = np.loadtxt(input_filename_beam, unpack=True)  

beam_x_arcsec, beam_y, beam_err_down, beam_err_up = data_beam

beam_x = beam_x_arcsec / 0.17076006194998467 * 1e3*nc.pc  # in cm

beam_y *= 1e3 # from Jy to mJy
beam_err_up *= 1e3 # from Jy to mJy
beam_err_down *= 1e3 # from Jy to mJy

beam_y /= np.trapz(2*np.pi*beam_y, beam_x)


#
params_obs = dict([("line_frequency", 1900*1e9),
                   ("redshift", 5.96),
                   ("line_FWHM", 296*nc.km),
                   ("sersic_effective_radius", 1.1),
                   ("sersic_index", 1.2),
                   ("exp_effective_radius", 3.3)])

err = [err_down,err_up]


observational_data_fuji = obs_data(x_data=fuji_x, y_data=fuji_y, x_beam=beam_x, y_beam=beam_y, err=err, params_obs=params_obs)

