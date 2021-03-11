# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 10:42:59 2021

@author: anna
"""

import os
import numpy as np
from astropy.io import fits
from astropy.table import Table

import mydir
import natconst as nc 

from data_classes import obs_data


"""
ALPINE
"""
#
##define the input file
#input_filename = os.path.join(mydir.script_dir, "input_data", "ALPINE_merged_catalogs.fits")
#
#
##First, you can open the file and check the information in it
#fits_table_hdu = fits.open(input_filename)
#
##...like the number of HDU, in this case 5
#fits_table_hdu.info()
#
#
#table_header = fits_table_hdu[1].header
#table_data = fits_table_hdu[1].data
#
#
#evt_data = Table(table_data)
#print("Data contains %d events" % len(evt_data))


"""
FUJIMOTO+19
"""

input_filename = os.path.join(mydir.script_dir, "input_data", "radial_cii_ALMA-ALL.dat")

data = np.loadtxt(input_filename, unpack=True)  

fuji_x_arcsec, fuji_y, err_down, err_up = data

fuji_x = fuji_x_arcsec / 0.17076006194998467 * nc.pc  # in cm

fuji_y *= 1e3 # from Jy to mJy
err_up *= 1e3 # from Jy to mJy
err_down *= 1e3 # from Jy to mJy


input_filename_beam = os.path.join(mydir.script_dir, "input_data", "radial_psf_ALMA-ALL.dat")

data_beam = np.loadtxt(input_filename, unpack=True)  

beam_x_arcsec, beam_y, beam_err_down, beam_err_up = data

beam_x = beam_x_arcsec / 0.17076006194998467 * nc.pc  # in cm

beam_y *= 1e3 # from Jy to mJy
beam_err_up *= 1e3 # from Jy to mJy
beam_err_down *= 1e3 # from Jy to mJy

beam_y /= np.trapz(2*np.pi*beam_y, beam_x)


params_obs = dict([("line_frequency", 1900*1e9),
                   ("redshift", 5.96),
                   ("line_FWHM", 296*nc.km),
                   ("sersic_effective_radius", 1.1),
                   ("sersic_index", 1.2),
                   ("exp_effective_radius", 3.3)])

err = [err_down,err_up]


observational_data = obs_data(x_data=fuji_x, y_data=fuji_y, x_beam=beam_x, y_beam=beam_y, err=err, params_obs=params_obs)



