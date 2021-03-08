# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 14:50:57 2021

@author: anna
"""


import numpy as np

import scipy.interpolate


def sersic(r):
    
    re= 1.1
    sersic_index = 1.2
    b=2*sersic_index-1/3
    I=7.8*np.exp(b*((1/re)**(1/sersic_index)-1))
    
    return I*np.exp(-b*((r/re)**(1/sersic_index)-1))

def halo(r):
    
    rn=3.3
    C=1.5*np.exp(1/rn)
    
    return C*np.exp(-r/rn)


def twod_making(profile, x_axis, nimage=1000):

    x_ext = np.linspace(-x_axis.max(), x_axis.max(), nimage)
    
    y_ext = x_ext
    
    xy, yx = np.meshgrid(x_ext,y_ext)
    
    x_1d = np.linspace(-x_axis.max()*np.sqrt(2), x_axis.max()*np.sqrt(2), nimage)
    
    profile_ext_pos = np.interp(x_1d[x_1d>0.], x_axis, profile, right=0.)
    
    profile_ext_neg = profile_ext_pos[::-1]
    
    profile_ext = np.concatenate((profile_ext_neg, profile_ext_pos))
        
    function_ext = scipy.interpolate.interp1d(x_1d, profile_ext)
    
    z = function_ext(np.sqrt(xy**2 + yx**2))
    
    return x_ext, y_ext, z


def create_beam_from_data(r_beam_raw, beam_raw):
    
    h_data = np.linspace(0.3,10,nimage)

    h_data_ext = np.linspace(0,10,nimage)

    beam = np.interp(h_data, beam_x, beam_y)

    #cut_low = 0.
    #mask    = beam < cut_low
    #beam[mask] =  cut_low


    beam_ext = np.interp(h_data_ext, h_data, beam)
    normalization_beam = np.sqrt(np.trapz(2*np.pi * h_data_ext * beam_ext**2, h_data_ext))

    beam /= normalization_beam

    return beam


