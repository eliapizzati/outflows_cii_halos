# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 14:50:57 2021

@author: anna
"""


import numpy as np

import scipy.interpolate
import natconst as nc


def sersic(r, r_e, sersic_index, central):
        
    b = 2 * sersic_index - 1/3
    
    I_e = central * np.exp(b*((1/r_e)**(1/sersic_index)-1))
    
    return I_e*np.exp(-b*((r/r_e)**(1/sersic_index)-1.))

def halo(r, r_n, central):
        
    C = central * np.exp(1/r_n)
    
    return C*np.exp(-r/r_n)


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




