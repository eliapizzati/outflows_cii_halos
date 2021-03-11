# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 15:01:16 2021

@author: anna
"""

import numpy as np

class obs_data():
    
    def __init__(self, x_data, y_data, err, x_beam, y_beam, params_obs):
        
        self.x = x_data
        
        self.data = y_data
        
        self.err = err
        
        self.err_down =  err[0]
        self.err_up = err[1]
        
        self.x_beam = x_beam 
        
        self.beam = y_beam
        
        self.params_obs = params_obs
        
    
    
    def check_nans(self):
        
        var = np.asarray(self.data)
        
        return np.isnan(np.sum(var))
    



