"""
This script defines the class for the observational data
"""

import numpy as np
import matplotlib.pyplot as plt
import natconst as nc

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
    


        
    def plot(self, ax=None, size=14):
        """
        Plot the data

        Parameters
        ----------
        ax: matplotlib axis
        size: int

        Returns
        -------
        None.
        """
        
        if ax is None:
            
            fig, ax = plt.subplots(1, 1, sharex=True, figsize=(1.5*8.27,1.5*4.))
    
            ax.set_xlabel("b [kpc]", size=size)
            
            ax.set_ylabel("data [mJy/arcsec^2]", size=size)
                          
            ax.tick_params(labelsize=size)
            #ax.set_xlim(0.,10)
            ax.set_ylim(1e-4,1e1)
            
            ax.set_yscale("log")
    
            ax.errorbar(self.x/(1000*nc.pc), self.data, self.err, label="data",\
                        markerfacecolor='C3',markeredgecolor='C3', marker='o',linestyle='', ecolor = 'C3')

            ax.plot(self.x_beam/(1000*nc.pc), self.beam*self.data[0]/self.beam[0], label="beam", linestyle="-", color="gray")
            
            if "name" in self.params_obs:
                ax.set_title(self.params_obs["name"], size=size)
            
        elif ax is not None:
                
            ax.errorbar(self.x/(1000*nc.pc), self.data, self.err, label="beam",\
                         marker='o',linestyle='')

            ax.plot(self.x_beam/(1000*nc.pc), self.beam*self.data[0]/self.beam[0], label="beam", linestyle="--", color="gray")
            
        return
    
    def print_values(self):
        """
        Print the values of the data
        Returns
        -------
        None.
        """

        print(self.params_obs)
        
        return