# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 18:12:43 2021

@author: anna
"""


import numpy as np
import matplotlib.pyplot as plt
import natconst as nc
import os 

import mydir
import emcee
import corner

from scipy.interpolate import interp1d
import matplotlib



matplotlib.rcParams.update({
        "font.size": 16.0,
        "font.family": 'sans-serif',
#        "font.sans-serif": ['Helvetica'],
        "axes.titlesize": 16.,
        "axes.labelsize": 16.,
        "xtick.labelsize": 16.,
        "ytick.labelsize": 16.,
        "xtick.major.size": 6.0,
        "ytick.major.size": 6.0,
        "xtick.minor.size": 3.0,
        "ytick.minor.size": 3.0,
        "legend.fontsize": 13.0,
        "legend.frameon": False
 #       "figure.dpi": 200,
 #       "savefig.dpi": 200,
#        "text.usetex": True
})



plot1 = False  # emission chains
plot2 = True      # corners
plot3 = False   # violins
plot4 = False  # final trends



params = dict([("DM_model", "NFW"),
           ("beta", 1.0), 
           ("SFR", 50.),
           ("f_esc_ion", 0.), 
           ("f_esc_FUV", 0.), 
           ("v_c", 200.),
           ("redshift", 5.),
           ("Zeta", 1.0),
           ("alfa", 1.0),
           ("R_in", 0.3)])    



if plot1:
    """
    # emission chains
    """
    
    from mcmc import get_emission_fast, get_other_params
    from mcmc import h, grid
    
    from load_data import obs_data_list

    folder_data = "data_emcee"
    
    folder_plot = "plot_emcee"

    nwalkers= 96
    nsteps = 1e4

    sample_step = int(40 * (nsteps/1e4))
    walker_step = int(12 * (nwalkers/96))


    fig, axs = plt.subplots(4, 2, sharey=True, sharex=True, figsize=(1.3*8.27,1.3*12.))

    for ax in axs:
        ax.set_xlabel("b [kpc]")
        ax.set_ylabel(r"flux [mJy/arcsec$^2$]")
        ax.set_yscale('log')
        ax.set_ylim((0.01,12))    
        ax.set_xlim((0.3,10))
    


    for data,  data_counter in zip(obs_data_list, range(len(obs_data_list))):
        

        filename = "{}_{:.0f}".format(data.params_obs["name_short"], nsteps)

    
        path = os.path.join(mydir.data_dir, folder_data, "{}.h5".format(filename))

        reader = emcee.backends.HDFBackend(path)

        samples = reader.get_chain()
        samples_flat = reader.get_chain(flat=True)

        ndim = 3


        other_params = get_other_params(data.params_obs["redshift"], data.params_obs["line_FWHM"])
    
        beam_interp = np.interp(h, data.x_beam/1e3/nc.pc, data.beam, right=0.)
    
        beam_interp[beam_interp<0.] = 0.
    
        beam_func = interp1d(h, beam_interp, \
                         fill_value = (beam_interp[0], 0.), bounds_error=False)
    
        beam_2d = beam_func(np.sqrt(grid[0]**2 + grid[1]**2))
        f_beam = np.fft.fft2(beam_2d)

                
        counter = 0
        for  walker in samples[::sample_step]:
            for theta in walker[::walker_step]:
                counter += 1
                print("computing emission for theta =", theta, "\t", "iteration number {}/{}"\
                  .format(counter, int( nsteps*nwalkers/sample_step/walker_step)))

                #if theta[0]>1.15:# and theta[1]>50.:

                intensity = get_emission_fast(theta, data, other_params, h, grid, f_beam)

                ax[data_counter].plot(h, intensity, alpha=0.1, color="gray")

        
        alpine = ax[data_counter].errorbar(data.x/(1000*nc.pc), data.data, yerr=data.err, \
            markerfacecolor='maroon',markeredgecolor='maroon', marker='o',\
            linestyle='', ecolor = 'maroon')
    

        ax[data_counter].legend(loc="lower center", ncol=8, fontsize="small")
        
            
        plt.show()
            
    
if plot2:
    """
    #corners
    """
    from load_data import obs_data_list
    data = obs_data_list[0]
    
    folder_data = "data_emcee"
    
    folder_plot = "plot_emcee"

    nwalkers= 96
    nsteps = 1e4

    sample_step = int(40 * (nsteps/1e4))
    walker_step = int(12 * (nwalkers/96))
    
    ndim = 3
    
    labels = ["beta", "SFR", "v_c"]

    filename = "{}_{:.0f}".format(data.params_obs["name_short"], nsteps)

    path = os.path.join(mydir.data_dir, folder_data, "{}.h5".format(filename))

    reader = emcee.backends.HDFBackend(path)

    samples = reader.get_chain()
    samples_flat = reader.get_chain(flat=True)

    kwargs = dict(
            bins=50, smooth=0.9, labels=labels, #label_kwargs=dict(fontsize=14),
            #title_kwargs=dict(fontsize=14), 
            color='#0072C1',
            truth_color='C3', quantiles=[0.16, 0.5, 0.84],
            levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
            plot_density=False, plot_datapoints=True, fill_contours=True,
            max_n_ticks=3)

    fig,ax = plt.subplots(ndim,ndim,figsize=(1.3*8.27/2,1.3*12./3))

    corner.corner(samples_flat, **kwargs, fig=fig, labelpad=-0.5)

    axes = np.array(fig.axes).reshape((ndim, ndim))



    true_values = [None, data.params_obs["SFR"], data.params_obs["v_c"], None]
    
    err_ups = [None, data.params_obs["SFR"]+data.params_obs["SFR_err_up"], data.params_obs["v_c"]+data.params_obs["v_c_err_up"], None ]
    err_downs = [None, data.params_obs["SFR"]-data.params_obs["SFR_err_down"], data.params_obs["v_c"]-data.params_obs["v_c_err_down"], None ]
    
    # Loop over the diagonal
    for i in range(ndim):
        
        ax = axes[i, i]
        if i == 1 or i == 2:
            ax.axvline(true_values[i], color=kwargs["truth_color"], linestyle='--', linewidth=1.3)
            ax.axvspan(err_downs[i], err_ups[i], alpha=0.2, color='red')
    
        if i == 3:
            ax.set_xlim(-500,100)
    
        
    # Loop over the histograms
    for yi in range(ndim):
        for xi in range(yi):
            
            
            ax = axes[yi, xi]
            
            if yi == 3:
                ax.set_ylim(-500,200)
            
            
            if xi == 1 or xi == 2:
                ax.axvline(true_values[xi], color=kwargs["truth_color"], linestyle='--', linewidth=1.3)
                ax.axvspan(err_downs[xi], err_ups[xi], alpha=0.2, color='red')
                if xi == 1:
                    ax.set_xlim(1.,250.)
                if xi == 2:
                    ax.set_xlim(100.,450.)
    
    
            if yi == 1 or yi == 2:
                ax.axhline(true_values[yi], color=kwargs["truth_color"], linestyle='--', linewidth=1.3)
                ax.axhspan(err_downs[yi], err_ups[yi], alpha=0.2, color='red')
                if yi == 1:
                    ax.set_ylim(1.,250.)
                if yi == 2:
                    ax.set_ylim(100.,450.)

    plt.subplots_adjust(left=0.17,  # the left side of the subplots of the figure
                        right=0.98,  # the right side of the subplots of the figure
                        bottom=0.17,  # the bottom of the subplots of the figure
                        top=0.93,  # the top of the subplots of the figure
                        wspace=0.05,  # the amount of width reserved for space between subplots,
                        # expressed as a fraction of the average axis width
                        hspace=0.05)  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height