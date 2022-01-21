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
from scipy.optimize import curve_fit
import matplotlib

from collections import namedtuple





matplotlib.rcParams.update({
        "font.size": 16.0,
        "font.family": 'sans-serif',
#        "font.sans-serif": ['Helvetica'],
        "axes.titlesize": 15.,
        "axes.labelsize": 14.,
        "xtick.labelsize": 14.,
        "ytick.labelsize": 14.,
        "xtick.major.size": 4.0,
        "ytick.major.size": 4.0,
        "xtick.minor.size": 2.0,
        "ytick.minor.size": 2.0,
        "legend.fontsize": 13.0,
        "legend.frameon": False
 #       "figure.dpi": 200,
 #       "savefig.dpi": 200,
#        "text.usetex": True
})


plot0 = False # single emission chain
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


if plot0:


    from mcmc import get_emission_fast, get_other_params
    from mcmc import h, grid

    from load_data import obs_data_list

    folder_data = "data_emcee"

    folder_plot = "plot_emcee"

    nwalkers = 96
    nsteps = 1e5

    thin = 10
    discard = 4000

    sample_step = int(12 * (nsteps / 1e4))
    walker_step = int(12 * (nwalkers / 96))

    fig, ax = plt.subplots(figsize=(1.3*8.27/1.6,1.2*12./2.2))

    ax.set_yscale('log')
    ax.set_ylim((0.003, 12))
    ax.set_xlim((0.3, 16))
    ax.tick_params(length=4, axis="both", which="major")
    ax.tick_params(length=2, axis="both", which="minor")

    ax.set_xlabel("b [kpc]")
    ax.set_ylabel('flux [mJy/arcsec$^2$]')

    data = obs_data_list[0]

    filename = "{}_{:.0f}_final".format(data.params_obs["name_short"], nsteps)

    path = os.path.join(mydir.data_dir, folder_data, "{}.h5".format(filename))

    reader = emcee.backends.HDFBackend(path)

    samples = reader.get_chain(thin=thin, discard=discard)
    samples_flat = reader.get_chain(flat=True, thin=thin, discard=discard)

    print("###################################################################")

    print("postprocessing an MCMC with the following params:")
    print("n steps = {}".format(nsteps))
    print("n walkers = {}".format(nwalkers))
    print("data object = {}".format(data.params_obs["name_short"]))
    print("filename = {}".format(filename))

    print("###################################################################")

    ndim = 3

    other_params = get_other_params(data.params_obs["redshift"], data.params_obs["line_FWHM"])

    beam_interp = np.interp(h, data.x_beam / 1e3 / nc.pc, data.beam, right=0.)

    beam_interp[beam_interp < 0.] = 0.

    beam_func = interp1d(h, beam_interp, \
                         fill_value=(beam_interp[0], 0.), bounds_error=False)

    beam_2d = beam_func(np.sqrt(grid[0] ** 2 + grid[1] ** 2))
    f_beam = np.fft.fft2(beam_2d)

    counter = 0
    for walker in samples[::sample_step]:
        for theta in walker[::walker_step]:
            counter += 1
            print("computing emission for theta =", theta, "\t", "iteration number {}/{}" \
                      .format(counter, int(nsteps * nwalkers / sample_step / walker_step)))

            # if theta[0]>1.15:# and theta[1]>50.:
            intensity = get_emission_fast(theta, data, other_params, h, grid, f_beam)

            ax.plot(h, intensity, alpha=0.1, color="gray")

        alpine = ax.errorbar(data.x / (1000 * nc.pc), data.data, yerr=data.err, \
                                                 markerfacecolor='maroon', markeredgecolor='maroon', marker='o', \
                                                 linestyle='', ecolor='maroon')

        # axs_flat[data_counter].legend(loc="lower center", ncol=8, fontsize="small")

    plt.subplots_adjust(left=0.15,  # the left side of the subplots of the figure
                        right=0.97,  # the right side of the subplots of the figure
                        bottom=0.25,  # the bottom of the subplots of the figure
                        top=0.9,  # the top of the subplots of the figure
                        wspace=0.05,  # the amount of width reserved for space between subplots,
                        # expressed as a fraction of the average axis width
                        hspace=0.27)  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height

    plt.show()


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
    nsteps = 1e5
    
    thin = 10
    discard = 4000


    sample_step = int(12 * (nsteps/1e5))
    walker_step = int(12 * (nwalkers/96))


    fig, axs = plt.subplots(4, 2, sharey=True, sharex=True, figsize=(1.3*8.27,1.3*12.))

    axs_flat = axs.flatten()
    for ax in axs_flat:
         #ax.set_xlabel("b [kpc]")
         #ax.set_ylabel(r"flux [mJy/arcsec$^2$]")
         ax.set_yscale('log')
         ax.set_ylim((0.003,12))
         ax.set_xlim((0.3,16))

    axs[3,0].set_xlabel("b [kpc]")
    axs[3,1].set_xlabel("b [kpc]")
    fig.text(0.011, 0.53, 'flux [mJy/arcsec$^2$]', va='center', rotation='vertical')

    


    for data,  data_counter in zip(obs_data_list, range(len(obs_data_list))):

        axs_flat[data_counter].set_title("{}".format(data.params_obs["name_short"]))

        filename = "{}_{:.0f}_final".format(data.params_obs["name_short"], nsteps)

    
        path = os.path.join(mydir.data_dir, folder_data, "{}.h5".format(filename))

        reader = emcee.backends.HDFBackend(path)

        samples = reader.get_chain(thin=thin, discard=discard)
        samples_flat = reader.get_chain(flat=True, thin=thin, discard=discard)

        print("###################################################################")

        print("postprocessing an MCMC with the following params:")
        print("n steps = {}".format(nsteps))
        print("n walkers = {}".format(nwalkers))
        print("data object = {}".format(data.params_obs["name_short"]))
        print("filename = {}".format(filename))

        print("###################################################################")

        ndim = 3

        other_params = get_other_params(data.params_obs["redshift"], data.params_obs["line_FWHM"])

        beam_interp = np.interp(h, data.x_beam / 1e3 / nc.pc, data.beam, right=0.)

        beam_interp[beam_interp < 0.] = 0.

        beam_func = interp1d(h, beam_interp, \
                             fill_value=(beam_interp[0], 0.), bounds_error=False)

        beam_2d = beam_func(np.sqrt(grid[0] ** 2 + grid[1] ** 2))
        f_beam = np.fft.fft2(beam_2d)

        counter = 0
        for  walker in samples[::sample_step]:
            for theta in walker[::walker_step]:
                counter += 1
                print("computing emission for theta =", theta, "\t", "iteration number {}/{}"\
                  .format(counter, int( nsteps*nwalkers/sample_step/walker_step)))

                #if theta[0]>1.15:# and theta[1]>50.:
                intensity = get_emission_fast(theta, data, other_params, h, grid, f_beam)

                axs_flat[data_counter].plot(h, intensity, alpha=0.1, color="gray")


        alpine = axs_flat[data_counter].errorbar(data.x/(1000*nc.pc), data.data, yerr=data.err, \
            markerfacecolor='maroon',markeredgecolor='maroon', marker='o',\
            linestyle='', ecolor = 'maroon')
    

        #axs_flat[data_counter].legend(loc="lower center", ncol=8, fontsize="small")

    plt.subplots_adjust(left=0.1,  # the left side of the subplots of the figure
                        right=0.98,  # the right side of the subplots of the figure
                        bottom=0.08,  # the bottom of the subplots of the figure
                        top=0.95,  # the top of the subplots of the figure
                        wspace=0.05,  # the amount of width reserved for space between subplots,
                        # expressed as a fraction of the average axis width
                        hspace=0.27)  # the amount of height reserved for space between subplots,
    # expressed as a fraction of the average axis height

    plt.show()
            
    
if plot2:
    """
    #corners
    """
    from load_data import obs_data_list
    
    data = obs_data_list[7]
    
    folder_data = "data_emcee"
    
    folder_plot = "plot_emcee"

    nwalkers= 96
    nsteps = 1e5
    
    thin = 10
    discard =5000

    
    ndim = 3
    
    labels = [r"log( $\eta$ )", "SFR [M$_\odot$ yr$^{-1}$]", r"$v_c^{\mathrm{(vir)}}$ [km s$^{-1}$]"]

    filename = "{}_{:.0f}_final".format(data.params_obs["name_short"], nsteps)

    path = os.path.join(mydir.data_dir, folder_data, "{}.h5".format(filename))

    reader = emcee.backends.HDFBackend(path)

    samples = reader.get_chain(thin=thin, discard=discard)
    samples_flat = reader.get_chain(flat=True, thin=thin, discard=discard)

    #tau = reader.get_autocorr_time(thin=thin, discard=discard)

    print("###################################################################")

    print("postprocessing an MCMC with the following params:")
    print("n steps = {}".format(nsteps))
    print("n walkers = {}".format(nwalkers))
    print("data object = {}".format(data.params_obs["name_short"]))
    print("filename = {}".format(filename))

    print("###################################################################")

    kwargs = dict(
            bins=50, smooth=0.9, labels=labels, label_kwargs=dict(labelpad=0),
            #title_kwargs=dict(fontsize=14), 
            color='#0072C1',
            truth_color='C3', quantiles=[0.16, 0.5, 0.84],
            levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
            plot_density=False, plot_datapoints=True, fill_contours=True,
            max_n_ticks=3)

    fig,ax = plt.subplots(ndim,ndim,figsize=(1.3*8.27/1.6,1.3*12./2.2))

    corner.corner(samples_flat, **kwargs, fig=fig)

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
    
        summary = namedtuple('summary', ['median', 'lower', 'upper', 'string'])
        quantiles=kwargs['quantiles']
        quants_to_compute = np.array([quantiles[0], 0.5, quantiles[2]])
        quants = np.percentile(samples_flat.T[i], quants_to_compute * 100)
        summary.median = quants[1]
        summary.plus = quants[2] - summary.median
        summary.minus = summary.median - quants[0]
    
        if i == 0:
            ax.set_title(r"$ %.2f_{-%.2f}^{+%.2f} $ "%(summary.median, summary.minus, summary.plus))
    
        if i == 1:
            ax.set_title(r"$ %.0f_{-%.0f}^{+%.0f} \ \mathrm{M}_\odot\mathrm{yr}^{-1}$"%(summary.median, summary.minus, summary.plus))
           
        if i == 2:
            ax.set_title(r"$ %.0f_{-%.0f}^{+%.0f} \ \mathrm{km s}^{-1}$"%(summary.median, summary.minus, summary.plus))


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

    fig.text(0.705, 0.825, "results for""\n""VC5110377875",\
             fontsize = 16, bbox = {'facecolor': 'lightgray', 'alpha': 0.8, \
            'boxstyle': "round", 'edgecolor': 'none', 'pad': 0.5}, \
            horizontalalignment = 'center', verticalalignment = 'center')

    #fig.text(0.71, 0.77, "{}".format(data.params_obs["name_short"]), ha='center', size=19)

    # plt.subplots_adjust(left=0.13,  # the left side of the subplots of the figure
    #                     right=0.98,  # the right side of the subplots of the figure
    #                     bottom=0.13,  # the bottom of the subplots of the figure
    #                     top=0.93,  # the top of the subplots of the figure
    #                     wspace=0.05,  # the amount of width reserved for space between subplots,
    #                     # expressed as a fraction of the average axis width
    #                     hspace=0.05)  # the amount of height reserved for space between subplots,
    # # expressed as a fraction of the average axis height
    #
    plt.show()
    
    
if plot3:
    """
    # violins
    """
    from load_data import obs_data_list

    datas = []
    names_plot = []

    for data in obs_data_list:

        datas.append(data)
        names_plot.append(data.params_obs["name_short"])

    
    folder_data = "data_emcee"

    folder_plot = "plot_emcee"

    thin = 300
    discard =4000

    nwalkers = 96
    nsteps = 1e5

    log_betas = []
    log_probs = []
    #sfrs = []
    #vcs = []

    for data in datas[::]:

        filename = "{}_{:.0f}_final".format(data.params_obs["name_short"], nsteps)
    
        print("###################################################################")
        
        print("postprocessing an MCMC with the following params:")
        print("n steps = {}".format( nsteps))
        print("n walkers = {}".format( nwalkers))
        print("data object = {}".format( data.params_obs["name_short"]))
        print("filename = {}".format( filename))
                     
        print("###################################################################")
    
    
        path = os.path.join(mydir.data_dir, folder_data, "{}.h5".format(filename))
        
        reader = emcee.backends.HDFBackend(path)
        
        samples_flat = reader.get_chain(flat=True, thin=thin, discard=discard)
        
        log_prob_samples = reader.get_log_prob(flat=True, thin=thin, discard=discard)
    
        log_prior_samples = reader.get_blobs(flat=True, thin=thin, discard=discard)

        all_samples = np.concatenate(
                (samples_flat, log_prob_samples[:, None],\
                log_prior_samples[:, None]), \
                axis=1)
        
        mask = np.logical_and(all_samples[:,3]>-20, all_samples[:,0]<1.3)
        log_betas.append(all_samples[:,0][mask])
        #sfrs.append(all_samples[:,1][all_samples[:,3]>-20])
        #vcs.append(all_samples[:,2][all_samples[:,3]>-20])
        log_probs.append(all_samples[:,3][mask]-all_samples[:,4][mask])

    betas = 10**np.asarray(log_betas)
    logprobs = np.asarray(log_probs)
    #sfrs = np.asarray(sfrs)
    #vcs = np.asarray(vcs)

    ndim = 3
        
    fig, [ax_betas, ax_probs] = plt.subplots(1,2, figsize=(1.3*8.27,1.3*3.2), sharey=True)
    
    parts = ax_betas.violinplot(betas, positions=None, vert=False, widths=0.8, showmeans=False, showextrema=False, showmedians=False,\
                  quantiles=[[0.16,0.84]]*len(betas), \
                   points=50, bw_method=None)
    
    i=0
    for pc in parts['bodies']:
        pc.set_facecolor('C{}'.format(i))
        pc.set_edgecolor('gray')
        pc.set_alpha(1)
        i=i+1
    
    parts = ax_probs.violinplot(logprobs, positions=None, vert=False, widths=0.8, showmeans=False, showextrema=False, showmedians=False,\
                  quantiles=[[0.16,0.84]]*len(betas), \
                   points=50, bw_method=None)
    
    i=0
    for pc in parts['bodies']:
        pc.set_facecolor('C{}'.format(i))
        pc.set_edgecolor('gray')
        pc.set_alpha(1)
        i=i+1
    
    ax_betas.set_xlabel("$\eta$")
    ax_betas.set_xlim(1.,20.)
    
    ax_probs.set_xlim(-15,1)
    ax_probs.set_xlabel(r"log( $\mathcal{L}$(d|$\theta$,m) )")
    
    ax_betas.set_yticks(np.arange(1,len(names_plot)+1))
    ax_betas.set_yticklabels(names_plot, rotation=35.)
    
    ax_probs.axvline(-8, color = "gray", linestyle="--")
    
    plt.subplots_adjust(left = 0.17,  # the left side of the subplots of the figure
        right = 0.95,   # the right side of the subplots of the figure
        bottom = 0.17,  # the bottom of the subplots of the figure
        top = 0.95,     # the top of the subplots of the figure
        wspace = 0.1,  # the amount of width reserved for space between subplots,
        # expressed as a fraction of the average axis width
        hspace = 0.1)  # the amount of height reserved for space between subplots,
                      # expressed as a fraction of the average axis height
    
    

if plot4:
    """
    # final trends
    """    
    from load_data import obs_data_list

    
    def power(x, a,b):
        return a*x**b
    
    def muratov_fit(x):
        return 3.6*x**(-0.35)

    def fluetsch_agn_fit(x):
        return 0.85*x**(0.76)
    
    
    likelihood_means = []
    beta_means = []
    mstars = []
    sfrs = []
    vcs = []
    sigma_betas_up = []
    sigma_betas_down = []
    sigma_mstars_up = []
    sigma_mstars_down = []
    sigma_sfrs_up = []
    sigma_sfrs_down = []
    

    datas = []
    names_plot = []

    for data in obs_data_list:

        datas.append(data)
        names_plot.append(data.params_obs["name_short"])

    
    folder_data = "data_emcee"

    folder_plot = "plot_emcee"

    thin = 300
    discard =4000

    nwalkers = 96
    nsteps = 1e5

    log_betas = []
    log_probs = []
    #sfrs = []
    #vcs = []

    for data in datas[::]:

        filename = "{}_{:.0f}_final".format(data.params_obs["name_short"], nsteps)
    
        print("###################################################################")
        
        print("postprocessing an MCMC with the following params:")
        print("n steps = {}".format( nsteps))
        print("n walkers = {}".format( nwalkers))
        print("data object = {}".format( data.params_obs["name_short"]))
        print("filename = {}".format( filename))
                     
        print("###################################################################")
    
    
        path = os.path.join(mydir.data_dir, folder_data, "{}.h5".format(filename))
        
        reader = emcee.backends.HDFBackend(path)
        
        samples_flat = reader.get_chain(flat=True, thin=thin, discard=discard)
        
        log_prob_samples = reader.get_log_prob(flat=True, thin=thin, discard=discard)

        log_prior_samples = reader.get_blobs(flat=True, thin=thin, discard=discard)

        all_samples = np.concatenate(
            (samples_flat, log_prob_samples[:, None], \
             log_prior_samples[:, None]), \
            axis=1)

        mask = np.logical_and(all_samples[:, 3] > -20, all_samples[:, 0] < 1.3)
        log_betas.append(all_samples[:, 0][mask])
        # sfrs.append(all_samples[:,1][all_samples[:,3]>-20])
        # vcs.append(all_samples[:,2][all_samples[:,3]>-20])
        log_probs.append(all_samples[:, 3][mask] - all_samples[:, 4][mask])


    betas = 10**np.asarray(log_betas)
    logprobs = np.asarray(log_probs)
    #sfrs = np.asarray(sfrs)
    #vcs = np.asarray(vcs)


    
    for data, beta, likelihood in zip(datas, betas, logprobs):
        print(data.params_obs["name_short"], np.median(beta),\
              np.percentile(beta, 84)-np.median(beta), np.median(beta)-np.percentile(beta, 16),\
              likelihood.max())
        #data.params_obs.update(beta_best_fit=beta.mean())
        #data.params_obs.update(beta_uncertainty=np.std(beta))    
        #data.params_obs.update(likelihood_best_fit=likelihood.mean())

        likelihood_means.append(likelihood.mean())
        beta_means.append(np.median(beta))
        mstars.append(data.params_obs["M_star"]/1e10)
        sfrs.append(data.params_obs["SFR"])
        vcs.append(data.params_obs["v_c"])
        sigma_betas_up.append(np.percentile(beta, 84)-np.median(beta))
        sigma_betas_down.append(np.median(beta)-np.percentile(beta,16))
        sigma_mstars_up.append((data.params_obs["M_star_err_up"])/1e10)
        sigma_mstars_down.append((data.params_obs["M_star_err_down"])/1e10)
        sigma_sfrs_up.append((data.params_obs["SFR_err_up"]))
        sigma_sfrs_down.append((data.params_obs["SFR_err_down"]))
      
        
    likelihood_means = np.asarray(likelihood_means)
    beta_means = np.asarray(beta_means)
    mstars = np.asarray(mstars)
    sfrs = np.asarray(sfrs)
    vcs = np.asarray(vcs)
    sigma_betas_up = np.asarray(sigma_betas_up)
    sigma_betas_down = np.asarray(sigma_betas_down)

    sigma_mstars_up = np.asarray(sigma_mstars_up)
    sigma_mstars_down = np.asarray(sigma_mstars_down)
    sigma_sfrs_up = np.asarray(sigma_sfrs_up)
    sigma_sfrs_down = np.asarray(sigma_sfrs_down)

    
    #beta_means = np.asarray([5.5,8.4,6.4,5.3,5.9,4.0,8.2,4.3])

    

    fig, [ax_vc, ax_sfr] = plt.subplots(1,2,sharey=True,figsize=(1.3*8.27,1.3*3.2))
    
    cmap_rend_col = matplotlib.cm.get_cmap('viridis_r')
    
    SFR_min = 5.
    SFR_max = 110.
    
    norm = matplotlib.colors.Normalize(vmin=SFR_min, vmax=SFR_max)

    mask = [True,False,True,True,True,True,True,True]

    ax_vc.scatter(mstars[mask], beta_means[mask],\
               marker='o', color = cmap_rend_col((sfrs[mask]-SFR_min)/(SFR_max-SFR_min)),\
               s=500*(1./abs(likelihood_means[mask])) )
    

    ax_vc.errorbar(mstars[mask], beta_means[mask], \
                     [sigma_betas_down[mask],sigma_betas_up[mask]],\
                     [sigma_mstars_down[mask],sigma_mstars_up[mask]],
                      ecolor="gray", linestyle="", barsabove=False, alpha=0.2)

    
    from matplotlib.ticker import ScalarFormatter
    for axis in [ax_vc.xaxis, ax_vc.yaxis]:
        axis.set_major_formatter(ScalarFormatter())
        
    yticks = [2.,3.,4.,6.,10.,20.]
    xticks_mstar = [0.3,0.4,0.6,1.,2.0]
    xticks_sfr = [10.,20.,50.,100.]
    
    #ax.set_ylim(1.3,6.0)
    ax_vc.set_xlabel(r"M$_{\rm star}$ [10$^{10}$ M$_\odot$]")
    ax_vc.set_ylabel("$\eta$", labelpad=-4.)
    ax_vc.set_xscale("log")
    ax_vc.set_yscale("log")
    ax_vc.set_yticks(yticks)
    ax_vc.set_yticklabels(yticks)
    ax_vc.set_xticks(xticks_mstar)
    ax_vc.set_xticklabels(xticks_mstar)
    
    cax_vc = plt.axes([0.07, 0.17, 0.015,0.78])
    
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, cax=cax_vc, orientation='vertical')
    cb.set_label("SFR [M$_\odot$ yr$^{-1}$]", rotation=90., labelpad=-8)
   
    cax_vc.yaxis.set_ticks_position('left')
    cax_vc.yaxis.set_label_position('left')

    rel_error_betas = (sigma_betas_down[mask]+sigma_betas_up[mask])/2/beta_means[mask]
    rel_error_mstars = (sigma_mstars_down[mask]+sigma_mstars_up[mask])/2/mstars[mask]

    #mask = beta_means < 10.
    popt, pcov = curve_fit(power, mstars[mask], beta_means[mask], [3.6,-0.35],np.sqrt(rel_error_betas**2 + rel_error_mstars**2)*beta_means[mask])

        
    a, b = popt
    sigma_a, sigma_b = np.sqrt(pcov.diagonal())

    print("best fit for beta vs mstars: y = ({:.2f}\pm{:.2f}) * x ** ({:.2f}\pm{:.2f})".format(a,sigma_a,b,sigma_b))
    xaxis = np.linspace(0.3,2.)
    
    fit, = ax_vc.plot(xaxis, power(xaxis, a ,b), color="C0", linestyle="-")
    
    i=0
    while i < 50:
        a_sample = np.random.uniform(low=a-sigma_a, high = a+ sigma_a)
        b_sample = np.random.uniform(low=b-sigma_b, high = b+ sigma_b)

        a_sample_m = np.random.uniform(low=3.6-0.2*3.6, high = 3.6+ 0.2*3.6)
        b_sample_m = np.random.uniform(low=-0.35-0.02, high = -0.35+ 0.02)

        ax_vc.plot(xaxis, power(xaxis, a_sample_m ,b_sample_m), color="C1", alpha=0.015, linewidth=10.)
        ax_vc.plot(xaxis, power(xaxis, a_sample ,b_sample), color="C0", alpha=0.015, linewidth=5.)
        i+=1
        
    muratov, = ax_vc.plot(xaxis, muratov_fit(xaxis), color="C1", linestyle="-")
    
    ax_vc.legend([fit, muratov], ["This work", "Muratov+15"], loc="lower left")##
    

    
    cmap_rend_col = matplotlib.cm.get_cmap('inferno_r')

    
    v_c_min = 180.
    v_c_max = 260.
    
    norm = matplotlib.colors.Normalize(vmin=v_c_min, vmax=v_c_max)
    

    ax_sfr.scatter(sfrs[mask], beta_means[mask],\
               marker='o', color = cmap_rend_col((vcs[mask]-v_c_min)/(v_c_max-v_c_min)),\
               s=500*(1./abs(likelihood_means[mask])) )
    

    ax_sfr.errorbar(sfrs[mask], beta_means[mask], \
                     [sigma_betas_down[mask],sigma_betas_up[mask]],\
                     [sigma_sfrs_down[mask],sigma_sfrs_up[mask]],
                      ecolor="gray", linestyle="", barsabove=False, alpha=0.2)


    
    ax_sfr.set_xlabel(r"SFR [M$_\odot$ yr$^{-1}$]")
    #ax_sfr.set_ylabel("beta")
    ax_sfr.set_xscale("log")
    ax_sfr.set_yscale("log")
    ax_sfr.set_yticks(yticks)
    ax_sfr.set_yticklabels(yticks)
    ax_sfr.set_xticks(xticks_sfr)
    ax_sfr.set_xticklabels(xticks_sfr)

    
    cax_sfr = plt.axes([0.9, 0.17, 0.015,0.78])
    
    
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, cax=cax_sfr, orientation='vertical')
    cb.set_label(r"$v_c^{\mathrm{(vir)}}$ [km s$^{-1}$]", rotation=90., labelpad=2)

    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    
    #mask = beta_means < 10.
    mask = [True,False,True,True,True,True,True,True]

    rel_error_betas = (sigma_betas_down[mask] + sigma_betas_up[mask]) / 2 / beta_means[mask]
    rel_error_sfrs = (sigma_sfrs_down[mask] + sigma_sfrs_up[mask]) / 2 / sfrs[mask]

    print(np.sqrt(rel_error_betas ** 2 + rel_error_sfrs ** 2))
    # mask = beta_means < 10.
    popt, pcov = curve_fit(power, sfrs[mask], beta_means[mask], [3.6, -0.35],
                           np.sqrt(rel_error_betas ** 2 + rel_error_sfrs ** 2) * beta_means[mask])


    a, b = popt
    sigma_a, sigma_b = np.sqrt(pcov.diagonal())
    
    print("best fit for beta vs sfrs: y = ({:.2f}\pm{:.2f}) * x ** ({:.2f}\pm{:.2f})".format(a,sigma_a,b,sigma_b))
    xaxis = np.linspace(10.,150.)
    
    ax_sfr.plot(xaxis, power(xaxis, a ,b), color="C0", linestyle="-")
    
    i=0
    while i < 100:
        a_sample = np.random.uniform(low=a-sigma_a, high = a+ sigma_a)
        b_sample = np.random.uniform(low=b-sigma_b, high = b+ sigma_b)

        a_sample_m = np.random.uniform(low=3.6-0.2*3.6, high = 3.6+ 0.2*3.6)
        b_sample_m = np.random.uniform(low=-0.35-0.02, high = -0.35+ 0.02)

        ax_sfr.plot(xaxis, power(xaxis, a_sample ,b_sample), color="C0", alpha=0.015, linewidth=5.)
        i+=1

    ax_sfr.fill_between(xaxis, power(xaxis, a-sigma_a ,b-sigma_b), power(xaxis, a+sigma_a ,b+sigma_b), color="C0", alpha=0.1)



    plt.subplots_adjust(left = 0.16,  # the left side of the subplots of the figure
                        right = 0.885,   # the right side of the subplots of the figure
                        bottom = 0.17,  # the bottom of the subplots of the figure
                        top = 0.95,     # the top of the subplots of the figure
                        wspace = 0.05,  # the amount of width reserved for space between subplots,
                        # expressed as a fraction of the average axis width
                        hspace = 0.1)  # the amount of height reserved for space between subplots,
                      # expressed as a fraction of the average axis height


    plt.show()
