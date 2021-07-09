# -*- coding: utf-8 -*-
"""
Created on Mon May 17 14:49:56 2021

@author: anna
"""



import os
import matplotlib.pyplot as plt
import numpy as np
import emcee

from matplotlib import ticker, cm, colors
from scipy.optimize import curve_fit

import plot_config as pltc
import mydir
import natconst as nc

from load_data import obs_data_list, names, names_CII_halo, names_wo_CII_halo,\
                      names_other,  observational_data_fuji, names_CII_halo_short


violin = False
fit = True

names_list = names_CII_halo_short

datas = []
names_plot = []

for data in obs_data_list:

    if data.params_obs["name_short"] in names_list:
        datas.append(data)
        names_plot.append(data.params_obs["name_short"])


folder_data = "data_emcee"
    
if not os.path.exists(os.path.join(mydir.data_dir, folder_data)):
    os.mkdir(os.path.join(mydir.data_dir, folder_data))


folder_plot = "plot_emcee"

if not os.path.exists(os.path.join(mydir.plot_dir, folder_plot)):
    os.mkdir(os.path.join(mydir.plot_dir, folder_plot))


nwalkers = 96
nsteps = 1e4

log_betas = []
log_probs = []
sfrs = []
vcs = []

for data in datas[::]:
    filename = "{}_{:.0f}_new_priors".format(data.params_obs["name_short"], nsteps)

    print("###################################################################")
    
    print("postprocessing an MCMC with the following params:")
    print("n steps = {}".format( nsteps))
    print("n walkers = {}".format( nwalkers))
    print("data object = {}".format( data.params_obs["name_short"]))
    print("filename = {}".format( filename))
                 
    print("###################################################################")


    path = os.path.join(mydir.data_dir, folder_data, "{}.h5".format(filename))
    
    reader = emcee.backends.HDFBackend(path)
    
    samples_flat = reader.get_chain(flat=True)
    
    log_prob_samples = reader.get_log_prob(flat=True)

    all_samples = np.concatenate(
            (samples_flat, log_prob_samples[:, None]),\
            #log_prior_samples[:, None]), \
            axis=1)

    log_betas.append(all_samples[:,0][all_samples[:,3]>-20])
    sfrs.append(all_samples[:,1][all_samples[:,3]>-20])
    vcs.append(all_samples[:,2][all_samples[:,3]>-20])
    log_probs.append(all_samples[:,3][all_samples[:,3]>-20])

betas = 10**np.asarray(log_betas)
logprobs = np.asarray(log_probs)
sfrs = np.asarray(sfrs)
vcs = np.asarray(vcs)

ndim = 3

if violin:
        
    fig, [ax_betas, ax_probs] = plt.subplots(1,2, figsize=(13,7), sharey=True)
    
    parts = ax_betas.violinplot(betas, positions=None, vert=False, widths=0.8, showmeans=False, showextrema=False, showmedians=False,\
                  #quantiles=[[0.25,0.75]]*len(betas), \
                   points=100, bw_method=None)
    
    i=0
    for pc in parts['bodies']:
        pc.set_facecolor('C{}'.format(i))
        pc.set_edgecolor('gray')
        pc.set_alpha(1)
        i=i+1
    
    parts = ax_probs.violinplot(logprobs, positions=None, vert=False, widths=0.8, showmeans=False, showextrema=False, showmedians=False,\
                  #quantiles=[[0.25,0.75]]*len(betas), \
                   points=100, bw_method=None)
    
    i=0
    for pc in parts['bodies']:
        pc.set_facecolor('C{}'.format(i))
        pc.set_edgecolor('gray')
        pc.set_alpha(1)
        i=i+1
    
    ax_betas.set_xlabel("beta")
    ax_betas.set_xlim(1.,20.)
    
    ax_probs.set_xlim(-15,1)
    ax_probs.set_xlabel("log probability")
    
    ax_betas.set_yticks(np.arange(1,len(names_plot)+1))
    ax_betas.set_yticklabels(names_plot)
    
    plt.subplots_adjust(left = 0.15,  # the left side of the subplots of the figure
        right = 0.95,   # the right side of the subplots of the figure
        bottom = 0.1,  # the bottom of the subplots of the figure
        top = 0.95,     # the top of the subplots of the figure
        wspace = 0.1,  # the amount of width reserved for space between subplots,
        # expressed as a fraction of the average axis width
        hspace = 0.1)  # the amount of height reserved for space between subplots,
                      # expressed as a fraction of the average axis height
    
    

    fig, [ax_sfrs, ax_vcs] = plt.subplots(1,2, figsize=(13,7), sharey=True)
    
    parts = ax_sfrs.violinplot(sfrs, positions=None, vert=False, widths=0.8, showmeans=False, showextrema=False, showmedians=False,\
                  #quantiles=[[0.25,0.75]]*len(betas), \
                   points=100, bw_method=None)
    
    i=0
    for pc in parts['bodies']:
        pc.set_facecolor('C{}'.format(i))
        pc.set_edgecolor('gray')
        pc.set_alpha(1)
        i=i+1
    
    parts = ax_vcs.violinplot(vcs, positions=None, vert=False, widths=0.8, showmeans=False, showextrema=False, showmedians=False,\
                  #quantiles=[[0.25,0.75]]*len(betas), \
                   points=100, bw_method=None)
    
    i=0
    for pc in parts['bodies']:
        pc.set_facecolor('C{}'.format(i))
        pc.set_edgecolor('gray')
        pc.set_alpha(1)
        i=i+1
    
    ax_sfrs.set_xlabel("SFR [msun/yr]")
    ax_sfrs.set_xlim(1.,250.)
    
    ax_vcs.set_xlim(100.,450.)
    ax_vcs.set_xlabel("vc [km/s]")
    
    ax_sfrs.set_yticks(np.arange(1,len(names_plot)+1))
    ax_sfrs.set_yticklabels(names_plot)
    
    plt.subplots_adjust(left = 0.15,  # the left side of the subplots of the figure
        right = 0.95,   # the right side of the subplots of the figure
        bottom = 0.1,  # the bottom of the subplots of the figure
        top = 0.95,     # the top of the subplots of the figure
        wspace = 0.1,  # the amount of width reserved for space between subplots,
        # expressed as a fraction of the average axis width
        hspace = 0.1)  # the amount of height reserved for space between subplots,
                      # expressed as a fraction of the average axis height
        


if fit:
    
    
    def power(x, a,b):
        return a*x**b
    
    def muratov_fit(x):
        return 3.6*x**(-0.35)
    
    
    likelihood_means = []
    beta_means = []
    mstars = []
    sfrs = []
    vcs = []
    sigma_betas = []
    sigma_mstars_up = []
    sigma_mstars_down = []
    sigma_sfrs_up = []
    sigma_sfrs_down = []
    

    
    for data, beta, likelihood in zip(datas, betas, logprobs):
        print(data.params_obs["name_short"], beta.mean(), 10**likelihood.mean())
        data.params_obs.update(beta_best_fit=beta.mean())
        data.params_obs.update(beta_uncertainty=np.std(beta))    
        data.params_obs.update(likelihood_best_fit=likelihood.mean())

        likelihood_means.append(likelihood.mean())
        beta_means.append(np.median(beta))
        mstars.append(data.params_obs["M_star"]/1e10)
        sfrs.append(data.params_obs["SFR"])
        vcs.append(data.params_obs["v_c"])
        sigma_betas.append(np.std(beta))
        sigma_mstars_up.append((data.params_obs["M_star_err_up"])/1e10)
        sigma_mstars_down.append((data.params_obs["M_star_err_down"])/1e10)
        sigma_sfrs_up.append((data.params_obs["SFR_err_up"]))
        sigma_sfrs_down.append((data.params_obs["SFR_err_down"]))
      
        
    likelihood_means = np.asarray(likelihood_means)
    beta_means = np.asarray(beta_means)
    mstars = np.asarray(mstars)
    sfrs = np.asarray(sfrs)
    vcs = np.asarray(vcs)
    sigma_betas = np.asarray(sigma_betas)
    sigma_mstars_up = np.asarray(sigma_mstars_up)
    sigma_mstars_down = np.asarray(sigma_mstars_down)
    sigma_sfrs_up = np.asarray(sigma_sfrs_up)
    sigma_sfrs_down = np.asarray(sigma_sfrs_down)

    
    # PLOT 4 -- BETAS VS MSTAR MCMC
    
    
    fig, ax = plt.subplots(figsize=(8.27,5.))
    
    cmap_rend_col = cm.get_cmap('viridis_r')
    
    SFR_min = 5.
    SFR_max = 110.
    
    norm = colors.Normalize(vmin=SFR_min, vmax=SFR_max)

    cii = ax.scatter(mstars, beta_means,\
               marker='o', color = cmap_rend_col((sfrs-SFR_min)/(SFR_max-SFR_min)),\
               s=500*(1./abs(likelihood_means)) )
    

    ax.errorbar(mstars, beta_means, \
                     [sigma_betas/2,sigma_betas/2],\
                     [sigma_mstars_down,sigma_mstars_up],
                      ecolor="gray", linestyle="", barsabove=False, alpha=0.2)

    
    from matplotlib.ticker import ScalarFormatter
    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_major_formatter(ScalarFormatter())
        
    yticks = [3.,4.,6.,10.,20.]
    xticks_mstar = [0.3,0.4,0.6,1.,2.0]
    xticks_sfr = [10.,20.,50.,100.]

    #ax.set_ylim(1.3,6.0)
    ax.set_xlabel("M_star [M_sun]/1e10")
    ax.set_ylabel("beta")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)
    ax.set_xticks(xticks_mstar)
    ax.set_xticklabels(xticks_mstar)
    
    plt.subplots_adjust(left = 0.15,  # the left side of the subplots of the figure
        right = 0.98,   # the right side of the subplots of the figure
        bottom = 0.15,  # the bottom of the subplots of the figure
        top = 0.9,     # the top of the subplots of the figure
        wspace = 0.1,  # the amount of width reserved for space between subplots,
        # expressed as a fraction of the average axis width
        hspace = 0.1)  # the amount of height reserved for space between subplots,
                      # expressed as a fraction of the average axis height
    
    
    cmap = cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, orientation='vertical')
    cb.set_label(r'SFR', rotation=90.)
    
    
    mask = beta_means < 10.
    popt, pcov = curve_fit(power, mstars[mask], beta_means[mask], [3.6,-0.35],np.sqrt(sigma_betas[mask]**2 + (sigma_mstars_down[mask]+sigma_mstars_up[mask])**2/4))

        
    a, b = popt
    sigma_a, sigma_b = np.sqrt(pcov.diagonal())

    print("best fit for beta vs mstars: y = ({:.2f}\pm{:.2f}) * x ** ({:.2f}\pm{:.2f})".format(a,sigma_a,b,sigma_b))
    xaxis = np.linspace(0.3,2.)
    
    fit, = ax.plot(xaxis, power(xaxis, a ,b), color="C0", linestyle="-")
    
    i=0
    while i < 50:
        a_sample = np.random.uniform(low=a-sigma_a, high = a+ sigma_a)
        b_sample = np.random.uniform(low=b-sigma_b, high = b+ sigma_b)

        a_sample_m = np.random.uniform(low=3.6-0.2*3.6, high = 3.6+ 0.2*3.6)
        b_sample_m = np.random.uniform(low=-0.35-0.02, high = -0.35+ 0.02)

        ax.plot(xaxis, power(xaxis, a_sample_m ,b_sample_m), color="C1", alpha=0.015, linewidth=10.)
        ax.plot(xaxis, power(xaxis, a_sample ,b_sample), color="C0", alpha=0.015, linewidth=5.)
        i+=1
        
    muratov, = ax.plot(xaxis, muratov_fit(xaxis), color="C1", linestyle="-")
    
    plt.legend([cii, fit, muratov], ["CII", "best-fit", "m+15"], loc="lower left")##
    
    # PLOT 5 -- betas vs SFR mcmc
    
    fig, ax = plt.subplots(figsize=(8.27,5.))
    
    cmap_rend_col = cm.get_cmap('viridis_r')
    
    v_c_min = 180.
    v_c_max = 260.
    
    norm = colors.Normalize(vmin=v_c_min, vmax=v_c_max)
    
    plt.show()
    
    cii = ax.scatter(sfrs, beta_means,\
               marker='o', color = cmap_rend_col((vcs-v_c_min)/(v_c_max-v_c_min)),\
               s=500*(1./abs(likelihood_means)) )
    

    ax.errorbar(sfrs, beta_means, \
                     [sigma_betas/2,sigma_betas/2],\
                     [sigma_sfrs_down,sigma_sfrs_up],
                      ecolor="gray", linestyle="", barsabove=False, alpha=0.2)


    
    ax.set_xlabel("SFR")
    ax.set_ylabel("beta")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)
    ax.set_xticks(xticks_sfr)
    ax.set_xticklabels(xticks_sfr)

    
    plt.subplots_adjust(left = 0.15,  # the left side of the subplots of the figure
        right = 0.98,   # the right side of the subplots of the figure
        bottom = 0.15,  # the bottom of the subplots of the figure
        top = 0.9,     # the top of the subplots of the figure
        wspace = 0.1,  # the amount of width reserved for space between subplots,
        # expressed as a fraction of the average axis width
        hspace = 0.1)  # the amount of height reserved for space between subplots,
                      # expressed as a fraction of the average axis height
    
    
    cmap = cm.ScalarMappable(norm=norm, cmap=cmap_rend_col)
    cmap.set_array([])
    
    cb = fig.colorbar(cmap, orientation='vertical')
    cb.set_label(r'    $v_c$', rotation=0.)
    
    
    plt.legend([cii,], ["CII",])##
    
    mask = beta_means < 10.
    popt, pcov = curve_fit(power, sfrs[mask], beta_means[mask], [3.6,-0.35],np.sqrt(sigma_betas[mask]**2 + (sigma_sfrs_down[mask]+sigma_sfrs_up[mask])**2/4))
    
    a, b = popt
    sigma_a, sigma_b = np.sqrt(pcov.diagonal())
    
    print("best fit for beta vs sfrs: y = ({:.2f}\pm{:.2f}) * x ** ({:.2f}\pm{:.2f})".format(a,sigma_a,b,sigma_b))
    xaxis = np.linspace(10.,150.)
    
    ax.plot(xaxis, power(xaxis, a ,b), color="C0", linestyle="-")
    
    ax.fill_between(xaxis, power(xaxis, a-sigma_a ,b-sigma_b), power(xaxis, a+sigma_a ,b+sigma_b), color="C0", alpha=0.05)

    plt.show()