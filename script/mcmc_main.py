"""
This script is the core of the mcmc analysis.
It defines the likelihood function and the priors and starts the sampling routines.
"""

import os
import numpy as np
import scipy
import logging

from multiprocessing import Pool

import time
import emcee
from scipy.interpolate import interp1d


import natconst as nc

import my_dir
from load_data import obs_data_list, names, names_CII_halo, observational_data_fuji
from mcmc_fast_emission import get_emission_fast, get_other_params
from mcmc_likelihood_and_priors import log_likelihood, log_probability

import gnedincooling as gc

gc.frtinitcf(0, os.path.join(my_dir.script_dir, "input_data", "cf_table.I2.dat"))

# preliminar parameters (TO CHECK EVERY TIME)

stored_data_loc = "quasar" # can be either quasar (for the machine is quasar in leiden) or mac (for mac laptop) or github (for github)
                        # or linux (for linux laptop)
parallel = True
optimization = False
model = "old"

resume_old_chain = True

# selecting the object to fit

data_counter = int(input("which data object?"))
data = obs_data_list[data_counter]

nwalkers = 48
nsteps = int(input("number of steps?"))

filename = "{}_{:.0f}_updated_{}".format(data.params_obs["name_short"], nsteps, model)


if stored_data_loc == "quasar":
    path = os.path.join("/data2/pizzati/projects/outflows/data_mcmc", "{}.h5".format(filename))
elif stored_data_loc == "mac" or stored_data_loc == "local": #default local is the mac laptop
    path = os.path.join("/Users/eliapizzati/projects/outflows/data_mcmc", "{}.h5".format(filename))
elif stored_data_loc == "github":
    folder = "data_emcee"
    if not os.path.exists(os.path.join(my_dir.data_dir, folder)):
        os.mkdir(os.path.join(my_dir.data_dir, folder))
    path = os.path.join(my_dir.data_dir, folder, "{}.h5".format(filename))
else:
    raise ValueError("stored_data_loc not recognized")


# parameters for the integration part

rmax = 30
h_resol = 500
r_resol = 500

cut = 45.
integrator = "RK45"

# setting up the grid and the other observational parameters

h = np.linspace(0.3, rmax, h_resol)
h_ext = np.linspace(-rmax, rmax, 2 * h_resol)
grid = np.meshgrid(h_ext, h_ext)

redshift = data.params_obs["redshift"]

FWHM_vel = data.params_obs["line_FWHM"]

beam_interp = np.interp(h, data.x_beam / 1e3 / nc.pc, data.beam, right=0.)

beam_interp[beam_interp < 0.] = 0.

beam_func = interp1d(h, beam_interp, \
                     fill_value=(beam_interp[0], 0.), bounds_error=False)

beam_2d = beam_func(np.sqrt(grid[0] ** 2 + grid[1] ** 2))
f_beam = np.fft.fft2(beam_2d)

other_params = get_other_params(redshift, FWHM_vel, r_resol, cut, integrator)

# setting up the initial guess for the parameters (if the optimization is not performed)
if not optimization:
    beta_best_fits = [5.5, 7.5, 6.4, 5.3, 5.9, 4.0, 8.2, 4.3]
    data.params_obs.update(beta_best_fit=beta_best_fits[data_counter])


# setting up the log file

logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', \
                    datefmt='%m/%d/%Y %I:%M:%S %p', \
                    level=logging.INFO, \
                    handlers=[logging.FileHandler(os.path.join(my_dir.log_dir, "{}.log".format(filename))), \
                              logging.StreamHandler()])

logging.info("###################################################################")

logging.info("running an MCMC with the following params:")
logging.info("n steps = {}".format(nsteps))
logging.info("n walkers = {}".format(nwalkers))
logging.info("parallelization = {}".format(parallel))
logging.info("data object = {}".format(data.params_obs["name_short"]))
logging.info("filename = {}".format(filename))

logging.info("###################################################################")


if optimization:
    # setting up the initial guess for the parameters (if the optimization is performed)

    # find optimal starting points for each walker
    chi2_func = lambda *args: -2 * log_probability(*args)[0]

    bounds = [(0.1, 1.), (0.,2.4), (2.,2.6)]


    args = (data, other_params, h, grid, f_beam)

    seed = 125546
    rand = np.random.RandomState(seed)

    result_opt = scipy.optimize.differential_evolution(chi2_func, bounds=bounds, popsize=25, recombination=0.7,
                                                       disp=True, polish=True, args=args, seed= rand)

    print(result_opt.x)
    print(result_opt.success)
    print("chi2 max", -2 * log_likelihood(result_opt.x, *args))

    ndim = len(result_opt.x)

    pos = [
        [np.clip(result_opt.x[i] + 1e-2 * (bounds[i][1] - bounds[i][0]) * np.random.randn(1)[0], bounds[i][0],
                 bounds[i][1])
         for i in range(ndim)] for i in range(nwalkers)]

    print(pos)

else:

    theta_input = [np.log10(data.params_obs["beta_best_fit"]), data.params_obs["log_SFR"], data.params_obs["log_v_c"]]

    ndim = len(theta_input)

    pos = theta_input + np.asarray([0.2, 0.2, 0.2]) * np.random.randn(nwalkers, ndim)

    # pos += np.asarray([0.5, 0., 0.]) * np.random.rand(nwalkers, ndim)

    pos[2][pos[2] < 1.] = 1. + np.abs(pos[2][pos[2] < 1.])

    pos[pos < 0.] = 1e-3

if resume_old_chain:
    ndim = 3
    pos = None

    if optimization == True:
        raise UserWarning("resuming old chain; optimization part will be ignored!")

backend = emcee.backends.HDFBackend(path)
if not resume_old_chain:
    backend.reset(nwalkers, ndim)


if parallel:

    with Pool() as pool:

        sampler = emcee.EnsembleSampler(nwalkers=nwalkers, ndim=ndim, log_prob_fn=log_probability, \
                                        args=(data, other_params, h, grid, f_beam), backend=backend, pool=pool)
        sampler.run_mcmc(pos, nsteps, progress=True)

else:

    sampler = emcee.EnsembleSampler(nwalkers=nwalkers, ndim=ndim, log_prob_fn=log_probability, \
                                    args=(data, other_params, h, grid, f_beam), backend=backend)
    sampler.run_mcmc(pos, nsteps, progress=True)

samples = sampler.get_chain()

print("Sampling done!")
print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))

print(
    "Mean autocorrelation time: {0:.3f} steps".format(
        np.mean(sampler.get_autocorr_time())))






