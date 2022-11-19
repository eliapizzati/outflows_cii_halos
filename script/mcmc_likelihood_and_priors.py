

import numpy as np
import time

import natconst as nc

from mcmc_fast_emission import get_emission_fast

def log_likelihood(theta, data, other_params, h, grid, f_beam):
    """
    Computes the log likelihood of the model given the data and the parameters

    Parameters
    ----------
    theta: array
        array of the parameters
    data:
        data to be fitted
    other_params: dict
        dictionary of the other parameters
    h: array

    grid: array

    f_beam: array

    Returns
    -------
    log_likelihood: float
        log likelihood of the model given the data and the parameters
    """

    intensity_convolved = get_emission_fast(theta, data, other_params, h, grid, f_beam)

    emission_predicted = np.interp(data.x / 1e3 / nc.pc, h, intensity_convolved)

    res = emission_predicted - data.data

    residuals = 2 * res / (data.err_down + data.err_up)

    chi2 = np.sum(residuals ** 2)

    return - chi2 / 2.


def log_prior_uniform(theta, data):
    """
    defines the priors for the set of parameters theta

    Parameters
    ==========
    theta: array
        parameters of the model
    data: class
        contains the observational data

    Returns
    =======
    priors value: float

    """

    log_beta, log_SFR, log_v_c = theta

    if 0.0 < log_beta < 1.5 and 0. < log_SFR < 2.4 and 2. < log_v_c < 2.6:
        return 0.0

    return -np.inf


def log_prior_gaussian(theta, data):
    """
    defines the priors for the set of parameters theta

    Parameters
    ==========
    theta: array
        parameters of the model
    data: class
        contains the observational data

    Returns
    =======
    priors value: float

    """

    log_beta, log_SFR, log_v_c = theta

    if 0.0 < log_beta < 1.5 and 0. < log_SFR < 2.4 and 2. < log_v_c < 2.6:
        prior = 0.0

        prior += - 2 * (log_SFR - data.params_obs["log_SFR"]) ** 2 / (
                data.params_obs["log_SFR_err_up"] + data.params_obs["log_SFR_err_down"]) ** 2
        prior += - 2 * (log_v_c - data.params_obs["log_v_c"]) ** 2 / (
                data.params_obs["log_v_c_err_up"] + data.params_obs["log_v_c_err_down"]) ** 2

        return prior
    else:
        return -np.inf


def log_prior_SFR_gaussian(theta, data):
    """
    defines the priors for the set of parameters theta

    Parameters
    ==========
    theta: array
        parameters of the model
    data: class
        contains the observational data

    Returns
    =======
    priors value: float

    """

    log_beta, log_SFR, log_v_c = theta

    if 0.0 < log_beta < 1.5 and 0. < log_SFR < 2.4 and 2. < log_v_c < 2.6:
        prior = 0.0

        prior += - 2 * (log_SFR - data.params_obs["log_SFR"]) ** 2 / (
                data.params_obs["log_SFR_err_up"] + data.params_obs["log_SFR_err_down"]) ** 2

        return prior
    else:
        return -np.inf


def log_probability(theta, data, other_params, h, grid, f_beam):
    """
    defines the probability as a combination of priors and likelihood

    Parameters
    ==========
    theta: array
        Parameters of the model
    data: class
        Data to be fitted
    other_params: dict
        Other parameters of the model
    h: array
        Array of the radial distance
    grid: array
        2d array of the radial distances (created with np.meshgrid)
    f_beam: array
        2d fft array of the beam profile

    Returns
    =======
    probability: float
        Probability of the model given the data
    priors: float
        Priors of the model
    """

    lp = log_prior_gaussian(theta, data)

    if not np.isfinite(lp):
        return -np.inf, -np.inf

    return lp + log_likelihood(theta, data, other_params, h, grid, f_beam), lp

