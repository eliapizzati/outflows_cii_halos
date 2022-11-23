
To make the gnedincooling package working, I have to run the following command in the _other_scripts_ folder:

_f2py3 -c -m gnedincooling frt_cf3m.F_

and then copy the output file in the main script folder.


Just a small documentation of the scripts in the main folder.

N.B. There are three places where I can store the mcmc data (they are too heavy to be
stored in the github repository). I have to change the path in the scripts to make them work.
The three places are: 
 - the mac laptop (keywork: "mac")
 - the linux laptop (keywork: "linux")
 - the quasar machine in leiden (keyword: "quasar"; folder /data2/pizzati/projects/...)
 - the github folder (keyword: "github")

On top of that, there are some environments that I have to activate to run the scripts.
These are: conda env outflows_new on the mac, outenv python3 on the quasar machine, and 
outflows_new (?) on the linux laptop.

There are two main workflows: the first one builds on many classes and it is suitable for 
every task that is not computationally intensive. The second one is a faster version of the pipeline 
used for the mcmc routine. It is not as flexible as the first one, but it is faster.

The first workflow is based on the following scripts:

 - model_classes.py: contains the classes used to build the model. 
                   It is the core of the pipeline.

 - data_classes.py: contains the classes used to store the observational data. 
                  It is the core of the pipeline.

 - cmb_suppression.py: contains the functions used to compute the cmb suppression.

 - radiation_fields.py: contains the functions used to compute the radiation fields.

 - load_data.py: loads the observational data.

 - _sol_modules.py_ (requires gnedincooling): contains the functions used to compute
            the solutions of the differential equations and find the n,v,T profiles

 - _post_sol_modules.py_: contains the functions used to compute
            the abundances and the final emission in CII (everything after the
            solution of the differential equations)

 - _main_pipeline_for_ALPINE.py_: contains the main function to run the pipeline. It is the main script
               to run the pipeline for ALPINE data

 - _main_profile_saver.py_: similar to main_pipeline_for_ALPINE but just for saving data files

 - _main_profile_reader.py_: script to read the profiles saved in main_profile_saver


The second workflow is based on the following scripts (all starting with _mcmc_):

 - _mcmc_fast_solver.py_ (requires gnedincooling): defines the system of differential equations

 - _mcmc_fast_emission.py_ (requires gnedincooling): contains the core function get_emission_fast that integrates
                 the system of differential equations and computes the emission

 - _mcmc_likelihood_and_priors.py_ (requires gnedincooling): contains the functions used to compute the likelihood and the priors

 - _mcmc_main.py_ (requires gnedincooling): contains the main function to run the mcmc routine

 - _mcmc_post_analysis.py_ (requires gnedincooling): contains the functions used to analyze the results of the mcmc routine

 - _mcmc_test_switchers.py_ (requires gnedincooling): creates the switchers to study emission as a function of the parameters


Other scritps that are doing other useful stuff:

 - _my_utils.py_ : contains some useful functions for the analysis of the data

 - _info.py_: contains the information about the parameters
 
 - _my_dir.py_: contains the paths to the data folders

 - _natconst.py_: contains the natural constants

 - _plot_config.py_: contains the configuration for the plots

Then there are the scripts that are specifically created to plot stuff
(e.g. for the paper, presentations, etc.):

 - _plotter_alpine.py_ (requires gnedincooling): plots the emission for alpine systems


 - _plotter_non_detections_analysis.py_ : plots the final figure for the paper 
                                where we analyze the non-detections 

 - _plotter_mcmc.py_ (requires gnedincooling): plots the results of the mcmc routine,
                            used for the plots related to mcmc in the paper

 - _plotter_paper.py_(requires gnedincooling): plots the first figures for the paper
                                       (i.e., NFW profile, cmb emission)

 - _plotter_ppt.py_(requires gnedincooling): plots the figures for the powerpoint presentation
                  (i.e., cc85,cooling function, cmb suppression theory, 
                   profiles, ionization)

 - _plotter_thesis_cosmo.py_ : used in the thesis to print some figures
                                for the introduction chapter

 - _plotter_thesis_model.py_ (requires gnedincooling): used in the thesis to print some figures
                                for the model chapter
 
 - _plotter_thesis_profiles.py_ (requires gnedincooling): used in the thesis to create some figures
                                for the model/results chapter

 - _plotter_vc_disk.py_(requires gnedincooling): analyzes the inclusion of the disk profiles in the graivty


Other scritps that I created and that are not part of the workflow are:

 - _cooling.py_ and _cooling_switchers.py_ (requires gnedincooling):  study the cooling function to check 
that everything is behaving well. The cooling function is defined in gnedincooling

 - _trial_spectrum.py_ (requires gnedincooling): tries to study the spectral shape of the CII line

 - _temp_fuji_alpine_comparison.py_: compares the two sets of data