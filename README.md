
To make the gnedincooling package working, I have to run the following command in the _other_scripts_ folder:

_f2py3 -c -m gnedincooling frt_cf3m.F_

and then copy the output file in the main script folder.


Just a small documentation of the scripts in the folder.



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

 - 

Other scritps that are doing other useful stuff:

 - _my_utils.py_ : contains some useful functions for the analysis of the data

 - _info.py_: contains the information about the parameters
 
 - _mydir.py_: contains the paths to the data folders

 - _natconst.py_: contains the natural constants

 - _plot_config.py_: contains the configuration for the plots


Other scritps that I created and that are not part of the workflow are:

 - _Cooling.py_ and _cooling_switchers.py_ study the cooling function to check 
that everything is behaving well. The cooling function is defined in gnedincooling

 - _trial_spectrum.py_: tries to study the spectral shape of the CII line

 - _temp_fuji_alpine_comparison.py_: compares the two sets of data