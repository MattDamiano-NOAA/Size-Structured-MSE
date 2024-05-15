# Size-Structured-MSE

# This repository contains the project files for the size-structured management strategy evaluation (MSE) tool that was developed for a NOAA Marine Fisheries Initiative (MARFIN) project during 2019-2023.

# The original goals of this project were to test an array of management strategies against commercial and recreational objectives for two Southeast marine stocks: South Atlantic black sea bass and Atlantic cobia, forecast quantities of interest such as catches by fleet, recruitment, abundance, and spawning stock biomass (SSB), and evaluate tradeoffs where they occur.

# This repository contains R files for the generalized size-structured operating models developed by Damiano et al. (2024*) and conditioned to mimic the most recent stock assessment for black sea bass (SEDAR 2023), and the size-structured stock assessment ADMB file developed by Cao et al. (2017). 

*https://doi.org/10.1016/j.fishres.2024.107028 (reach out for a copy of the paper if so desired)

# Model files and function:

# R files:
# Projection.R - central MSE file, calls all other files to run operating model and assessment model, run projections
# Functions.R - contains all operating model functions
# Input.data.R - sets up inputs to be used in functions to condition the operating model
# Data.gen.R - runs operating model functions, generates .dat files for the assessment model to fit
# Set.pars.R - allows user to turn priors/penalties on or off for assessment model-estimated quantities 
# Input.data.update.R - updates certain vectors, matrices, and arrays of inputs for projections
# Mse.proj.R - contains key function for mse projections: projects population dynamics, catches, and surveys
# Data.gen.proj.R - similarly, re-runs certain functions and updates .dat files for assessments in projections

# AD Model Builder files:
# NSLSAP01: .tpl file containing the size-structured assessment model. Must be compiled externally using a program like EMACS.

# For a complete description of how the operating models work, see Damiano et al. (2023) in Fisheries Research (in review), and for a complete description of the assessment model, see Cao et al.'s paper: Improving assessment of Pandalus stocks using a seasonal, size-structured assessment model with environmental variables. Part I: Model description and application https://doi.org/10.1139/cjfas-2016-0020

# Model notes and updates will be catalogued in this read.me 

This repository is a scientific product and is not official communication of the National Oceanic and
Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is
provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of
Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed
by all applicable Federal law. Any reference to specific commercial products, processes, or services by service
mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or
favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a
DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by
DOC or the United States Government.”
