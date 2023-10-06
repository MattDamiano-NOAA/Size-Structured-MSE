# Size-Structured-MSE

# This repository contains the project files for the size-structured management strategy evaluation (MSE) tool that was developed for a NOAA Marine Fisheries Initiative (MARFIN) project during 2019-2023.

# The original goals of this project were to test an array of management strategies against commercial and recreational objectives for two Southeast marine stocks: South Atlantic black sea bass and Atlantic cobia, forecast quantities of interest such as catches by fleet, recruitment, abundance, and spawning stock biomass (SSB), and evaluate tradeoffs where they occur.

# This repository contains R files for the generalized size-structured operating models developed by Damiano et al. (2023) and conditioned to mimic the most recent stock assessment for black sea bass (2023), and the size-structured stock assessment ADMB file developed by Cao et al. (2017). 

# Model files and function:
# Projection.R - central MSE file, calls all other files to run operating model and assessment model, run projections
# Functions.R - contains all operating model functions
# Input.data.R - sets up inputs to be used in functions to condition the operating model
# Data.gen.R - runs operating model functions, generates .dat files for the assessment model to fit
# Set.pars.R - allows user to turn priors/penalties on or off for assessment model-estimated quantities 
# Input.data.update.R - updates certain vectors, matrices, and arrays of inputs for projections
# Mse.proj.R - contains key function for mse projections: projects population dynamics, catches, and surveys
# Data.gen.proj.R - similarly, re-runs certain functions and updates .dat files for assessments in projections

# For a complete description of how the operating models work, see Damiano et al. (2023) in Fisheries Research (intent to submit Oct 2023), and for a complete description of the assessment model, see Cao et al.'s paper: Improving assessment of Pandalus stocks using a seasonal, size-structured assessment model with environmental variables. Part I: Model description and application https://doi.org/10.1139/cjfas-2016-0020

# Model notes and updates will be catalogued in this read.me 
