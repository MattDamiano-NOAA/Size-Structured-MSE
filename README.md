# Size-Structured-MSE for common dolphinfish (Coryphaena hippurus)

# This repository contains the project files for the size-structured management strategy evaluation (MSE) tool for Atlantic dolphinfish that is currently being developed for the National Marine Fisheries Service, Southeast Fisheries Science Center, and was partially funded by the NMFS-SeaGrant Population and Ecosystem Dynamics Fellowship during Sept. 2020 - March 2023. 

# This repository currently contains four R files
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
