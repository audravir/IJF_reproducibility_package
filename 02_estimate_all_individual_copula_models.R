rm(list=ls(all=TRUE))

# MCMC size. The size used in paper M = 50000
# beware that M4 and M5 are computationally heavy, they might take a while...

M    = 50000

# R will print the status report *for each model*:
# - the % already completed
# - how much time it takes to run 100 iterations
# - how much time has passed in total

load('data/FXdata.Rdata')

source('individual_models_and_functions/M1_vector_dcc.R')
source('individual_models_and_functions/M2_vectordcc_tcop.R')
source('individual_models_and_functions/M3_dcc_HEAVY_scalar_tcop.R')
source('individual_models_and_functions/M4_XM.R')
source('individual_models_and_functions/M5_CAW_iw.R')
