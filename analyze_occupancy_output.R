# analyze output of run_occupancy_models.R
library(cmdstanr)
library(posterior)

setwd('/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM')


bbs_bufferclip_fit <- readRDS('bbs_bufferclip_fit.RDS')
bbs_naive_fit <- readRDS('bbs_naive_fit.RDS')

# get posterior Z probabilities


# get the log predictive density for the true values of the 3rd visit (distinguish Q = 0, Q = 1, and Q = 0 but third visit = 1)

# get posterior alpha diversity for each site, visualize result, and compare to naive alpha