# Run full BBS occupancy analysis
library("brms")

# # The next three lines only need to be run (in order) once.
# source('/Users/jacobsocolar/Dropbox/Work/Code/Occupancy/biogeographicMSOM/data_prep/bbs_import.R')
# source('/Users/jacobsocolar/Dropbox/Work/Code/Occupancy/biogeographicMSOM/data_prep/range_maps.R')
# source('/Users/jacobsocolar/Dropbox/Work/Code/Occupancy/biogeographicMSOM/data_prep/format_for_occupancy.R')

##### Sort directories and core/thread info ####
if(grepl("bo1scm", getwd())) {
  # reminder: run this script from the *code* directory (which then gets switched- 
  # too much of a headache to run bash from data directory, so do it this way).
  # get cores/threads from bash
  cpu_info <- as.numeric(commandArgs(T))
  n_cores <- n_chains <- cpu_info[1]
  n_threads <- cpu_info[2]
  # set directory
  setwd("../data_bMSOM")
  code_dir <- "../biogeographicMSOM/"
} else {
  setwd('/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM')
  code_dir <- "/Users/jacob/Dropbox/Work/Code/Occupancy/biogeographicMSOM/"
  n_cores <- n_chains <- 4
  n_threads <- n_cores*2
}

# read data
flattened_data <- readRDS('flattened_data_2018.RDS')
flattened_data$elev_scaled <- as.numeric(scale(flattened_data$elev))
fd3 <- flattened_data[,c("N", "species", "site", "elev_scaled", "distance", "distance_scaled", "distance_transformed")]
fd3$elev_scaled2 <- fd3$elev_scaled^2
fd3$trials <- 5
# We consider two classes of model:  
# The first class, defined by `bbs_naive <- cmdstan_model(...)`, contains just elevation as an occupancy 
#   covariate. 
# The second class, defined by `bbs_distance <- cmdstan_model(...)`, treats (some a priori function of)
#   distance-from-range as a continuous covariate on occupancy. By transforming the distance-from-range
#   in the data object itself, the single stan model can handle a wide variety of distance-from-range
#   models.

prior1 <- c(set_prior("logistic(0, 1)", class = "Intercept"), 
            set_prior("normal(0, 10)", class = "sd"), 
            set_prior("normal(0, 10)", class = "sd", dpar = "zi"), 
            set_prior("normal(0, 10)", class = "b", dpar = "zi"))
# naive model (no clip)
bbs_naive2 <- brm(bf(N | trials(trials) ~ (1 |b| species), 
                 zi ~ elev_scaled + elev_scaled2 + 
                   (1 |b| species) + 
                   (0 + elev_scaled || species) +
                   (0 + elev_scaled2 || species)), 
                 data = fd3, family = zero_inflated_binomial(), 
                 prior = prior1, cores = 4, backend = 'cmdstanr')

saveRDS(bbs_naive2, "/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/stan_outputs/brms_bbs/naive2.RDS")

bbs_dist2 <- brm(bf(N | trials(trials) ~ (1 |b| species), 
                   zi ~ elev_scaled + elev_scaled2 + distance_transformed + 
                     (1 |b| species) + 
                     (0 + elev_scaled || species) +
                     (0 + elev_scaled2 || species) +
                     (0 + distance_transformed || species)), 
                data = fd3, family = zero_inflated_binomial(), 
                prior = prior1, backend = 'cmdstanr', cores = 4)
saveRDS(bbs_dist2, "/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/stan_outputs/brms_bbs/dist2.RDS")

##### clipping by distance #####
hist(fd3$distance)
hist(fd3$distance[fd3$N > 0])
max(fd3$distance[fd3$N > 0])
max(fd3$distance)
sum(fd3$distance > max(fd3$distance[fd3$N>0]))

fd4 <- fd3[fd3$distance < 400000, ]
bbs_dist_clip2 <- brm(bf(N | trials(trials) ~ (1 |b| species), 
                   zi ~ elev_scaled + elev_scaled2 + distance_transformed + 
                     (1 |b| species) + 
                     (0 + elev_scaled || species) +
                     (0 + elev_scaled2 || species) +
                     (0 + distance_transformed || species)), 
                data = fd4, family = zero_inflated_binomial(), 
                prior = prior1, backend = 'cmdstanr', cores = 4)
saveRDS(bbs_dist_clip2, "/Users/jacob/Dropbox/Work/Occupancy/biogeographicMSOM/stan_outputs/brms_bbs/dist_clip_400K2.RDS")
