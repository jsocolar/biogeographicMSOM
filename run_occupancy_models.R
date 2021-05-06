# Run full BBS occupancy analysis
library("cmdstanr")
library("posterior")
library("loo")
opencl_options <- list(stan_opencl = TRUE, opencl_platform_id = 0, opencl_device_id = 1)

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
    setwd('/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM')
    code_dir <- "/Users/jacobsocolar/Dropbox/Work/Code/Occupancy/biogeographicMSOM/"
    n_cores <- n_chains <- 4
    n_threads <- n_cores*2
}

# read data
flattened_data <- readRDS('flattened_data_2018.RDS')
fd2 <- flattened_data[order(flattened_data$Q, decreasing = T),]

# We consider two classes of model:  
# The first class, defined by `bbs_bufferclip <- cmdstan_model(...)`, contains no occupancy covariates. 
#   Distance relationships are handled by passing data that contain only a subset of species-point 
#   combinations. By passing only species-point combinations within some buffer distance of the species' 
#   ranges (e.g. 150 km), we clip the data to the buffered species ranges. By passing unclipped data
#   (equivalent to clipped data with an arbitrarily large buffer), we recover the naive model that lacks
#   range information.
# The second class, defined by `bbs_distance <- cmdstan_model(...)`, treats (some a priori function of)
#   distance-from-range as a continuous covariate on occupancy. By transforming the distance-from-range
#   in the data object itself, the single stan model can handle a wide variety of distance-from-range
#   models.

##### Class 1 #####
bbs_bufferclip <- cmdstan_model(paste0(code_dir, "stan_files/bufferclip_BBS_v2.stan"), cpp_options = opencl_options)

# naive model (no clip)
naive_data_stan <-list(
    n_visit = 5, 
    n_species = max(fd2$sp_id),
    n_tot = nrow(fd2),
    id_sp = fd2$sp_id,
    Q = fd2$Q,
    det_data = rowSums(fd2[,c("v1", "v2", "v3", "v4", "v5")])
)

bbs_naive_fit <- bbs_bufferclip$sample(data = naive_data_stan,
                                       iter_warmup = 1000, iter_sampling = 1000,
                                       chains = 4, parallel_chains = 4, save_warmup = T, 
                                       output_dir = "/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/stan_outputs")

saveRDS(bbs_naive_fit, 'bbs_naive_fit.RDS')

dsum <- posterior::summarise_draws(bbs_naive_fit$draws(c("mu_b0", "mu_d0", "sigma_b0", "sigma_d0")))
dsum

posterior::summarise_draws(draws[,!grepl("log_lik", names(draws))])

bbs_naive_loo <- loo(bbs_naive_fit$draws("log_lik"), r_eff=relative_eff(bbs_naive_fit$draws("log_lik")))

saveRDS(bbs_naive_loo, 'bbs_naive_loo.RDS')

##### Class 2 #####
bbs_distance2 <- cmdstan_model(paste0(code_dir, "stan_files/distance_BBS_v2m1.stan"), cpp_options = opencl_options)
# linear distance
dist_data_stan <-list(
    n_visit = 5, 
    n_species = max(fd2$sp_id),
    n_tot = nrow(fd2),
    id_sp = fd2$sp_id,
    Q = fd2$Q,
    dist = fd2$distance_transformed,
    elev = (fd2$elev - mean(fd2$elev))/sd(fd2$elev),
    det_data = rowSums(fd2[,c("v1", "v2", "v3", "v4", "v5")])
)

bbs_dist_fit <- bbs_distance2$sample(data = dist_data_stan, iter_warmup = 1000, iter_sampling = 1000,
                                       chains = 4, parallel_chains = 4, save_warmup = T, 
                                       output_dir = "/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/stan_outputs")
a <- Sys.time()
bdf_draws <- bbs_dist_fit$draws()
b <- Sys.time() - a
test <- read_cmdstan_csv("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/stan_outputs/distance_BBS_v2m1-202104231740-1-6fe406.csv")
saveRDS(bbs_dist_fit, 'bbs_dist_fit.RDS')
bdf_summary <- parsummarise_draws(bbs_dist_fit$draws(), n_cores = 4, n_chunks = 5000)



bbs_lindist_fit2 <- bbs_distance2$sample(data = lindist_data_stan, iter_warmup = 1000, iter_sampling = 1000,
                                       chains = 4, parallel_chains = 4, save_warmup = T, 
                                       output_dir = "/Users/JacobSocolar/Desktop")
saveRDS(bbs_lindist_fit2, 'bbs_lindist_fit2.RDS')


bbs_lindist_fit$summary()
bbs_lindist_fit2$summary()
