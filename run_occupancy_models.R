# Run full BBS occupancy analysis
library(cmdstanr)
library(posterior)

# # The next three lines only need to be run once.
# source('/Users/jacobsocolar/Dropbox/Work/Code/Occupancy/biogeographicMSOM/data_prep/bbs_import.R')
# source('/Users/jacobsocolar/Dropbox/Work/Code/Occupancy/biogeographicMSOM/data_prep/range_maps.R')
# source('/Users/jacobsocolar/Dropbox/Work/Code/Occupancy/biogeographicMSOM/data_prep/format_for_occupancy.R')

setwd('/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM')
flattened_data <- readRDS('flattened_data.RDS')

# We consider three classes of model:  
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
# The third class (not yet implemented) is not yet well thought through. The idea is to allow model-based
#   estimation of parameters of the function that transforms distance-from-range. For example, we might
#   want to estimate a model whether the linear predictor looks includes a term like
#       b1*(are we out of range) + b2*(distance-from-range)
#   which corresponds to transforming distances to
#       0 [if distance == 0]
#       b1/b2 + distance_from_range [if distance > 0]
#   and then estimating the coefficient b2 from data in the model with a term like
#       b2*transformed_distance
#   Each specific version of this third class will require it's own dedicated Stan model. It's unclear to
#   me whether we'll ultimately want to try none, few, or many possibile functional forms of this sort.

##### Class 1 #####
bbs_bufferclip <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/Occupancy/biogeographicMSOM/stan_files/bufferclip_BBS.stan",
                                threads = T, force_recompile = T)

# 150 km buffer clip
bufferclip_data <- flattened_data[flattened_data$distance_updated < 150000, ]
bufferclip_data_stan <- list(n_species = max(flattened_data$sp_id),
                             n_visit = 4,
                             n_pt = length(unique(flattened_data$site)),
                             n_tot = nrow(bufferclip_data),
                             id_sp = bufferclip_data$sp_id,
                             Q = bufferclip_data$Q,
                             vis_cov1 = matrix(data = rep(c(1,2,4,5), each = nrow(bufferclip_data)), ncol=4),
                             det_data = bufferclip_data[,c(1,2,4,5)],
                             grainsize = 1)

n_cores <- n_chains <- 4
set_num_threads(2)
bbs_bufferclip_fit <- bbs_bufferclip$sample(data = bufferclip_data_stan, num_warmup = 1000, num_samples = 2000,
                                            num_chains = n_chains, num_cores = n_cores)
saveRDS(bbs_bufferclip_fit, 'bbs_bufferclip_fit.RDS')

# naive model (no clip)
naive_data_stan <-list(n_species = max(flattened_data$sp_id),
                       n_visit = 4,
                       n_pt = length(unique(flattened_data$site)),
                       n_tot = nrow(flattened_data),
                       id_sp = flattened_data$sp_id,
                       Q = flattened_data$Q,
                       vis_cov1 = matrix(data = rep(c(1,2,4,5), each = nrow(flattened_data)), ncol=4),
                       det_data = flattened_data[,c(1,2,4,5)],
                       grainsize = 1)

n_cores <- n_chains <- 4
set_num_threads(2)
bbs_naive_fit <- bbs_bufferclip$sample(data = naive_data_stan, num_warmup = 1000, num_samples = 2000,
                                            num_chains = n_chains, num_cores = n_cores)
saveRDS(bbs_naive_fit, 'bbs_naive_fit.RDS')