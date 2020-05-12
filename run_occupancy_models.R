# Run full BBS occupancy analysis
library(cmdstanr)

# source('/Users/jacobsocolar/Dropbox/Work/Code/biogeographicMSOM/bbs_import.R')
# source('/Users/jacobsocolar/Dropbox/Work/Code/biogeographicMSOM/range_maps.R')
# source('/Users/jacobsocolar/Dropbox/Work/Code/biogeographicMSOM/format_for_occupancy.R')
# 
# 
# bufferclip_data <- flattened_data[flattened_data$distance_updated < 150000, ]
# n_cores <- n_chains <- 1
# threads_per_chain <- 3
# bufferclip_data_stan <- list(n_species = nrow(species),
#                              n_visit = 4,
#                              n_pt = nrow(sites),
#                              n_tot = nrow(bufferclip_data),
#                              id_sp = bufferclip_data$sp_id,
#                              Q = bufferclip_data$Q,
#                              vis_cov1 = matrix(data = rep(c(1,2,4,5), each = nrow(bufferclip_data)), ncol=4),
#                              det_data = bufferclip_data[,c(1,2,4,5)],
#                              grainsize = ceiling(nrow(bufferclip_data)/threads_per_chain))
# 
# saveRDS(bufferclip_data_stan, file = '/Users/jacobsocolar/Dropbox/Work/biogeographicMSOM/bufferclip_data_stan.RDS')

bufferclip_data_stan <- readRDS('/Users/jacobsocolar/Dropbox/Work/biogeographicMSOM/bufferclip_data_stan.RDS')

bbs_bufferclip <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/biogeographicMSOM/stan_files/bufferclip_BBS.stan",
                                threads = T)
set_num_threads(3)
bbs_bufferclip_fit <- bbs_bufferclip$sample(data = bufferclip_data_stan, num_warmup = 1000, num_samples = 1000,
                                            num_chains = n_chains, num_cores = n_cores)
