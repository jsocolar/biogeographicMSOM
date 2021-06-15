data_package <- readRDS("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_data/wandes_bsd9_package.RDS")
wandes_mod <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/Occupancy/biogeographicMSOM/stan_files/occupancy_v9_wandes.stan",
                         cpp_options = list(stan_threads = TRUE))
wandes_samples <- wandes_mod$sample(data = data_package$data, 
                               chains = 1,
                               parallel_chains = 1,
                               threads_per_chain = 4,
                               refresh = 1,
                               iter_sampling = 1000,
                               iter_warmup = 1000,
                               save_warmup = T,
                               step_size = .0015,
                               max_treedepth = 9,
                               output_dir = "/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_stan_outputs",
                              # inv_metric = bird_stan_data9_1_package$inv_metric,
                               adapt_engaged = T)

wandes_csv <- read_cmdstan_csv("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_stan_outputs/occupancy_v9_wandes-202105052209-1-868b48.csv")

data_package_1 <- data_package
data_package_1$inv_metric <- wandes_csv$inv_metric

saveRDS(data_package_1, "/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_data/wandes_bsd9_package_1.RDS")
#############
wandes_mod <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/Occupancy/biogeographicMSOM/stan_files/occupancy_v9_wandes.stan",
                            cpp_options = list(stan_threads = TRUE))

data_package_1 <- readRDS("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_data/wandes_bsd9_package_1.RDS")
wandes_samples <- wandes_mod$sample(data = data_package_1$data, 
                                    chains = 1,
                                    parallel_chains = 1,
                                    threads_per_chain = 4,
                                    refresh = 1,
                                    iter_sampling = 1000,
                                    iter_warmup = 1000,
                                    save_warmup = T,
                                    step_size = .0015,
                                    max_treedepth = 9,
                                    output_dir = "/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/wandes_stan_outputs",
                                    inv_metric = data_package_1$inv_metric,
                                    adapt_engaged = T)
