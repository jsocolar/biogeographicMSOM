library("cmdstanr"); library("dplyr"); library("posterior")
num_chains <- 1

set.seed(101)
# Define data size ----
n_point <- 20
n_species <- 10
n_visit <- 4

# Define covariates
vis_cov1 <- matrix(data = NA, nrow = n_point, ncol = max(n_visit))
for(i in 1:n_point){
  vis_cov1[i, ] <- runif(n_visit) - .5
}

# Define parameters ----
# Hyperparameters
occ.hyper <- list(b0 = c(0, .5))
b0 <- rnorm(n_species, occ.hyper$b0[1], occ.hyper$b0[2])

det.hyper <- list(d0 = c(-2, .5), d1 = c(0, 1))
d0 <- rnorm(n_species, det.hyper$d0[1], det.hyper$d0[2])
d1 <- rnorm(n_species, det.hyper$d1[1], det.hyper$d1[2])

# Simulate parameters from hyperparameters
logit.occ <- psi <- matrix(NA, nrow = n_point, ncol = n_species)
for(i in 1:n_point){
  for(k in 1:n_species){
    logit.occ[i, k] <- b0[k]
    psi[i, k] <- boot::inv.logit(logit.occ[i, k])
  }
}

logit.det <- theta <- array(NA, dim = c(n_point, max(n_visit), n_species))
for(i in 1:n_point){
  for(j in 1:n_visit){
    for(k in 1:n_species){
      logit.det[i, j, k] <- d0[k] + d1[k]*vis_cov1[i,j]
      theta[i, j, k] <- boot::inv.logit(logit.det[i, j, k])
    }
  }
}

# Simulate data ----
Z <- matrix(NA, nrow = n_point, ncol = n_species)
for(i in 1:n_point){
  for(k in 1:n_species){
    Z[i, k] <- rbinom(1, 1, psi[i, k])
  }
}

det_data <- array(NA, dim = c(n_point, max(n_visit), n_species))
for(i in 1:n_point){
  for(j in 1:n_visit){
    for(k in 1:n_species){
      det_data[i,j,k] <- Z[i, k] * rbinom(1, 1, theta[i,j,k])
    }
  }
}

Q <- apply(det_data, c(1,3), function(x){return(as.numeric(sum(x) > 0))})

stan.data_J <- list(n_point = n_point, n_species = n_species, 
                    n_visit = n_visit,
                    det_data = det_data, 
                    Q = Q, 
                    vis_cov1 = vis_cov1)

# Format for Stan ----
det_df <- apply(det_data, 3, as_data_frame) %>%
  bind_rows(., .id="id_species") %>%
  rename(det_1 = V1, det_2 = V2, det_3 = V3, det_4 = V4) %>%
  mutate(id_species = as.numeric(id_species),
         id_point = rep(1:n_point, n_species), 
         vis_cov1_1 = vis_cov1[id_point, 1], 
         vis_cov1_2 = vis_cov1[id_point, 2], 
         vis_cov1_3 = vis_cov1[id_point, 3],
         vis_cov1_4 = vis_cov1[id_point, 4]) %>%
  mutate(Q = rowSums(select(.,det_1, det_2, det_3, det_4)) > 0) %>%
  arrange(desc(Q), id_species, id_point)

stan_data <- list(
  n_visit = n_visit, 
  n_species = length(unique(det_df$id_species)),
  n_pt = length(unique(det_df$id_point)),
  n_tot = nrow(det_df),
  id_sp = det_df$id_species,
  vis_cov1 = as.matrix(det_df[,paste0("vis_cov1_", 1:n_visit)]),
  det_data = as.matrix(det_df[,paste0("det_", 1:n_visit)]), 
  Q = det_df$Q, 
  grainsize = 200)

num_chains <- 1

setwd('/Users/jacobsocolar/Dropbox/Work/Code/biogeographicMSOM')

# Unthreaded version
mod_R <- cmdstan_model("stan_files/occupancy_problems_reproducible.stan", threads=F)
time_R <- system.time(samps_R <- mod_R$sample(data = stan_data,
                                              num_chains = num_chains,
                                              num_cores = num_chains))

# Threaded version
mod_P <- cmdstan_model("stan_files/occupancy_problems_reproducible.stan", threads=T)

# Run 1 thread
# WORKS
set_num_threads(1)
stan_data$grainsize <- 200
time_1 <- system.time(samps_1 <- mod_P$sample(data = stan_data,
                                              num_chains = num_chains,
                                              num_cores = num_chains))

stan_data$grainsize <- 50
time_1.1 <- system.time(samps_1.1 <- mod_P$sample(data = stan_data,
                                              num_chains = num_chains,
                                              num_cores = num_chains))
# DOES NOT WORK
stan_data$grainsize <- 49
time_1.2 <- system.time(samps_1.2 <- mod_P$sample(data = stan_data,
                                              num_chains = num_chains,
                                              num_cores = num_chains))

# Run 4 thread ----
# WORKS
set_num_threads(4)
stan_data$grainsize <- 50
time_4 <- system.time(samps_4 <- mod_P$sample(data = stan_data,
                                              num_chains = num_chains,
                                              num_cores = num_chains))
# DOES NOT WORK
stan_data$grainsize <- 49
time_4.1 <- system.time(samps_4.1 <- mod_P$sample(data = stan_data,
                                              num_chains = num_chains,
                                              num_cores = num_chains))

# Run 5 thread ----
# WORKS
set_num_threads(5)
stan_data$grainsize <- 50
time_5 <- system.time(samps_5 <- mod_P$sample(data = stan_data,
                                              num_chains = num_chains,
                                              num_cores = num_chains))
# DOES NOT WORK
stan_data$grainsize <- 49
time_5 <- system.time(samps_5 <- mod_P$sample(data = stan_data,
                                              num_chains = num_chains,
                                              num_cores = num_chains))

stan_data$grainsize <- 40
time_5 <- system.time(samps_5 <- mod_P$sample(data = stan_data,
                                              num_chains = num_chains,
                                              num_cores = num_chains))
