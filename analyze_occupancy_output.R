# analyze output of run_occupancy_models.R

library(cmdstanr)
library(posterior)

setwd('/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM')

bbs_bufferclip_fit <- readRDS('bbs_bufferclip_fit.RDS')
bbs_naive_fit <- readRDS('bbs_naive_fit.RDS')
flattened_data <- readRDS('flattened_data.RDS')
warbler_2018_array <- readRDS('warbler_2018_array.RDS')

##### regenerate stan data objects (from run_occupancy_models.R) #####
bufferclip_data <- flattened_data[flattened_data$distance_updated < 150000, ]
bufferclip_data_stan <- list(n_species = max(flattened_data$sp_id),
                             n_visit = 4,
                             n_pt = length(unique(flattened_data$site)),
                             n_tot = nrow(bufferclip_data),
                             id_sp = bufferclip_data$sp_id,
                             Q = bufferclip_data$Q,
                             vis_cov1 = matrix(data = rep(c(1,2,4,5), each = nrow(bufferclip_data)), ncol=4),
                             det_data = bufferclip_data[,c(1,2,4,5)],    # We withold visit 3 for validation
                             grainsize = 1)

naive_data_stan <-list(n_species = max(flattened_data$sp_id),
                       n_visit = 4,
                       n_pt = length(unique(flattened_data$site)),
                       n_tot = nrow(flattened_data),
                       id_sp = flattened_data$sp_id,
                       Q = flattened_data$Q,
                       vis_cov1 = matrix(data = rep(c(1,2,4,5), each = nrow(flattened_data)), ncol=4),
                       det_data = flattened_data[,c(1,2,4,5)],
                       grainsize = 1)


##### functions get posterior lp for 3rd-visit withhold #####
# bufferclip_Z: get posterior Z matrix for model fit with the bufferclip_BBS.stan code file.  
#     Input: fitted model (cmdstanR output); the data file supplied to cmdstanR
#     Output: posterior Z matrix; flattened version of the withheld detection covariate; and unflattened arrays of logit_psi and logit_theta
bufferclip_Z <- function(bufferclip_stanfit, stan_data){
  draws <- as_draws_df(bufferclip_stanfit$draws())
  logit_psi <- array(data = NA, dim = c(nrow(draws), stan_data$n_species))
  logit_theta <- array(data = NA, dim = c(nrow(draws), stan_data$n_species, 5))
  for(k in 1:stan_data$n_species){
    logit_psi[ , k] <- draws[["mu_b0"]] + draws[["sigma_b0"]] * draws[[paste0("b0_raw[", k, "]")]]
    for(j in 1:5){
      logit_theta[ , k, j] <- draws[["mu_d0"]] + draws[["sigma_d0"]] * draws[[paste0("d0_raw[", k, "]")]] +
        j * (draws[["mu_d1"]] + draws[["sigma_d1"]] * draws[[paste0("d1_raw[", k, "]")]])
    }
  }
  Z <- t3f <- array(data = NA, dim = c(stan_data$n_tot, nrow(draws)))
  for(i in 1:stan_data$n_tot){
    k <- stan_data$id_sp[i]
    if(stan_data$Q[i] == 1){
      Z[i, ] <- 1
    }else{
      psi <- plogis(logit_psi[ , k])
      prod_1m_theta <- (1 - plogis(logit_theta[ , k, 1])) *
        (1 - plogis(logit_theta[ , k, 2])) *
        (1 - plogis(logit_theta[ , k, 4])) *
        (1 - plogis(logit_theta[ , k, 5]))
      Z[i, ] <- (psi * prod_1m_theta)/(psi * prod_1m_theta + (1 - psi))
    }
    t3f[i, ] <- plogis(logit_theta[ , k, 3])
  }
  return(list(Z = Z, theta_v3_flattened = t3f, logit_psi = logit_psi, logit_theta = logit_theta))
}


# bufferclip_lp_v3: get posterior lp for 3rd visit withhold from bufferclip_Z object
#     Input:  output from bufferclip_Z(); logical vector of which rows from flattened_data to include in assessment of posterior lp;
#             logical vector of which rows from the posterior Z-matrix to include in assessment of posterior lp
#       NOTE: flattened_data (see line 9) is not passed as a function argument but must exist in the global environment!!
#     Output: log probability of observed visit-3 data withhold for the desired rows; log probability contibutions of rows with each combination of Q and 3rd visit detection
bufferclip_lp_v3 <- function(bufferclip_Z_object, flattened_include, z_include){
  z_include <- as.logical(z_include)
  flattened_include <- as.logical(flattened_include)
  v3_include <- flattened_data$v3[flattened_include]
  theta3_include <- bufferclip_Z_object$theta_v3_flattened[z_include, ]
  Z_include <- bufferclip_Z_object$Z[z_include, ]
  Z_x_theta3 <- theta3_include*Z_include
  Q_include <- as.numeric(flattened_data$Q[flattened_include])
  v3_array <- matrix(v3_include, length(v3_include), ncol(Z_x_theta3))
  lp_v3 <- log(Z_x_theta3) * v3_array  + log(1 - Z_x_theta3) * (1 - v3_array)
  return(list(all = lp_v3, 
              Q0v0 = lp_v3[Q_include == 0 & v3_include == 0, ],
              Q0v1 = lp_v3[Q_include == 0 & v3_include == 1, ],
              Q1v0 = lp_v3[Q_include == 1 & v3_include == 0, ],
              Q1v1 = lp_v3[Q_include == 1 & v3_include == 1, ]))
}

# Get posterior Z
Z_150km <- bufferclip_Z(bbs_bufferclip_fit, bufferclip_data_stan)
Z_naive <- bufferclip_Z(bbs_naive_fit, naive_data_stan)

# Get posterior log probability of 3rd visit withhold
include_150 <- flattened_data$distance_updated < 150000
lp_v3_150 <- bufferclip_lp_v3(Z_150km, include_150, rep(1, nrow(Z_150km$Z)))
lp_v3_naive <- bufferclip_lp_v3(Z_naive, include_150, include_150)


##### get pointwise posterior alpha-diversity #####
par_alpha <- function(j){
  nr <- nrow(warbler_2018_array$sites)
  breaks <- floor((0:n_cores)*nr/n_cores)
  
  alpha_predicted <- data.frame(site = unique(flattened_data$site)[(breaks[j]+1):breaks[j+1]], clipped_median = NA, clipped_l95 = NA, clipped_u95 = NA,
                                naive_median = NA, naive_l95 = NA, naive_u95 = NA,
                                naive_unclipped_median = NA, naive_unclipped_l95 = NA, naive_unclipped_u95 = NA,
                                clipped_minus_naive_median = NA, clipped_minus_naive_l95 = NA, clipped_minus_naive_u95 = NA,
                                clipped_minus_naive_unclipped_median = NA, clipped_minus_naive_unclipped_l95 = NA, clipped_minus_naive_unclipped_u95 = NA,
                                naive_minus_naive_unclipped_median = NA, naive_minus_naive_unclipped_l95 = NA, naive_minus_naive_unclipped_u95 = NA)
  fd150 <- flattened_data[include_150, ]
  
  get_alpha <- function(x){
    sampled_z <- apply(x, 2, function(y){rbinom(nrow(x), 1, y)})
    if(nrow(x) == 1){
      return(sampled_z)
    }else{
      return(colSums(sampled_z))
    }
  }
  # input posterior Z all species potentially present at a single site (rows are species, columns are posterior iterations)
  # return posterior distribution of alpha diversity. This implementation does not simply sum the probabilities in Z, but rather
  #      it propagates the full binomial sampling uncertainty in the true alpha-diversity of the site.
  
  for(i in 1:(breaks[j+1] - breaks[j])){
    site <- alpha_predicted$site[i]
    
    # get Z for the site from the different models
    site_Z_clipped <- Z_150km$Z[fd150$site == site, , drop = F]
    site_Z_naive <- Z_naive$Z[include_150,][fd150$site == site, , drop = F]
    site_Z_naive_unclipped <- Z_naive$Z[flattened_data$site == site, , drop = F]
    
    # get the alpha diversity posteriors from Z
    site_alpha_clipped <- get_alpha(site_Z_clipped)
    site_alpha_naive <- get_alpha(site_Z_naive)
    site_alpha_naive_unclipped <- get_alpha(site_Z_naive_unclipped)
    
    # Summarize the median and 95% CIs of these posteriors, as well as of the posterior alpha diversity difference
    alpha_predicted$clipped_median[i] <- quantile(site_alpha_clipped, .5)
    alpha_predicted$clipped_l95[i] <- quantile(site_alpha_clipped, .025)
    alpha_predicted$clipped_u95[i] <- quantile(site_alpha_clipped, .975)
    
    alpha_predicted$naive_median[i] <- quantile(site_alpha_naive, .5)
    alpha_predicted$naive_l95[i] <- quantile(site_alpha_naive, .025)
    alpha_predicted$naive_u95[i] <- quantile(site_alpha_naive, .975)
    
    alpha_predicted$naive_unclipped_median[i] <- quantile(site_alpha_naive_unclipped, .5)
    alpha_predicted$naive_unclipped_l95[i] <- quantile(site_alpha_naive_unclipped, .025)
    alpha_predicted$naive_unclipped_u95[i] <- quantile(site_alpha_naive_unclipped, .975)
    
    c_m_n <- site_alpha_clipped - site_alpha_naive
    alpha_predicted$clipped_minus_naive_median[i] <- quantile(c_m_n, .5)
    alpha_predicted$clipped_minus_naive_l95[i] <- quantile(c_m_n, .025)
    alpha_predicted$clipped_minus_naive_u95[i] <- quantile(c_m_n, .975)
    
    c_m_nu <- site_alpha_clipped - site_alpha_naive_unclipped
    alpha_predicted$clipped_minus_naive_unclipped_median[i] <- quantile(c_m_nu, .5)
    alpha_predicted$clipped_minus_naive_unclipped_l95[i] <- quantile(c_m_nu, .025)
    alpha_predicted$clipped_minus_naive_unclipped_u95[i] <- quantile(c_m_nu, .975)
    
    n_m_nu <- site_alpha_naive - site_alpha_naive_unclipped
    alpha_predicted$naive_minus_naive_unclipped_median[i] <- quantile(n_m_nu, .5)
    alpha_predicted$naive_minus_naive_unclipped_l95[i] <- quantile(n_m_nu, .025)
    alpha_predicted$naive_minus_naive_unclipped_u95[i] <- quantile(n_m_nu, .975)
  }
  return(alpha_predicted)
}

n_cores <- 7 # note that n_cores is referenced internally in par_alpha and so must be defined explicitly on its own line, as here

# The below takes tens of minutes to execute, even parallelized across 7 (virtual) cores on my machine.
alpha_predicted_output_list <- parallel::mclapply(1:7, par_alpha, mc.cores = n_cores)

alpha_predicted_output <- do.call(rbind, alpha_predicted_output_list)

saveRDS(alpha_predicted_output, "alpha_predicted_output.RDS")


##### Results summary and plotting #####

# summary of lp by model (colored histograms) and Q and v3 (2x2 grid)
d_lp150 <- density(colSums(lp_v3_150$all), bw = 2)
d_lpnaive <- density(colSums(lp_v3_naive$all), bw = 2)

plot(d_lp150, xlim = c(-9895, -9417), xlab = 'log-probability of observed data', ylab = 'posterior density', main = '')
polygon(d_lp150, col = 'coral')
polygon(d_lpnaive, col = 'blue')

# map of points colored by alpha diversity under the three models, and alpha-diversity difference between the models (pairs plot arrangement?)



# maps of occupancy probabilities for a couple of focal species (big range, small range, small range and nearly undetected)


summary(colSums(lp_v3_150$all))
summary(colSums(lp_v3_naive$all))

summary(colSums(lp_v3_150$Q0v0))
summary(colSums(lp_v3_naive$Q0v0))

summary(colSums(lp_v3_150$Q0v1))
summary(colSums(lp_v3_naive$Q0v1))

summary(colSums(lp_v3_150$Q1v0))
summary(colSums(lp_v3_naive$Q1v0))

summary(colSums(lp_v3_150$Q1v1))
summary(colSums(lp_v3_naive$Q1v1))

dim(lp_v3_150$Q1v0)
all.equal(lp_v3_150$Q1v0[1,], lp_v3_150$Q1v0[4000, ])
lp_v3_150$Q1v0[1:10,1:10]

