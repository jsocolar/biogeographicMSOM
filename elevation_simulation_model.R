# things to try: bimodal species richness, range breadth correlated to maximum occupancy prob, 


library(cmdstanr)
library(posterior)
library(truncnorm)
opencl_options <- list(stan_opencl = TRUE, opencl_platform_id = 0, opencl_device_id = 1)

data_sim <- function(n_sp, n_pt, n_visit = 4, rich_grad = 1, mean_occupancy_core, sd_psi_int, mean_detection, sd_theta_int,
                     mean_lquadratic, sd_lquadratic, pasture = 'n'){
  a <- 2*rich_grad*n_sp/(41 * (1 + rich_grad))
  b <- a/rich_grad
  elev_mids <- (sample(0:40, n_sp, replace = T, prob = seq(a, b, length.out = 41)) + runif(n_sp, -.5, .5))/10 - 2
  mean_psi_int <- boot::logit(mean_occupancy_core)
  psi_int <- rnorm(n_sp, mean_psi_int, sd_psi_int)
  mean_theta_int <- boot::logit(mean_detection)
  theta_int <- rnorm(n_sp, mean_theta_int, sd_theta_int)
  lquadratic <- rtruncnorm(n_sp, b = 0, mean=mean_lquadratic, sd=sd_lquadratic)
  pt_locations <- runif(n_pt, -1, 1)
  if(pasture!='n'){
    pasture_sp <- rnorm(n_sp, pasture[1], pasture[2])
  }
  
  occupancy_probs <- data.frame(elevation = seq(-3, 3, length.out = 600))
  lower <- upper <- rep(NA, n_sp)
  for(i in 1:n_sp){
    sp_psi_logit <-  psi_int[i] + lquadratic[i]*(occupancy_probs$elevation - elev_mids[i])^2
    occupancy_probs <- cbind(occupancy_probs, boot::inv.logit(sp_psi_logit))
    names(occupancy_probs)[i+1] <- paste0("sp_",i)
    lower[i] <- min(occupancy_probs$elevation[occupancy_probs[,i+1] > (.05 * max(occupancy_probs[,i+1]))])
    upper[i] <- max(occupancy_probs$elevation[occupancy_probs[,i+1] > (.05 * max(occupancy_probs[,i+1]))])
  }
  
  simulated_data <- data.frame(point = paste0("point_", 1:n_pt), elev = pt_locations, species = "sp_1", sp_id = 1, pasture = rep(c(1,0), n_pt)[1:n_pt])
  simulated_data$elev2 <- simulated_data$elev^2
  simulated_data$relev_unscaled <- simulated_data$elev - elev_mids[1]
  simulated_data$relev2_unscaled <- simulated_data$relev_unscaled^2
  if(pasture == 'n'){
    simulated_data$psi <- boot::inv.logit(psi_int[1] + lquadratic[1]*simulated_data$relev2_unscaled)
  }else{
    simulated_data$psi <- boot::inv.logit(psi_int[1] + lquadratic[1]*simulated_data$relev2_unscaled + (pasture_sp[1]+pasture[3]*elev_mids[1])*simulated_data$pasture)
  }
  simulated_data$Z <- rbinom(n_pt, 1, simulated_data$psi)
  simulated_data$theta <- boot::inv.logit(theta_int[1])
  simulated_data$n_det <- rbinom(n_pt, n_visit, simulated_data$theta*simulated_data$Z)
  simulated_data$species_in_data <- sum(simulated_data$n_det) > 0
  simulated_data$lower <- lower[1]
  simulated_data$upper <- upper[1]

  for(i in 2:n_sp){
    sim_data <- data.frame(point = paste0("point_", 1:n_pt), elev = pt_locations, species = paste0("sp_", i), sp_id = i, pasture = rep(c(1,0,0), n_pt)[1:n_pt])
    sim_data$elev2 <- sim_data$elev^2
    sim_data$relev_unscaled <- sim_data$elev - elev_mids[i]
    sim_data$relev2_unscaled <- sim_data$relev_unscaled^2
    if(pasture == 'n'){
      sim_data$psi <- boot::inv.logit(psi_int[i] + lquadratic[i]*sim_data$relev2_unscaled)
    }else{
      sim_data$psi <- boot::inv.logit(psi_int[i] + lquadratic[i]*sim_data$relev2_unscaled + (pasture_sp[i]+pasture[3]*elev_mids[i])*sim_data$pasture)
    }
    sim_data$Z <- rbinom(n_pt, 1, sim_data$psi)
    sim_data$theta <- boot::inv.logit(theta_int[i])
    sim_data$n_det <- rbinom(n_pt, n_visit, sim_data$theta*sim_data$Z)
    sim_data$species_in_data <- sum(sim_data$n_det) > 0
    sim_data$lower <- lower[i]
    sim_data$upper <- upper[i]
    simulated_data <- rbind(simulated_data, sim_data)
  }
  return(list(occupancy_probs = occupancy_probs, lower = lower, upper = upper, simulated_data = simulated_data, n_sp = n_sp, pasture_sp = pasture_sp))
}

elev_pasture_model <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/Occupancy/biogeographicMSOM/stan_files/elevation_pasture.stan",
                                    cpp_options = opencl_options)


b_breadth <- n_breadth <- b_center <- n_center <- actual <- b_lower <- b_upper <- n_lower <- n_upper <-list()
b_mse <- n_mse <- vector()
for(j in 1:10){
  set.seed(j+100)
  data1 <- data_sim(n_sp = 100, n_pt = 40, n_visit = 4, rich_grad = 1, mean_occupancy_core = .5, sd_psi_int = 1, mean_detection = .5, sd_theta_int = .5,
                    mean_lquadratic = -50, sd_lquadratic = 10, pasture = c(0,1,0))
  plot(data1$occupancy_probs$sp_1 ~ data1$occupancy_probs$elevation, type = "l",
       ylim = c(0, max(as.matrix(data1$occupancy_probs[,2:(data1$n_sp+1)]))),
       main = "elevational distributions",
       xlab = "elevation",
       ylab = "occupancy probability")
  for(i in 2:data1$n_sp){
    lines(data1$occupancy_probs$elevation, data1$occupancy_probs[,i+1])
  }
  plot(rowSums(as.matrix(data1$occupancy_probs[,2:(data1$n_sp+1)])) ~ data1$occupancy_probs$elevation, xlim = c(-1,1), ylim = c(0, 55), type = "l")
  lines(rowSums(as.matrix(data1$occupancy_probs[,1+unique(data1$simulated_data$sp_id[data1$simulated_data$n_det>0])])) ~ data1$occupancy_probs$elevation)
  pts <- unique(data1$simulated_data$point)
  elevs <- richness <- exp_richness <- pastures <- vector()
  for(i in 1:length(pts)){
    pt <- pts[i]
    richness[i] <- sum(data1$simulated_data$Z[data1$simulated_data$point == pt])
    exp_richness[i] <- sum(data1$simulated_data$psi[data1$simulated_data$point == pt])
    elevs[i] <- unique(data1$simulated_data$elev[data1$simulated_data$point == pt])
    pastures[i] <- unique(data1$simulated_data$pasture[data1$simulated_data$point == pt])
  }
  points(elevs, richness, col = c("black", "red")[pastures+1])
  
  ##########
  
  simulated_data_reduced_prelim <- data1$simulated_data[data1$simulated_data$species_in_data == 1,]
  simulated_data_reduced <- simulated_data_reduced_prelim[order(simulated_data_reduced_prelim$n_det, decreasing = T), ]
  simulated_data_reduced$sp_id_reduced <- as.numeric(as.factor(simulated_data_reduced$species))
  simulated_data_reduced$elevMedian <- (simulated_data_reduced$upper + simulated_data_reduced$lower)/2
  simulated_data_reduced$elevBreadth <- simulated_data_reduced$upper - simulated_data_reduced$lower
  simulated_data_reduced$relev <- 2*(simulated_data_reduced$elev - simulated_data_reduced$elevMedian)/simulated_data_reduced$elevBreadth
  simulated_data_reduced$relev2 <- simulated_data_reduced$relev^2
  n_species <- length(unique(simulated_data_reduced$species))
  
  ##########
  naive_data_reduced <- list(n_visit = 4, n_species = n_species,
                             n_tot = nrow(simulated_data_reduced), id_sp = simulated_data_reduced$sp_id_reduced,
                             Q = as.integer(simulated_data_reduced$n_det > 0), elev = simulated_data_reduced$elev,
                             elev2 = simulated_data_reduced$elev2, pasture = simulated_data_reduced$pasture,
                             det_data = simulated_data_reduced$n_det)
  naive_samples <- elev_pasture_model$sample(data = naive_data_reduced, save_warmup = T, parallel_chains = 4)
  naive_draws <- naive_samples$draws()
  
  bMSOM_data_reduced <- list(n_visit = 4, n_species = n_species,
                             n_tot = nrow(simulated_data_reduced), id_sp = simulated_data_reduced$sp_id_reduced,
                             Q = as.integer(simulated_data_reduced$n_det > 0), elev = simulated_data_reduced$relev,
                             elev2 = simulated_data_reduced$relev2, pasture = simulated_data_reduced$pasture,
                             det_data = simulated_data_reduced$n_det)
  bMSOM_samples <- elev_pasture_model$sample(data = bMSOM_data_reduced, save_warmup = T, parallel_chains = 4)
  bMSOM_draws <- bMSOM_samples$draws()
  
  # bsum <- bMSOM_samples$summary()
  # print(bsum, n=1000)
  b_lower[[j]] <- b_upper[[j]] <- n_lower[[j]] <- n_upper[[j]] <- actual[[j]] <- vector()
  for(i in 1:n_species){
    b_lower[[j]][i] <- quantile(as.numeric(bMSOM_draws[,,paste0("b3_pasture[",i,"]")]), .05)
    b_upper[[j]][i] <- quantile(as.numeric(bMSOM_draws[,,paste0("b3_pasture[",i,"]")]), .95)
    n_lower[[j]][i] <- quantile(as.numeric(naive_draws[,,paste0("b3_pasture[",i,"]")]), .05)
    n_upper[[j]][i] <- quantile(as.numeric(naive_draws[,,paste0("b3_pasture[",i,"]")]), .95)
    actual[[j]][i] <- data1$pasture_sp[unique(simulated_data_reduced$sp_id[simulated_data_reduced$sp_id_reduced == i])]
  }
  
  length(b_upper[[j]])
  b_breadth[[j]] <- b_upper[[j]] - b_lower[[j]]
  b_center[[j]] <- (b_upper[[j]] + b_lower[[j]])/2
  n_breadth[[j]] <- n_upper[[j]] - n_lower[[j]]
  n_center[[j]] <- (n_upper[[j]] + n_lower[[j]])/2
  sum(b_breadth[[j]] > n_breadth[[j]])
  print(sum(actual[[j]] < n_lower[[j]]))
  print(sum(actual[[j]] > n_upper[[j]]))
  print(sum(actual[[j]] < b_lower[[j]]))
  print(sum(actual[[j]] > b_upper[[j]]))
  
  b_mse[j] <- mean((b_center[[j]] - actual[[j]])^2)
  n_mse[j] <- mean((n_center[[j]] - actual[[j]])^2)
  print(b_mse)
  print(n_mse)
  
  plot(b_breadth[[j]] ~ n_breadth[[j]])
  abline(0,1)
}


# bMSOM_samples$summary()

# ##########
# elev_model <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/Occupancy/biogeographicMSOM/stan_files/elevation_traditional.stan",
#                             cpp_options = opencl_options)
# naive_data_reduced <- list(n_visit = 4, n_species = n_species,
#                            n_tot = nrow(simulated_data_reduced), id_sp = simulated_data_reduced$sp_id_reduced,
#                            Q = as.integer(simulated_data_reduced$n_det > 0), elev = simulated_data_reduced$elev,
#                            elev2 = simulated_data_reduced$elev2, det_data = simulated_data_reduced$n_det)
# 
# naive_samples <- elev_model$sample(data = naive_data_reduced, save_warmup = T, parallel_chains = 4)
# 
# bMSOM_data_reduced <- list(n_visit = 4, n_species = n_species,
#                            n_tot = nrow(simulated_data_reduced), id_sp = simulated_data_reduced$sp_id_reduced,
#                            Q = as.integer(simulated_data_reduced$n_det > 0), elev = simulated_data_reduced$relev,
#                            elev2 = simulated_data_reduced$relev2, 
#                            det_data = simulated_data_reduced$n_det)
# 
# bMSOM_samples <- elev_model$sample(data=bMSOM_data_reduced, save_warmup = T, parallel_chains = 4)
# 
# 
# naive_draws <- as_draws_df(naive_samples$draws())
# bMSOM_draws <- as_draws_df(bMSOM_samples$draws())
# 
# occupancy_probs_fitted <- array(data = 0, dim = c(n_species, nrow(data1$occupancy_probs), 4000))
# for(i in 1:n_species){
#   mub0 <- t(as.matrix(naive_draws[, "mu_b0"]))
#   b0i <- t(as.matrix(naive_draws[, paste0("b0[", i, "]")]))
# 
#   logit_occupancy_fitted <- b0i[rep(1, nrow(data1$occupancy_probs)),] +
#         as.matrix(data1$occupancy_probs$elevation) %*% t(as.matrix(naive_draws[, paste0("b1_elev[", i, "]")]))  +
#         data1$occupancy_probs$elevation^2 %*% t(as.matrix(naive_draws[, paste0("b2_elev2[", i, "]")]))
# 
#   occupancy_probs_fitted[i,,] <- boot::inv.logit(logit_occupancy_fitted)
# }
# 
# richnesses <- apply(occupancy_probs_fitted, c(2,3), sum)
# richnesses_q10 <- apply(richnesses, 1, function(x) quantile(x, .1))
# richnesses_q90 <- apply(richnesses, 1, function(x) quantile(x, .9))
# 
# 
# occupancy_probs_fitted_b <- array(data = 0, dim = c(n_species, nrow(data1$occupancy_probs), 4000))
# for(i in 1:n_species){
#   mub0 <- t(as.matrix(bMSOM_draws[, "mu_b0"]))
#   b0i <- t(as.matrix(bMSOM_draws[, paste0("b0[", i, "]")]))
#   lower <- unique(simulated_data_reduced$lower[simulated_data_reduced$sp_id_reduced == i])
#   upper <- unique(simulated_data_reduced$upper[simulated_data_reduced$sp_id_reduced == i])
#   
#   
#   logit_occupancy_fitted <- b0i[rep(1, nrow(data1$occupancy_probs)),] +
#     (2*(as.matrix(data1$occupancy_probs$elevation) - lower)/(upper - lower) - 1) %*% t(as.matrix(bMSOM_draws[, paste0("b1_elev[", i, "]")]))  +
#     (2*(as.matrix(data1$occupancy_probs$elevation) - lower)/(upper - lower) - 1)^2 %*% t(as.matrix(bMSOM_draws[, paste0("b2_elev2[", i, "]")]))
#   
#   occupancy_probs_fitted_b[i,,] <- boot::inv.logit(logit_occupancy_fitted)
# }
# 
# richnesses_b <- apply(occupancy_probs_fitted_b, c(2,3), sum)
# richnesses_q10_b <- apply(richnesses_b, 1, function(x) quantile(x, .1))
# richnesses_q90_b <- apply(richnesses_b, 1, function(x) quantile(x, .9))
# 
# 
# plot(rowSums(as.matrix(data1$occupancy_probs[,2:(data1$n_sp+1)])) ~ data1$occupancy_probs$elevation, xlim = c(-1,1), ylim = c(0,25), type = "n")
# lines(rowSums(as.matrix(data1$occupancy_probs[,1+unique(data1$simulated_data$sp_id[data1$simulated_data$n_det>0])])) ~ data1$occupancy_probs$elevation)
# lines(richnesses_q10 ~ data1$occupancy_probs$elevation, type = "l", col = "red")
# lines(richnesses_q90 ~ data1$occupancy_probs$elevation, type = "l", col = "red")
# lines(richnesses_q10_b ~ data1$occupancy_probs$elevation, type = "l", col = "blue")
# lines(richnesses_q90_b ~ data1$occupancy_probs$elevation, type = "l", col = "blue")
# 
# 
# pts <- unique(data1$simulated_data$point)
# elevs <- richness <- vector()
# for(i in 1:length(pts)){
#   pt <- pts[i]
#   richness[i] <- sum(data1$simulated_data$Z[data1$simulated_data$point == pt])
#   elevs[i] <- unique(data1$simulated_data$elev[data1$simulated_data$point == pt])
# }
# 
# points(elevs, richness)
# 
# 
# # #############
# # augmented_elev_model <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/Occupancy/biogeographicMSOM/stan_files/elevation_dataAugmented.stan",
# #                                       cpp_options = opencl_options)
# # 
# # n_pseudospecies <- 200
# # 
# # naive_data_augmented <- list(n_visit = 4, n_species = n_species, S = n_pseudospecies,
# #                            n_tot = nrow(simulated_data_reduced),
# #                            id_sp = simulated_data_reduced$sp_id_reduced,
# #                            Q = as.integer(simulated_data_reduced$n_det > 0), 
# #                            elev = simulated_data_reduced$elev,
# #                            elev2 = simulated_data_reduced$elev2, 
# #                            det_data = simulated_data_reduced$n_det,
# #                            n_pt = n_pt, pt_elevs = pt_locations, pt_elevs2 = pt_locations^2)
# # 
# # augmented_samples <- augmented_elev_model$sample(data = naive_data_augmented, save_warmup = T, 
# #                                                  chains = 4, parallel_chains = 4,
# #                                                  output_dir = "/Users/JacobSocolar/Desktop")
# # 
# 
# ###############
# augmented_elev_model_looped <- cmdstan_model("/Users/jacobsocolar/Dropbox/Work/Code/Occupancy/biogeographicMSOM/stan_files/data_augmented_elevation_looped.stan",
#                                       cpp_options = opencl_options)
# 
# X <- matrix(data = 0, nrow = n_species, ncol = n_pt)
# uniquept <- unique(simulated_data_reduced$point)[order(unique(simulated_data_reduced$point))]
# for(i in 1:n_species){
#   for(j in 1:n_pt){
#     X[i,j] <- simulated_data_reduced$n_det[simulated_data_reduced$sp_id_reduced == i & simulated_data_reduced$point == uniquept[j]]
#   }
# }
# 
# elevs <- simulated_data_reduced[!duplicated(simulated_data_reduced$point),c("point", "elev", "elev2")]
# elevs_order <- elevs[order(elevs$point), ]
# 
# looped_data <- list(J = n_pt, K = n_visit, n = n_species, 
#                     x = X, S = n_species+200, elev = elevs_order$elev,
#                     elev2 = elevs_order$elev2)
# 
# looped_samples <- augmented_elev_model_looped$sample(data = looped_data, chains = 4, parallel_chains = 4,
#                                                      save_warmup = T, 
#                                                      output_dir = "/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/simulated_elevation")
# looped_draws <- as_draws_df(looped_samples$draws())
# 
# lsumm <- looped_samples$summary()
# print(lsumm, n=200)
# lsumm[lsumm$variable=="Omega",]
# 
# 
# 
# ###########
# occupancy_probs_fitted <- array(data = 0, dim = c(n_species, nrow(data1$occupancy_probs), 4000))
# for(i in 1:n_species){
#   b0i <- t(as.matrix(bMSOM_draws[, paste0("b0[", i, "]")]))
#   
#   logit_occupancy_fitted <-  b0i[rep(1, nrow(data1$occupancy_probs)),] +
#     as.matrix(data1$occupancy_probs$elevation) %*% t(as.matrix(bMSOM_draws[, paste0("b1_elev[", i, "]")]))  +
#     data1$occupancy_probs$elevation^2 %*% t(as.matrix(bMSOM_draws[, paste0("b2_elev2[", i, "]")]))
#   
#   occupancy_probs_fitted[i,,] <- boot::inv.logit(logit_occupancy_fitted)
# }
# 
# richnesses <- apply(occupancy_probs_fitted, c(2,3), sum)
# richnesses_q10 <- apply(richnesses, 1, function(x) quantile(x, .1))
# richnesses_q90 <- apply(richnesses, 1, function(x) quantile(x, .9))
# 
# xl <- 1
# n_sub <- 100
# occupancy_probs_aug_fitted <- array(data = 0, dim = c(200, n_sub, 4000))
# 
# for(i in 1:200){
#   b0i <- t(as.matrix(looped_draws[, paste0("b0[", i+n_species , "]")]))
#   
#   logit_occupancy_fitted <- b0i[rep(1, n_sub),] +
#     as.matrix(seq(-xl,xl,length.out = n_sub)) %*% t(as.matrix(looped_draws[, paste0("b1[", i+n_species , "]")]))  +
#     as.matrix(seq(-xl,xl,length.out = n_sub)^2) %*% t(as.matrix(looped_draws[, paste0("b2[", i+n_species, "]")]))
#   
#   occupancy_probs_aug_fitted[i,,] <- boot::inv.logit(logit_occupancy_fitted)
# }
# 
# richnesses_extra <- apply(occupancy_probs_aug_fitted, c(2,3), sum)
# richnesses_extra_q10 <- apply(richnesses_extra, 1, function(x) quantile(x, .1))
# richnesses_extra_q90 <- apply(richnesses_extra, 1, function(x) quantile(x, .9))
# 
# 
# plot(rowSums(boot::inv.logit(as.matrix(occupancy_probs[,2:(n_sp+1)]))) ~ occupancy_probs$elevation, xlim = c(-xl,xl), ylim = c(0,35), type = "l")
# lines(rowSums(boot::inv.logit(as.matrix(occupancy_probs[,1+unique(simulated_data$sp_id[simulated_data$n_det>0])]))) ~ occupancy_probs$elevation)
# lines(richnesses_q10 ~ data1$occupancy_probs$elevation, type = "l", col = "blue")
# lines(richnesses_q90 ~ data1$occupancy_probs$elevation, type = "l", col = "blue")
# 
# lines(richnesses_extra_q10 ~ seq(-xl,xl,length.out = n_sub), type = "l", col = "blue")
# lines(richnesses_extra_q90 ~ seq(-xl,xl,length.out = n_sub), type = "l", col = "blue")
# 
# pts <- unique(simulated_data$point)
# elevs <- richness <- vector()
# for(i in 1:length(pts)){
#   pt <- pts[i]
#   richness[i] <- sum(simulated_data$Z[simulated_data$point == pt])
#   elevs[i] <- unique(simulated_data$elev[simulated_data$point == pt])
# }
# 
# points(elevs, richness)
# 
# 
# 
# ##################
# 
# 
# occupancy_probs_new <- data.frame(elevation = seq(-4, 4, length.out = 200))
# for(i in 1:n_sp){
#   sp_probs <-  psi_int[i] + lquadratic[i]*(occupancy_probs_new$elevation - elev_mids[i])^2
#   occupancy_probs_new <- cbind(occupancy_probs_new, sp_probs)
#   names(occupancy_probs_new)[i+1] <- paste0("sp_",i)
# }
# 
# plot(boot::inv.logit(occupancy_probs_new$sp_1) ~ occupancy_probs_new$elevation, type = "l",
# #     ylim = c(0, max(boot::inv.logit(as.matrix(occupancy_probs_new[,2:(n_sp+1)])))),
#      main = "elevational distributions",
#      xlab = "elevation",
#      ylab = "occupancy probability",
#      ylim = c(0,.1))
# for(i in 2:n_sp){
#   lines(occupancy_probs_new$elevation, boot::inv.logit(occupancy_probs_new[,i+1]))
# }
# 
# for(i in 1:200){
#   lines(seq(-xl,xl,length.out = n_sub), occupancy_probs_aug_fitted[i,,8])
# }
# 
# 
# 
# 
# 
# 
# 
# unscaled_elev_data_reduced <- list(n_visit = 4, n_species = length(unique(simulated_data_reduced$species)),
#                            n_tot = nrow(simulated_data_reduced), id_sp = simulated_data_reduced$sp_id_reduced,
#                            Q = as.integer(simulated_data_reduced$n_det > 0), elev = simulated_data_reduced$relev_unscaled,
#                            elev2 = simulated_data_reduced$relev2_unscaled, det_data = simulated_data_reduced$n_det)
# 
# relev_unscaled_samples <- naive_model$sample(data = unscaled_elev_data_reduced)
# relev_unscaled_summary <- relev_unscaled_samples$summary()
# print(relev_unscaled_summary, n=1000)
# 
# relev_unscaled_samples$metadata()
# 
# head(simulated_data_reduced)

