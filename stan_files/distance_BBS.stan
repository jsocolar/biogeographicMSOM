functions{
  real partial_sum(int[,] det_slice, 
                   int start, int end, 
                   int n_visit,
                   vector b0, 
                   vector b1,
                   vector d0, 
                   vector d1,
                   vector site_cov1, // some function of distance-to-range
                   row_vector[] vis_cov1, // time (10-stop block)
                   int[] id_sp, 
                   int[] Q) {
    // indexing vars
    int len = end - start + 1;
    int r0 = start - 1;
    
    vector[len] lp;
    real logit_psi;
    row_vector[n_visit] logit_theta;
    for (r in 1:len) {
      // calculate psi & theta
      logit_psi = b0[id_sp[r0+r]] + b1[id_sp[r0+r]]*site_cov1[r0+r];
      logit_theta = d0[id_sp[r0+r]] + d1[id_sp[r0+r]]*vis_cov1[r0+r];
      // likelihood
      if (Q[r0 + r] == 1) 
        lp[r] = log_inv_logit(logit_psi) +
        bernoulli_logit_lpmf(det_slice[r] | logit_theta);
      else lp[r] = log_sum_exp(
        log_inv_logit(logit_psi) +
          log1m_inv_logit(logit_theta[1]) +
          log1m_inv_logit(logit_theta[2]) +
          log1m_inv_logit(logit_theta[3]) +
          log1m_inv_logit(logit_theta[4]),
        log1m_inv_logit(logit_psi));
    } 
    return sum(lp);
  }
}
data {
  // dimensions
  int<lower=1> n_visit; //fixed number of visits
  int<lower=1> n_species; //number of species
  int<lower=1> n_tot; // nrows in df
  
  // indexing variables
  int<lower=1> id_sp[n_tot]; // species ID
  int<lower=0, upper=1> Q[n_tot]; // detection/nondetection by point
  
  // data & covariates
  row_vector[n_visit] vis_cov1[n_tot]; // time (10-stop block)
  vector[n_tot] site_cov1;
  int<lower=0, upper=1> det_data[n_tot, n_visit]; // detection history
  
  // grainsize for reduce_sum 
  int<lower=1> grainsize;
}
parameters {
  real mu_b0;
  real<lower=0> sigma_b0;
  vector[n_species] b0_raw;
  
  real mu_b1;
  real<lower=0> sigma_b1;
  vector[n_species] b1_raw;
  
  real mu_d0;
  real<lower=0> sigma_d0;
  vector[n_species] d0_raw;
  
  real mu_d1;
  real<lower=0> sigma_d1;
  vector[n_species] d1_raw;
}
transformed parameters{
  // occupancy
  vector[n_species] b0 = mu_b0 + b0_raw * sigma_b0;
  vector[n_species] b1 = mu_b1 + b1_raw * sigma_b1;
  
  // detection
  vector[n_species] d0 = mu_d0 + d0_raw * sigma_d0;
  vector[n_species] d1 = mu_d1 + d1_raw * sigma_d1;
}
model {
  //Likelihood
  target += reduce_sum_static(partial_sum, det_data, grainsize, n_visit, b0, b1, d0, d1, site_cov1, vis_cov1, id_sp, Q);
  
  // Hyper-priors:
  mu_b0 ~ normal(0,10);
  mu_b1 ~ normal(0,10);
  mu_d0 ~ normal(0,10);
  mu_d1 ~ normal(0,10);
  
  sigma_b0 ~ normal(0,10);
  sigma_b1 ~ normal(0,10);
  sigma_d0 ~ normal(0,10);
  sigma_d1 ~ normal(0,10);
  
  //Random Effects
  b0_raw ~ normal(0, 1);
  b1_raw ~ normal(0, 1);
  d0_raw ~ normal(0, 1);
  d1_raw ~ normal(0, 1);
}
