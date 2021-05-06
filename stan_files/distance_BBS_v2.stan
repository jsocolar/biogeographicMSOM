// This contains a severe error (see below).  Keeping the file around so that it's clear what happened.

data {
  // dimensions
  int<lower=1> n_visit; //fixed number of visits
  int<lower=1> n_species; //number of species
  int<lower=1> n_tot; // number of species-sites
  
  // indexing variables
  int<lower=1> id_sp[n_tot]; // species ID
  int<lower=0, upper=1> Q[n_tot]; // observed detection or not

  // data & covariates
  vector[n_tot] elev;
  vector[n_tot] dist;
  int<lower=0, upper=n_visit> det_data[n_tot]; // detection history
}

transformed data{
  int<lower=1, upper=n_tot> nq = sum(Q); // number of Q==1
  int det_data1[nq] = segment(det_data, 1, nq);
  int id_sp1[nq] = segment(id_sp, 1, nq);
  int id_sp0[n_tot-nq] = segment(id_sp, nq+1, n_tot-nq);
}

parameters {
  real<lower=0, upper=1> psi_int; // inverse-logit of the logit-scale mean intercept for occupancy
  real<lower=0, upper=1> theta_int; // inverse-logit of the logit-scale mean intercept for detection
  
  real<lower=-1, upper=1> rho; // correlation between species-specific logit-scale intercepts for occupancy and detection
  
  real<lower=0> sigma_b0;
  vector<offset=logit(psi_int), multiplier=sigma_b0>[n_species] b0;
  
  real<lower=0> sigma_d0;
  vector<offset=logit(theta_int), multiplier=sigma_d0>[n_species] d0;
  
  real mu_b1_elev;
  real<lower=0> sigma_b1_elev;
  vector<offset=mu_b1_elev, multiplier=sigma_b1_elev>[n_species] b1_elev;
  
  real mu_b2_dist;
  real<lower=0> sigma_b2_dist;
  vector<offset=mu_b2_dist, multiplier=sigma_b2_dist>[n_species] b2_dist;
}

transformed parameters {
  real mu_b0 = logit(psi_int);
  real mu_d0 = logit(theta_int);
}

model {
  // Transformed parameters that I don't want to monitor
    vector[n_species] p_det1 = d0[id_sp1];    // per-visit logit detection probability at sites with a detection
    vector[n_species] log_p_no_det0 = log1m_inv_logit(d0[id_sp0]) * n_visit;  // log probability of no detections given occupancy at sites without a detection
    vector[n_species] logit_psi = b0[id_sp] + b1_elev[id_sp] .* elev + b2_dist[id_sp] .* dist;  // logit occupancy probability
    vector[nq] logit_psi1 = segment(logit_psi, 1, nq);  // logit occupancy probability at sites with a detection
    vector[n_tot-nq] logit_psi0 = segment(logit_psi, nq+1, n_tot-nq); // logit occupancy probability at sites without a detection
    
    // Instead of sampling [b0,d0] directly from a multivariate normal, we sample b0 from a normal and then d0 from the 
    // appropriate conditional distribution.
      real conditional_sigma_d0 = sigma_d0*sqrt(1-rho^2); // rho is the covariance divided by the square root of the product of the variances
                                                          // so squaring both sides we have that the conditional variance is 
                                                          // var(d0) - cov(b0,d0)^2/var(b0), which is correct.
      vector[n_species] conditional_mean_d0 = mu_d0 + (rho*sigma_d0/sigma_b0) * (b0 - mu_d0);
                                                          // and here the middle term is cov(b0,d0)/var(b0), which is again correct.
                                                          // THIS IS A SEVERE ERROR.  SHOULD BE b0 - mu_b0!!

  // Priors:
    mu_b1_elev ~ normal(0,20);
    mu_b2_dist ~ normal(0,20);
    
    sigma_b0 ~ normal(0,10);
    sigma_d0 ~ normal(0,10);
    
    sigma_b1_elev ~ normal(0,20);
    sigma_b2_dist ~ normal(0,20);

  //Random Effects
    b0 ~ normal(mu_b0, sigma_b0);
    d0 ~ normal(conditional_mean_d0, conditional_sigma_d0);

    b1_elev ~ normal(mu_b1_elev, sigma_b1_elev);
    b2_dist ~ normal(mu_b2_dist, sigma_b2_dist);

  //Likelihood
    // Species-Sites with a detection
      target += bernoulli_logit_lpmf(1 | logit_psi1);  // occupancy = 1 for sites with a detection
      target += binomial_logit_lpmf(det_data1 | n_visit, p_det1); // detection = detection history
    // Species-Sites with no detection
      for(i in 1:n_tot-nq){
        target += log_sum_exp(log_inv_logit(logit_psi0[i]) + log_p_no_det0[i], log1m_inv_logit(logit_psi0[i]));
      }
}
