data {
  // dimensions
  int<lower=1> n_visit; //fixed number of visits
  int<lower=1> n_species; //number of species
  int<lower=1> n_tot; // number of species-sites
  
  // indexing variables
  int<lower=1> id_sp[n_tot]; // species ID
  int<lower=0, upper=1> Q[n_tot]; // observed detection or not

  // data & covariates
  int<lower=0, upper=n_visit> det_data[n_tot]; // detection history
}

transformed data{
  int<lower=1, upper=n_tot> nq = sum(Q); // number of Q==1
  int det_data1[nq] = segment(det_data, 1, nq);
  int id_sp1[nq] = segment(id_sp, 1, nq);
  int id_sp0[n_tot-nq] = segment(id_sp, nq+1, n_tot-nq);
}

parameters {
  real mu_b0;
  real<lower=0> sigma_b0;
  vector<offset=mu_b0, multiplier=sigma_b0>[n_species] b0;
  
  real mu_d0;
  real<lower=0> sigma_d0;
  vector<offset=mu_d0, multiplier=sigma_d0>[n_species] d0;
}

model {
  // Transformed parameters that I don't want to monitor
    vector[nq] p_det1 = d0[id_sp1];    // per-visit logit detection probability at sites with a detection
    vector[n_tot - nq] log_p_no_det0 = log1m_inv_logit(d0[id_sp0]) * n_visit;  // log probability of no detections given occupancy at sites without a detection
    vector[nq] logit_psi1 = b0[id_sp1];  // logit occupancy probability at sites with a detection
    vector[n_tot-nq] logit_psi0 = b0[id_sp0]; // logit occupancy probability at sites without a detection
  // Priors:
    mu_b0 ~ logistic(0, 1);
    mu_d0 ~ logistic(0, 1);
    sigma_b0 ~ normal(0, 2);
    sigma_d0 ~ normal(0, 2);
  //Random Effects
    b0 ~ normal(mu_b0, sigma_b0);
    d0 ~ normal(mu_d0, sigma_d0);
  //Likelihood
    // Species-Sites with a detection
      target += bernoulli_logit_lpmf(1 | logit_psi1);  // occupancy = 1 for sites with a detection
      target += binomial_logit_lpmf(det_data1 | 5, p_det1); // detection = detection history
    // Species-Sites with no detection
      for(i in 1:n_tot-nq){
        target += log_sum_exp(log_inv_logit(logit_psi0[i]) + log_p_no_det0[i], log1m_inv_logit(logit_psi0[i]));
      }
}

generated quantities {
  vector[n_tot] log_lik;
  for (n in 1:nq) {
    log_lik[n] = bernoulli_logit_lpmf(1 | b0[id_sp[n]]) + binomial_logit_lpmf(det_data[n] | 5, d0[id_sp[n]]);
  }
  for (n in (nq+1):n_tot) {
    log_lik[n] = log_sum_exp(log_inv_logit(b0[id_sp[n]]) + log1m_inv_logit(d0[id_sp[n]]) * n_visit, log1m_inv_logit(b0[id_sp[n]]));
  }
}
