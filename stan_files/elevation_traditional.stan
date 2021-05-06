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
  vector[n_tot] elev2;
  int<lower=0, upper=n_visit> det_data[n_tot]; // detection history
}

transformed data{
  int<lower=1, upper=n_tot> nq = sum(Q); // number of Q==1
  int det_data1[nq] = segment(det_data, 1, nq);
  int id_sp1[nq] = segment(id_sp, 1, nq);
  int id_sp0[n_tot-nq] = segment(id_sp, nq+1, n_tot-nq);
}

parameters {
  real mu_b0; // logit-scale mean intercept for occupancy
  real mu_d0; // logit-scale mean intercept for detection
  
  real<lower=0> sigma_b0;
  vector[n_species] b0_raw;
  
  real<lower=0> sigma_d0;
  vector[n_species] d0_raw;
  
  real mu_b1_elev;
  real<lower=0> sigma_b1_elev;
  vector[n_species] b1_elev_raw;
  
  real mu_b2_elev2;
  real<lower=0> sigma_b2_elev2;
  vector[n_species] b2_elev2_raw;
}
transformed parameters{
  vector[n_species] b0 = mu_b0 + b0_raw*sigma_b0;
  vector[n_species] d0 = mu_d0 + d0_raw*sigma_d0;
  vector[n_species] b1_elev = mu_b1_elev + b1_elev_raw*sigma_b1_elev;
  vector[n_species] b2_elev2 = mu_b2_elev2 + b2_elev2_raw*sigma_b2_elev2;
}

model {
  // Transformed parameters that I don't want to monitor
    vector[nq] p_det1 = d0[id_sp1];    // per-visit logit detection probability at sites with a detection
    vector[n_tot - nq] log_p_no_det0 = log1m_inv_logit(d0[id_sp0]) * n_visit;  // log probability of no detections given occupancy at sites without a detection
    vector[n_tot] logit_psi = b0[id_sp] + b1_elev[id_sp] .* elev + b2_elev2[id_sp] .* elev2;  // logit occupancy probability
    vector[nq] logit_psi1 = segment(logit_psi, 1, nq);  // logit occupancy probability at sites with a detection
    vector[n_tot-nq] logit_psi0 = segment(logit_psi, nq+1, n_tot-nq); // logit occupancy probability at sites without a detection

  // Priors:
    mu_b0 ~ logistic(0,1);
    mu_d0 ~ logistic(0,1);
    sigma_b0 ~ normal(0,10);
    sigma_d0 ~ normal(0,10);
    
    mu_b1_elev ~ normal(0,200);
    mu_b2_elev2 ~ normal(0,200);
    sigma_b1_elev ~ normal(0,100);
    sigma_b2_elev2 ~ normal(0,100);

  //Random Effects
    b0_raw ~ std_normal();
    d0_raw ~  std_normal();

    b1_elev_raw ~ std_normal();
    b2_elev2_raw ~ std_normal();

  //Likelihood
    // Species-Sites with a detection
      target += bernoulli_logit_lpmf(1 | logit_psi1);  // occupancy = 1 for sites with a detection
      target += binomial_logit_lpmf(det_data1 | n_visit, p_det1); // detection = detection history
    // Species-Sites with no detection
      for(i in 1:n_tot-nq){
        target += log_sum_exp(log_inv_logit(logit_psi0[i]) + log_p_no_det0[i], log1m_inv_logit(logit_psi0[i]));
      }
}
