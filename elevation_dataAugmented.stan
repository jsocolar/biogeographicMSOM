data {
  // dimensions
  int<lower=1> n_visit; //fixed number of visits
  int<lower=1> n_species; //number of observed species
  int<lower=n_species> S; // superpopulation size
  int<lower=1> n_tot; // number of species-sites (excluding pseudospecies from S)
  int<lower=1> n_totS; // number of species-sites (including psuedospecies from S)
  
  // indexing variables
  int<lower=1> id_sp[n_totS]; // species ID
  int<lower=0, upper=1> Q[n_totS]; // observed detection or not

  // data & covariates
  vector[n_totS] elev;
  vector[n_totS] elev2;
  int<lower=0, upper=n_visit> det_data[n_totS]; // detection history
}

transformed data{
  int<lower=1, upper=n_tot> nq = sum(Q); // number of Q==1
  int det_data1[nq] = segment(det_data, 1, nq);
  int id_sp1[nq] = segment(id_sp, 1, nq);
  int id_sp0[n_tot-nq] = segment(id_sp, nq+1, n_tot-nq);
  int id_spN[n_totS - n_tot] = segment(id_sp, n_tot+1, n_totS - n_tot);
}

parameters {
  real mu_b0; // logit-scale mean intercept for occupancy
  real mu_d0; // logit-scale mean intercept for detection
  
  real<lower=0> sigma_b0;
  vector<offset=mu_b0, multiplier=sigma_b0>[n_species] b0;
  
  real<lower=0> sigma_d0;
  vector<offset=mu_d0, multiplier=sigma_d0>[n_species] d0;
  
  real mu_b1_elev;
  real<lower=0> sigma_b1_elev;
  vector<offset=mu_b1_elev, multiplier=sigma_b1_elev>[n_species] b1_elev;
  
  real mu_b2_elev2;
  real<lower=0> sigma_b2_elev2;
  vector<offset=mu_b2_elev2, multiplier=sigma_b2_elev2>[n_species] b2_elev2;
  
  real<lower=0, upper=1> Omega;
}

model {
  // Transformed parameters that I don't want to monitor
    vector[nq] p_det1 = mu_d0 + d0[id_sp1];    // per-visit logit detection probability at sites with a detection
    vector[n_tot - nq] log_p_no_det0 = log1m_inv_logit(mu_d0 + d0[id_sp0]) * n_visit;  // log probability of no detections given occupancy at sites without a detection
    vector[n_totS - n_tot] log_p_no_detN = log1m_inv_logit(mu_d0 + d0[id_spN]) * n_visit;  // log probability of no detections given occupancy at sites without a detection
    vector[n_totS] logit_psi = mu_b0 + b0[id_sp] + b1_elev[id_sp] .* elev + b2_elev2[id_sp] .* elev2;  // logit occupancy probability
    vector[nq] logit_psi1 = segment(logit_psi, 1, nq);  // logit occupancy probability at sites with a detection
    vector[n_tot-nq] logit_psi0 = segment(logit_psi, nq+1, n_tot-nq); // logit occupancy probability at sites without a detection
    vector[n_totS - n_tot] logit_psiN = segment(logit_psi, n_tot+1, n_totS-n_tot);

  // Priors:
    mu_b0 ~ logistic(0,1);
    mu_d0 ~ logistic(0,1);
    sigma_b0 ~ normal(0,10);
    sigma_d0 ~ normal(0,10);
    
    mu_b1_elev ~ normal(0,20);
    mu_b2_elev2 ~ normal(0,20);
    sigma_b1_elev ~ normal(0,20);
    sigma_b2_elev2 ~ normal(0,20);

  //Random Effects
    b0 ~ normal(mu_b0, sigma_b0);
    d0 ~ normal(mu_d0, sigma_d0);

    b1_elev ~ normal(mu_b1_elev, sigma_b1_elev);
    b2_elev2 ~ normal(mu_b2_elev2, sigma_b2_elev2);

  //Likelihood
    // Observed species
      target += bernoulli_lpmf(1 | Omega)*n_species;
      // Species-Sites with a detection
        target += bernoulli_logit_lpmf(1 | logit_psi1);  // occupancy = 1 for sites with a detection
        target += binomial_logit_lpmf(det_data1 | n_visit, p_det1); // detection = detection history
      // Species-Sites with no detection
        for(i in 1:n_tot-nq){
          target += log_sum_exp(log_inv_logit(logit_psi0[i]) + log_p_no_det0[i], log1m_inv_logit(logit_psi0[i]));
        }
    // Pseudo-species sites
      for(i in 1:n_totS - n_tot){
        target += log_sum_exp(bernoulli_lpmf(0 | Omega), // unavailable
                                bernoulli_lpmf(1 | Omega) + // available
                                log_sum_exp(log_inv_logit(logit_psiN[i]) + log_p_no_detN[i], log1m_inv_logit(logit_psiN[i])) // undetected
                              );
      }
}
