functions{
  
  // Compute the log-probability for the species-sites involving observed species at sites with no detections
  real lp_no_det(int n0, vector logit_psi0, vector log_p_no_det0) {
    // n0 number of species-points with no detection
    // logit_psi0 logit psi at those points
    // log_p_no_det0 log probability of no detections at those points
    
    real lp[n0];
    for(i in 1:n0){
      lp[i] = log_sum_exp(bernoulli_logit_lpmf(1 | logit_psi0[i]) + log_p_no_det0[i], bernoulli_logit_lpmf(0 | logit_psi0[i]));
    }
    return sum(lp);
  }

  // Compute the log-probability for the never-observed 
  real lp_no_detN(vector logit_psiN, real log_p_no_det, int n_pt) {
    // At each point, the species is either present but undetected, or not present
    vector[n_pt] lp;
    for(i in 1:n_pt){
      lp[i] = log_sum_exp(bernoulli_logit_lpmf(1 | logit_psiN[i]) + log_p_no_det, bernoulli_logit_lpmf(0 | logit_psiN[i]));
    }
    return sum(lp);
  }
}

data {
  // dimensions
  int<lower=1> n_visit; //fixed number of visits
  int<lower=1> n_species; //number of observed species
  int<lower=1> S; // number of augmented species
  int<lower=1> n_tot; // number of species-sites (excluding pseudospecies from S)

  // indexing variables
  int<lower=1> id_sp[n_tot]; // species ID
  int<lower=0, upper=1> Q[n_tot]; // observed detection or not

  // data & covariates
  vector[n_tot] elev;
  vector[n_tot] elev2;
  int<lower=0, upper=n_visit> det_data[n_tot]; // detection history
  
  int<lower=1> n_pt;
  vector[n_pt] pt_elevs;
  vector[n_pt] pt_elevs2;
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
  vector<offset=mu_b0, multiplier=sigma_b0>[S] b0_S;
  vector<offset=mu_d0, multiplier=sigma_d0>[S] d0_S;
  vector<offset=mu_b1_elev, multiplier=sigma_b1_elev>[S] b1_elev_S;
  vector<offset=mu_b2_elev2, multiplier=sigma_b2_elev2>[S] b2_elev2_S;
}

model {
  // Transformed parameters that I don't want to monitor
    vector[nq] p_det1 = d0[id_sp1];    // per-visit logit detection probability at sites with a detection
    vector[n_tot - nq] log_p_no_det0 = log1m_inv_logit(d0[id_sp0]) * n_visit;  // log probability of no detections given occupancy at sites without a detection
    vector[n_tot] logit_psi = b0[id_sp] + b1_elev[id_sp] .* elev + b2_elev2[id_sp] .* elev2;  // logit occupancy probability
    vector[nq] logit_psi1 = segment(logit_psi, 1, nq);  // logit occupancy probability at sites with a detection
    vector[n_tot-nq] logit_psi0 = segment(logit_psi, nq+1, n_tot-nq); // logit occupancy probability at sites without a detection
    
    vector[n_pt] logit_psiN;
    real logpnodetN;
    
  // Priors:
    mu_b0 ~ logistic(0,1);
    mu_d0 ~ logistic(0,1);
    sigma_b0 ~ normal(0,10);
    sigma_d0 ~ normal(0,10);
    
    mu_b1_elev ~ normal(0,20);
    mu_b2_elev2 ~ normal(0,20);
    sigma_b1_elev ~ normal(0,20);
    sigma_b2_elev2 ~ normal(0,20);
    
    Omega ~ beta(2,2);

  //Random Effects
    b0 ~ normal(mu_b0, sigma_b0);
    d0 ~ normal(mu_d0, sigma_d0);
    b1_elev ~ normal(mu_b1_elev, sigma_b1_elev);
    b2_elev2 ~ normal(mu_b2_elev2, sigma_b2_elev2);
    
    b0_S ~ normal(mu_b0, sigma_b0);
    d0_S ~ normal(mu_d0, sigma_d0);
    b1_elev_S ~ normal(mu_b1_elev, sigma_b1_elev);
    b2_elev2_S ~ normal(mu_b2_elev2, sigma_b2_elev2);

  //Likelihood
    // Observed species
        target += binomial_lpmf(n_species | n_species, Omega);   // All observed species must be available
      // Species-Sites with a detection
        target += bernoulli_logit_lpmf(1 | logit_psi1);  // occupancy = 1 for sites with a detection
        target += binomial_logit_lpmf(det_data1 | n_visit, p_det1); // detection = detection history
      // Species-Sites with no detection
        target += lp_no_det(n_tot-nq, logit_psi0, log_p_no_det0);

    // Pseudo-species
      for(i in 1:S){
        logit_psiN = b0_S[i] + b1_elev_S[i] * pt_elevs + b2_elev2_S[i] * pt_elevs2;
        logpnodetN = binomial_logit_lpmf(0 | n_visit, d0_S[i]);
        target += log_sum_exp( // Either:
                                  // Availabe and never detected on any point
                                    bernoulli_lpmf(1 | Omega) + // available
                                    lp_no_detN(logit_psiN, logpnodetN, n_pt), // never detected on any point
                                // Or
                                  // Unavailable                  
                                    bernoulli_lpmf(0 | Omega)
                              );
      }
}
