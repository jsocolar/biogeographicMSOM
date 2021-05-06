// Adapted from Carpenter 2017 https://mc-stan.org/users/documentation/case-studies/dorazio-royle-occupancy.html
functions {

  real lp_observed(int x, int K, real logit_psi, real logit_theta) {
    return log_inv_logit(logit_psi)
      + binomial_logit_lpmf(x | K, logit_theta);
  }

  real lp_unobserved(int K, real logit_psi, real logit_theta) {
    return log_sum_exp(lp_observed(0, K, logit_psi, logit_theta),
                       log1m_inv_logit(logit_psi));
  }

  real lp_never_observed(int J, int K, real logit_psi, real logit_theta,
                         real Omega) {
      real lp_unavailable = bernoulli_lpmf(0 | Omega);
      real lp_available = bernoulli_lpmf(1 | Omega)
        + J * lp_unobserved(K, logit_psi, logit_theta);
      return log_sum_exp(lp_unavailable, lp_available);
    }
}
data {
  int<lower=1> J;  // sites within region
  int<lower=1> K;  // visits to sites
  int<lower=1> n;  // observed species
  int<lower=0, upper=K> x[n,J];  // visits when species i was detected at site j
  int<lower=n> S;  // superpopulation size
}
parameters {
  real mu_b0; // logit-scale mean intercept for occupancy
  real mu_d0; // logit-scale mean intercept for detection
  
  real<lower=0> sigma_b0;
  vector<offset=mu_b0, multiplier=sigma_b0>[S] b0;
  
  real<lower=0> sigma_d0;
  vector<offset=mu_d0, multiplier=sigma_d0>[S] d0;
  
  real mu_b1_elev;
  real<lower=0> sigma_b1_elev;
  vector<offset=mu_b1_elev, multiplier=sigma_b1_elev>[S] b1_elev;
  
  real mu_b2_elev2;
  real<lower=0> sigma_b2_elev2;
  vector<offset=mu_b2_elev2, multiplier=sigma_b2_elev2>[S] b2_elev2;

  
  real<lower=0, upper=1> Omega;  // availability of species

  real<lower=-1,upper=1> rho_uv;  // correlation of (occupancy, detection)
  vector<lower=0>[2] sigma_uv;    // sd of (occupancy, detection)
  vector[2] uv[S];                // species-level (occupancy, detection)
}
transformed parameters {
  vector[S] logit_psi;    // log odds  of occurrence
  vector[S] logit_theta;  // log odds of detection
  for (i in 1:S)
    logit_psi[i] = b0[i] + b1_elev[i]*elev + b2_relev2[i]*elev2;
  for (i in 1:S)
    logit_theta[i] = d0[i];
}
model {
    vector[n_species] p_det1 = mu_d0 + d0[id_sp1];    // per-visit logit detection probability at sites with a detection
    vector[n_species] log_p_no_det0 = log1m_inv_logit(mu_d0 + d0[id_sp0]) * n_visit;  // log probability of no detections given occupancy at sites without a detection
    vector[n_species] logit_psi = mu_b0 + b0[id_sp] + b1_elev[id_sp] .* elev + b2_elev2[id_sp] .* elev2;  // logit occupancy probability
    vector[nq] logit_psi1 = segment(logit_psi, 1, nq);  // logit occupancy probability at sites with a detection
    vector[n_tot-nq] logit_psi0 = segment(logit_psi, nq+1, n_tot-nq); // logit occupancy probability at sites without a detection

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
  
  
  // priors
  alpha ~ cauchy(0, 2.5);
  beta ~ cauchy(0, 2.5);
  sigma_uv ~ cauchy(0, 2.5);
  (rho_uv + 1) / 2 ~ beta(2, 2);
  uv ~ multi_normal(rep_vector(0, 2), cov_matrix_2d(sigma_uv, rho_uv));
  Omega ~ beta(2,2);
  
  vector[S] logit_psi = mu_b0 + b0[id_sp] + b1_elev[id_sp] .* elev + b2_elev2[id_sp] .* elev2;  // logit occupancy probability


  // likelihood
  for (i in 1:n) {
    1 ~ bernoulli(Omega); // observed, so available
    for (j in 1:J) {
      if (x[i,j] > 0)
        target += lp_observed(x[i,j], K, logit_psi[i], logit_theta[i]);
      else
        target += lp_unobserved(K, logit_psi[i], logit_theta[i]);
    }
  }
  for (i in (n + 1):S)
    target += lp_never_observed(J, K, logit_psi[i], logit_theta[i], Omega);
}
// generated quantities {
//   real<lower=0,upper=S> E_N = S * Omega; // model-based expectation species
//   int<lower=0,upper=S> E_N_2;  // posterior simulated species
//   vector[2] sim_uv;
//   real logit_psi_sim;
//   real logit_theta_sim;
// 
//   E_N_2 = n;
//   for (i in (n+1):S) {
//     real lp_unavailable = bernoulli_lpmf(0 | Omega);
//     real lp_available = bernoulli_lpmf(1 | Omega)
//       + J * lp_unobserved(K, logit_psi[i], logit_theta[i]);
//     real Pr_available = exp(lp_available
//                             - log_sum_exp(lp_unavailable, lp_available));
//     E_N_2 = E_N_2 + bernoulli_rng(Pr_available);
//   }
// 
//   sim_uv = multi_normal_rng(rep_vector(0,2),
//                              cov_matrix_2d(sigma_uv, rho_uv));
//   logit_psi_sim = alpha + sim_uv[1];
//   logit_theta_sim = beta + sim_uv[2];
// }