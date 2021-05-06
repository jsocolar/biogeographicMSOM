functions {
  /* turn a vector into a matrix of defined dimension 
   * stores elements in row major order
   * Args: 
   *   X: a vector 
   *   N: first dimension of the desired matrix
   *   K: second dimension of the desired matrix 
   * Returns: 
   *   a matrix of dimension N x K 
   */ 
  matrix as_matrix(vector X, int N, int K) { 
    matrix[N, K] Y; 
    for (i in 1:N) {
      Y[i] = to_row_vector(X[((i - 1) * K + 1):(i * K)]); 
    }
    return Y; 
  } 
 /* compute correlated group-level effects
  * Args: 
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns: 
  *   matrix of scaled group-level effects
  */ 
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
}


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
  // Intercepts (logit scale)
    real mu_b0; // occupancy
    real mu_d0; // detection
  // Multivariate random intercept (logit scale)
    vector<lower=0>[2] int_sds; // standard deviations of random intercepts (logit scale)
    matrix[2,n_species] int_standardized; // standardized group-level effects 
    cholesky_factor_corr[2] L; // cholesky factor of correlation matrix
  // Random slopes (logit scale; uncorrelated)
    real mu_b1_elev;
    real<lower=0> sigma_b1_elev;
    vector<offset=mu_b1_elev, multiplier=sigma_b1_elev>[n_species] b1_elev;
  
    real mu_b2_dist;
    real<lower=0> sigma_b2_dist;
    vector<offset=mu_b2_dist, multiplier=sigma_b2_dist>[n_species] b2_dist;
}

transformed parameters {
  // random intercept terms (excluding means)
    matrix[n_species, 2] b0_d0;  
    // using vectors speeds up indexing in loops
      vector[n_species] b0;
      vector[n_species] d0;
    // compute actual group-level effects
      b0_d0 = scale_r_cor(int_standardized, int_sds, L);
      b0 = b0_d0[, 1];
      d0 = b0_d0[, 2];
}

model {
  // Transformed parameters that I don't want to monitor
    vector[n_species] p_det1 = mu_d0 + d0[id_sp1];    // per-visit logit detection probability at sites with a detection
    vector[n_species] log_p_no_det0 = log1m_inv_logit(mu_d0 + d0[id_sp0]) * n_visit;  // log probability of no detections given occupancy at sites without a detection
    vector[n_species] logit_psi = mu_b0 + b0[id_sp] + b1_elev[id_sp] .* elev + b2_dist[id_sp] .* dist;  // logit occupancy probability
    vector[nq] logit_psi1 = segment(logit_psi, 1, nq);  // logit occupancy probability at sites with a detection
    vector[n_tot-nq] logit_psi0 = segment(logit_psi, nq+1, n_tot-nq); // logit occupancy probability at sites without a detection

  // Priors:
    mu_b0 ~ logistic(0,1);
    mu_d0 ~ logistic(0,1);
    int_sds ~ normal(0,10);
    
    mu_b1_elev ~ normal(0,10);
    mu_b2_dist ~ normal(0,10);
    sigma_b1_elev ~ normal(0,10);
    sigma_b2_dist ~ normal(0,10);
    b1_elev ~ normal(mu_b1_elev, sigma_b1_elev);
    b2_dist ~ normal(mu_b2_dist, sigma_b2_dist);
    
    target += std_normal_lpdf(to_vector(int_standardized));
    target += lkj_corr_cholesky_lpdf(L | 1);

  //Likelihood
    // Species-Sites with a detection
      target += bernoulli_logit_lpmf(1 | logit_psi1);  // occupancy = 1 for sites with a detection
      target += binomial_logit_lpmf(det_data1 | n_visit, p_det1); // detection = detection history
    // Species-Sites with no detection
      for(i in 1:n_tot-nq){
        target += log_sum_exp(log_inv_logit(logit_psi0[i]) + log_p_no_det0[i], log1m_inv_logit(logit_psi0[i]));
      }
}




generated quantities {
  vector[n_tot] log_lik;
  for (n in 1:nq) {
    log_lik[n] = bernoulli_logit_lpmf(1 | mu_b0 + b0[id_sp[n]] + b1_elev[id_sp[n]] * elev[n] + b2_dist[id_sp[n]] * dist[n]) + binomial_logit_lpmf(det_data[n] | 5, d0[id_sp[n]]);
  }
  for (n in (nq+1):n_tot) {
    log_lik[n] = log_sum_exp(log_inv_logit(mu_b0 + b0[id_sp[n]] + b1_elev[id_sp[n]] * elev[n] + b2_dist[id_sp[n]] * dist[n]) + log1m_inv_logit(d0[id_sp[n]]) * n_visit, log1m_inv_logit(mu_b0 + b0[id_sp[n]] + b1_elev[id_sp[n]] * elev[n] + b2_dist[id_sp[n]] * dist[n]));
  }
}
