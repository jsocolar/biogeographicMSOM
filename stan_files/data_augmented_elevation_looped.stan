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
  vector[J] elev;
  vector[J] elev2;
}

parameters {
  real mu_b0;  //  site-level occupancy
  real<lower=0> sigma_b0;
  vector[S] b0_raw;
  
  real mu_b1;  //  site-level occupancy
  real<lower=0> sigma_b1;
  vector[S] b1_raw;
  
  real<upper=0> mu_b2;  //  site-level occupancy
  real<lower=0> sigma_b2;
  vector[S] b2_raw;
  
  real mu_d0;  //  site-level occupancy
  real<lower=0> sigma_d0;
  vector[S] d0_raw;

  real<lower=0, upper=1> Omega;  // availability of species
}
transformed parameters{
  vector[S] b0 = mu_b0 + b0_raw*sigma_b0;
  vector[S] b1 = mu_b1 + b1_raw*sigma_b1;
  vector[S] b2 = mu_b2 + b2_raw*sigma_b2;
  vector[S] d0 = mu_d0 + b0_raw*sigma_d0;
}

model {
  matrix[S,J] logit_psi;    // log odds  of occurrence
  vector[S] logit_theta;  // log odds of detection
  for (i in 1:S){
    for (j in 1:J){
      logit_psi[i,j] = b0[i] + b1[i]*elev[j] + b2[i]*elev2[j];
    }
    logit_theta[i] = d0[i];
  }
  
  
  // priors
    mu_b0 ~ logistic(0,1);
    mu_d0 ~ logistic(0,1);
    sigma_b0 ~ normal(0,10);
    sigma_d0 ~ normal(0,10);
    
    mu_b1 ~ normal(0,20);
    mu_b2 ~ normal(0,20);
    sigma_b1 ~ normal(0,20);
    sigma_b2 ~ normal(0,20);
    
    b0_raw ~ std_normal();
    d0_raw ~ std_normal();
    b1_raw ~ std_normal();
    b2_raw ~ std_normal();
    
    Omega ~ beta(2,2);

  // likelihood
  for (i in 1:n) {
   1 ~ bernoulli(Omega); // observed, so available
    for (j in 1:J) {
      if (x[i,j] > 0)
        target += lp_observed(x[i,j], K, logit_psi[i,j], logit_theta[i]);
      else
        target += lp_unobserved(K, logit_psi[i,j], logit_theta[i]);
    }
  }
  for (i in (n + 1):S){
    for(j in 1:J){
       target += lp_never_observed(J, K, logit_psi[i,j], logit_theta[i], Omega);
    }
  }
}