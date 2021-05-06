data {
  int<lower=1> N;
  int y[N];
  vector[N] elev;
  vector[N] elev2;
  vector[N] pasture;
  int<lower=1> n_species;
  int<lower=1, upper=n_species> id_sp[N];
}

parameters {
  real mu_b0; // logit-scale mean intercept for occupancy
  
  real<lower=0> sigma_b0;
  vector[n_species] b0_raw;
  
  real mu_b1_elev;
  real<lower=0> sigma_b1_elev;
  vector[n_species] b1_elev_raw;
  
  real mu_b2_elev2;
  real<lower=0> sigma_b2_elev2;
  vector[n_species] b2_elev2_raw;
  
  real mu_b3_pasture;
  real<lower=0> sigma_b3_pasture;
  vector[n_species] b3_pasture_raw;
}

transformed parameters {
  vector[n_species] b0 = mu_b0 + b0_raw*sigma_b0;
  vector[n_species] b1_elev = mu_b1_elev + b1_elev_raw*sigma_b1_elev;
  vector[n_species] b2_elev2 = mu_b2_elev2 + b2_elev2_raw*sigma_b2_elev2;
  vector[n_species] b3_pasture = mu_b3_pasture + b3_pasture_raw*sigma_b3_pasture;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
    vector[N] logit_prob = b0[id_sp] + b1_elev[id_sp] .* elev + b2_elev2[id_sp] .* elev2 + b3_pasture[id_sp] .* pasture;
    target += bernoulli_logit_lpmf(y | logit_prob);
  // Priors:
    mu_b0 ~ logistic(0,1);
    sigma_b0 ~ normal(0,10);
    
    mu_b1_elev ~ normal(0,20);
    mu_b2_elev2 ~ normal(0,20);
    mu_b3_pasture ~ normal(0,20);
    sigma_b1_elev ~ normal(0,20);
    sigma_b2_elev2 ~ normal(0,20);
    sigma_b3_pasture ~ normal(0,20);

  //Random Effects
    b0_raw ~ std_normal();
    b1_elev_raw ~ std_normal();
    b2_elev2_raw ~ std_normal();
    b3_pasture_raw ~ std_normal();
 
}

