functions{
    real partial_sum(int[,] det_slice, 
                     int start, int end, 
                     vector b0, 
                     vector d0, 
                     vector d5,
                     row_vector[] vis_cov1,
                     int[] id_sp, 
                     int[] Q) {
        // indexing vars
        int len = end - start + 1;
        int r0 = start - 1;
        
        vector[len] lp;
        vector[len] logit_psi;
        row_vector[4] logit_theta[len];
        for (r in 1:len) {
            // calculate psi & theta
            logit_psi[r] = b0[id_sp[r0+r]];
            logit_theta[r] = d0[id_sp[r0+r]] + d5[id_sp[r0+r]]*vis_cov1[r0+r];
            // likelihood
            if (Q[r0 + r] == 1) 
                lp[r] = log_inv_logit(logit_psi[r]) +
                    bernoulli_logit_lpmf(det_slice[r0 + r] | logit_theta[r]);
            else lp[r] = log_sum_exp(
                log_inv_logit(logit_psi[r]) +
                    log1m_inv_logit(logit_theta[r, 1]) +
                    log1m_inv_logit(logit_theta[r, 2]) +
                    log1m_inv_logit(logit_theta[r, 3]) +
                    log1m_inv_logit(logit_theta[r, 4]),
                log1m_inv_logit(logit_psi[r]));
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
    int<lower=1> id_sp[n_tot]; 
    int<lower=0, upper=1> Q[n_tot]; // species:cluster
    
    // data & covariates
    row_vector[n_visit] vis_cov1[n_tot]; // visit covariate 1
    int<lower=0, upper=1> det_data[n_tot, n_visit]; // detection history
    
    // grainsize for reduce_sum 
    int<lower=1> grainsize;
}
parameters {
    real mu_b0;
    real<lower=0> sigma_b0;
    vector[n_species] b0_raw;
    
    real mu_d0;
    real<lower=0> sigma_d0;
    vector[n_species] d0_raw;
    
    real mu_d5;
    real<lower=0> sigma_d5;
    vector[n_species] d5_raw;
}
transformed parameters{
    // occupancy
    vector[n_species] b0 = mu_b0 + b0_raw * sigma_b0;
    
    // detection
    vector[n_species] d0 = mu_d0 + d0_raw * sigma_d0;
    vector[n_species] d5 = mu_d5 + d5_raw * sigma_d5;
}
model {
    //Likelihood
    target += reduce_sum_static(partial_sum, det_data, grainsize, b0, d0, d5, vis_cov1, id_sp, Q);

    // Hyper-priors:
    mu_b0 ~ normal(0,10);
    mu_d0 ~ normal(0,10);
    mu_d5 ~ normal(0,10);
    
    sigma_b0 ~ normal(0,10);
    sigma_d0 ~ normal(0,10);
    sigma_d5 ~ normal(0,10);
    
    //Random Effects
    b0_raw ~ normal(0, 1);
    d0_raw ~ normal(0, 1);
    d5_raw ~ normal(0, 1);
}
