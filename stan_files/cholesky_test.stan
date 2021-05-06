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
  int<lower=1> n_species; //number of species
  vector[2] mus;
}

parameters {
  vector<lower=0>[2] int_sds; // standard deviations of random intercepts
  matrix[2,51] int_standardized; // standardized group-level effects
  cholesky_factor_corr[2] L; // cholesky factor of correlation matrix
}

transformed parameters {
  matrix[n_species, 2] b0_d0;  // actual group-level effects excluding means

  // compute actual group-level effects
  b0_d0 = scale_r_cor(int_standardized, int_sds, L);
}

model {
  // Transformed parameters that I don't want to monitor
  L ~ lkj_corr_cholesky(1);

  // Priors:
  int_sds ~ normal(2,.01);
  target += std_normal_lpdf(to_vector(int_standardized));

}

generated quantities{
  matrix[2,2] vcv = L*L';
}
