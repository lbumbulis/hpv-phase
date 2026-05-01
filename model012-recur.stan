
data {
  // Total number of person-time points ...
  int <lower=0> n01; // ... with 0 -> 1 transition
  int <lower=0> n10; // ... 1 -> 0 transition
  int <lower=0> n12; // ... 1 -> 2 transition
  int <lower=0> n00; // ... stay in state 0 (possible due to admin censoring)
  int <lower=0> n11; // ... stay in state 1
  
  // Sojourn times ...
  vector<lower=0>[n01] s01; // ... in state 0 before 0 -> 1 transition
  vector<lower=0>[n10] s10; // ... in state 1 before 1 -> 0 transition
  vector<lower=0>[n12] s12; // ... in state 1 before 1 -> 2 transition
  vector<lower=0>[n00] s00; // ... stayed in state 0 until end of study
  vector<lower=0>[n11] s11; // ... stayed in state 1 until end of study
  
  int <lower=1> max_M;    // maximum value of M considered for approximating NB
}

parameters {
  // Exponential rate parameter for 0->1
  real <lower=0> lambda1;
  // Gamma rate parameter for 1->0
  real <lower=0> lambda0;
  // Gamma rate parameter for 1->2
  real <lower=0> lambda2;
  
  // Negative Binomial log-mean and log-dispersion parameters for no. of phases M
  // 1->0
  real log_mu0;
  real log_phi0;
  // 1->2
  real log_mu2;
  real log_phi2;
}

transformed parameters {
  real mu0  = exp(log_mu0);
  real phi0 = exp(log_phi0);
  
  real mu2  = exp(log_mu2);
  real phi2 = exp(log_phi2);
}

model {
  // Priors
  lambda0  ~ gamma(1, 1);
  log_mu0  ~ normal(-2, 1);
  log_phi0 ~ normal(1.5, 0.75);
  
  lambda1  ~ gamma(1, 1);
  
  lambda2  ~ gamma(1, 1);
  log_mu2  ~ normal(0, sqrt(2));
  log_phi2 ~ normal(3, 0.5);
  
  // Log-likelihood for 0->1
  s01 ~ exponential(lambda1);
  target += exponential_lccdf(s00 | lambda1);
  
  // Log-likelihood for 1->0 and 1->2 transitions, and 1->1 (admin censored)
  // Note: M0 and M2 are independent, so the joint sum over (M0, M2) factors
  //       into separate terms for each of the competing gamma distributions.
  vector[max_M] log_M0_probs;
  vector[max_M] log_M2_probs;
  for (M in 1:max_M) {
    log_M0_probs[M] = neg_binomial_lpmf(M-1 | phi0, phi0/mu0);
    log_M2_probs[M] = neg_binomial_lpmf(M-1 | phi2, phi2/mu2);
  }
  
  // Density contributions (vectorized)
  vector[max_M] M_vals = linspaced_vector(max_M, 1, max_M);
  matrix[n10, max_M] logf0_C; // conditional on M
  matrix[n12, max_M] logf2_C;
  
  logf0_C = log(s10) * (M_vals-1)'           // outer product for terms depending on s10 and M
              + rep_matrix(-s10*lambda0, max_M)                // terms depending on s10 only
              + rep_matrix(                                    // terms depending on M only
                  (M_vals*log(lambda0) - lgamma(M_vals))', n10
                );
  logf2_C = log(s12) * (M_vals-1)'
              + rep_matrix(-s12*lambda2, max_M)
              + rep_matrix((M_vals*log(lambda2) - lgamma(M_vals))', n12);
  
  matrix[n10, max_M] logf0_M = logf0_C + rep_matrix(log_M0_probs', n10); // marginalize over M
  matrix[n12, max_M] logf2_M = logf2_C + rep_matrix(log_M2_probs', n12);
  
  // Survival function contributions (requires loop)
  // 1->0 transitions
  matrix[n10, max_M] logsurv2_10;
  // 1->2 transitions
  matrix[n12, max_M] logsurv0_12;
  // Admin censored in state 1
  matrix[n11, max_M] logsurv2_11;
  matrix[n11, max_M] logsurv0_11;
  
  // Pre-compute some values so the calculations can be vectorized
  vector[n10] s10_scaled  = s10 * lambda2;
  vector[n12] s12_scaled  = s12 * lambda0;
  vector[n11] s11_scaled2 = s11 * lambda2;
  vector[n11] s11_scaled0 = s11 * lambda0;
  
  for (M in 1:max_M) {
    for (ii in 1:n10) {
      logsurv2_10[ii, M] = gamma_lccdf(s10_scaled[ii] | M, 1) + log_M2_probs[M];
    }
    for (ii in 1:n12) {
      logsurv0_12[ii, M] = gamma_lccdf(s12_scaled[ii] | M, 1) + log_M0_probs[M];
    }
    for (ii in 1:n11) {
      logsurv2_11[ii, M] = gamma_lccdf(s11_scaled2[ii] | M, 1) + log_M2_probs[M];
      logsurv0_11[ii, M] = gamma_lccdf(s11_scaled0[ii] | M, 1) + log_M0_probs[M];
    }
  }
  
  for (ii in 1:n10) {
    target += log_sum_exp(logf0_M[ii]) + log_sum_exp(logsurv2_10[ii]);
  }
  for (ii in 1:n12) {
    target += log_sum_exp(logf2_M[ii]) + log_sum_exp(logsurv0_12[ii]);
  }
  for (ii in 1:n11) {
    target += log_sum_exp(logsurv2_11[ii]) + log_sum_exp(logsurv0_11[ii]);
  }
}

generated quantities {
  real prior_lambda  = gamma_rng(1, 1);
  
  real prior_log_mu0  = normal_rng(-2, 1);
  real prior_log_phi0 = normal_rng(1.5, 0.75);
  
  real prior_log_mu2  = normal_rng(0, sqrt(2));
  real prior_log_phi2 = normal_rng(3, 0.5);
}

