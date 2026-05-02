
functions { // for reduce_sum
  real partial_sum_10(array[] real s10_slice, int start, int end,
                      int max_M, vector log_M0_probs, vector log_M2_probs,
                      real lambda0, real lambda2) {
    real lp = 0;
    for (ii in 1:(end - start + 1)) {
      vector[max_M] logf0_row;
      vector[max_M] logsurv2_row;
      for (M in 1:max_M) {
        logf0_row[M]    = gamma_lpdf(s10_slice[ii]          | M, lambda0) + log_M0_probs[M];
        logsurv2_row[M] = gamma_lccdf(s10_slice[ii]*lambda2 | M, 1)       + log_M2_probs[M];
      }
      lp += log_sum_exp(logf0_row) + log_sum_exp(logsurv2_row);
    }
    return lp;
  }

  real partial_sum_12(array[] real s12_slice, int start, int end,
                      int max_M, vector log_M0_probs, vector log_M2_probs,
                      real lambda0, real lambda2) {
    real lp = 0;
    for (ii in 1:(end - start + 1)) {
      vector[max_M] logf2_row;
      vector[max_M] logsurv0_row;
      for (M in 1:max_M) {
        logf2_row[M]    = gamma_lpdf(s12_slice[ii]          | M, lambda2) + log_M2_probs[M];
        logsurv0_row[M] = gamma_lccdf(s12_slice[ii]*lambda0 | M, 1)       + log_M0_probs[M];
      }
      lp += log_sum_exp(logf2_row) + log_sum_exp(logsurv0_row);
    }
    return lp;
  }

  real partial_sum_11(array[] real s11_slice, int start, int end,
                      int max_M, vector log_M0_probs, vector log_M2_probs,
                      real lambda0, real lambda2) {
    real lp = 0;
    for (ii in 1:(end - start + 1)) {
      vector[max_M] logsurv2_row;
      vector[max_M] logsurv0_row;
      for (M in 1:max_M) {
        logsurv2_row[M] = gamma_lccdf(s11_slice[ii]*lambda2 | M, 1) + log_M2_probs[M];
        logsurv0_row[M] = gamma_lccdf(s11_slice[ii]*lambda0 | M, 1) + log_M0_probs[M];
      }
      lp += log_sum_exp(logsurv2_row) + log_sum_exp(logsurv0_row);
    }
    return lp;
  }
}

data {
  // Total number of person-time points ...
  int <lower=0> n01; // ... with 0 -> 1 transition
  int <lower=0> n10; // ... 1 -> 0 transition
  int <lower=0> n12; // ... 1 -> 2 transition
  int <lower=0> n00; // ... stay in state 0 (possible due to admin censoring)
  int <lower=0> n11; // ... stay in state 1
  
  // Sojourn times ...
  vector<lower=0>[n01] s01;     // ... in state 0 before 0 -> 1 transition
  vector<lower=0>[n00] s00;     // ... stayed in state 0 until end of study
  array[n10] real<lower=0> s10; // ... in state 1 before 1 -> 0 transition
  array[n12] real<lower=0> s12; // ... in state 1 before 1 -> 2 transition
  array[n11] real<lower=0> s11; // ... stayed in state 1 until end of study
  
  int <lower=1> max_M;    // maximum value of M considered for approximating NB
  
  int <lower=1> grainsize; // for reduce_sum
}

parameters {
  // Exponential log-rate parameter for 0->1
  real log_lambda1;
  
  // Log-mean sojourn time for 1->0 when M0=mu0, i.e. log(mu0/lambda0)
  real log_mean_s10;
  // Log-mean sojourn time for 1->2 when M2=mu2, i.e. log(mu2/lambda2)
  real log_mean_s12;
  
  // Negative Binomial log-mean and log-dispersion parameters for no. of phases M
  // 1->0
  real log_mu0;
  real log_phi0;
  // 1->2
  real log_mu2;
  real log_phi2;
}

transformed parameters {
  real lambda1 = exp(log_lambda1);
  
  real lambda0 = exp(log_mu0 - log_mean_s10);
  real lambda2 = exp(log_mu2 - log_mean_s12);
  
  real mu0  = exp(log_mu0);
  real phi0 = exp(log_phi0);
  
  real mu2  = exp(log_mu2);
  real phi2 = exp(log_phi2);
}

model {
  // Priors
  log_lambda1 ~ normal(0, 1);
  
  log_mean_s10 ~ normal(-1, 1);
  log_mu0      ~ normal(-0.5, 1);
  log_phi0     ~ normal(0, 0.75);
  
  log_mean_s12 ~ normal(1, 0.5);
  log_mu2      ~ normal(0, sqrt(2));
  log_phi2     ~ normal(3, 0.5);
  
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
  // 1->0
  target += reduce_sum(partial_sum_10, s10, grainsize,
                       max_M, log_M0_probs, log_M2_probs, lambda0, lambda2);
  // 1->2
  target += reduce_sum(partial_sum_12, s12, grainsize,
                       max_M, log_M0_probs, log_M2_probs, lambda0, lambda2);
  // 1->1
  target += reduce_sum(partial_sum_11, s11, grainsize,
                       max_M, log_M0_probs, log_M2_probs, lambda0, lambda2);
}

generated quantities {
  real prior_log_lambda1  = normal_rng(0, 1);

  real prior_log_mean_s10 = normal_rng(-1, 1);
  real prior_log_mu0      = normal_rng(-0.5, 1);
  real prior_log_phi0     = normal_rng(0, 0.75);

  real prior_log_mean_s12 = normal_rng(1, 0.5);
  real prior_log_mu2      = normal_rng(0, sqrt(2));
  real prior_log_phi2     = normal_rng(3, 0.5);
}

