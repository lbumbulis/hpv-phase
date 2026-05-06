
data {
  int <lower=0> N;        // number of subjects
  vector<lower=0>[N] s01; // sojourn times in state 0 before 0->1
  vector<lower=0>[N] s12; // sojourn times in state 1 before 1->2
  
  int <lower=1> max_M;    // maximum value of M considered for approximating NB
}

parameters {
  // Exponential rate parameter for 0->1
  real <lower=0> lambda1;
  
  // Gamma rate parameter for 1->2
  real <lower=0> lambda2;
  
  // Negative Binomial log-mean and log-dispersion parameters for no. of phases M
  real log_mu2;
  real log_phi2;
}

transformed parameters {
  real mu2  = exp(log_mu2);
  real phi2 = exp(log_phi2);
}

model {
  // Priors
  lambda1  ~ gamma(1, 1);
  lambda2  ~ gamma(1, 1);
  log_mu2  ~ normal(0, sqrt(2));
  log_phi2 ~ normal(3, 0.5);
  
  // Log-likelihood for 0->1
  s01 ~ exponential(lambda1);
  
  // Log-likelihood for 1->2
  vector[max_M] M_vals = linspaced_vector(max_M, 1, max_M);
  
  vector[max_M] log_M_probs;
  for (M in 1:max_M) {
    log_M_probs[M] = neg_binomial_lpmf(M-1 | phi2, phi2/mu2);
  }
  
  matrix[N, max_M] loglik_C; // conditional on M
  loglik_C = log(s12) * (M_vals-1)'                   // outer product for terms depending on s12 and M
              + rep_matrix(-s12*lambda2, max_M)                         // terms depending on s12 only
              + rep_matrix((M_vals*log(lambda2) - lgamma(M_vals))', N); // terms depending on M only
  
  matrix[N, max_M] loglik_M = loglik_C + rep_matrix(log_M_probs', N);   // marginalize over M
  
  for (ii in 1:N) {
    target += log_sum_exp(loglik_M[ii]);
  }
}

generated quantities {
  real prior_lambda  = gamma_rng(1, 1);
  real prior_log_mu2  = normal_rng(0, sqrt(2));
  real prior_log_phi2 = normal_rng(3, 0.5);
}

