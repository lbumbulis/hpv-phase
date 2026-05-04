
functions {
  // Build the transition intensity matrix Q
  matrix Q_fn(int M0, int M2, real lambda0, real lambda1, real lambda2) {
    int n1       = M0 * M2;
    int n_states = n1 + 2;  // state 0, state 1 phase grid, state 2
    matrix[n_states, n_states] Q = rep_matrix(0.0, n_states, n_states);
  
    // State 0 (row 1) -> first phase grid state (col 2)
    Q[1, 2] = lambda1;
  
    // Phase grid transitions (rows 2 to n1+1)
    for (k in 1:n1) {
      int m0 = (k-1) %/% M2 + 1;  // integer division: phase index in M0 direction
      int m2 = (k-1) % M2 + 1;    // remainder:        phase index in M2 direction
      int row = k+1;              // Q row (offset +1 for state 0)
  
      // Advance towards state 0
      if (m0 < M0) {
        Q[row, m0*M2 + m2 + 1] = lambda0; // -> next M0 phase, same M2 phase
      } else {
        Q[row, 1] = lambda0;              // M0 phases exhausted -> back to state 0
      }
      // Advance towards state 2
      if (m2 < M2) {
        Q[row, (m0-1)*M2 + m2 + 2] = lambda2; // -> same M0 phase, next M2 phase
      } else {
        Q[row, n_states] = lambda2;           // M2 phases exhausted -> state 2
      }
    }
  
    // Diagonal: negative row sums
    for (i in 1:n_states) {
      Q[i, i] = -sum(Q[i, ]);
    }
    return Q;
  }
  
  // Log-probability of a sequence of states preceded by state 0, containing
  // n1 1's, and ending in final_state, where the transition probability matrix
  // governing the transitions between consecutive states is given by P.
  real seq_log_prob(matrix P, int n1, int final_state) {
    int n_states = rows(P);
    
    // Calculate using a forward algorithm
    row_vector[n_states] alpha = rep_row_vector(0.0, n_states);
    alpha[1] = 1.0;  // sequences always start in state 0
  
    real log_prob = 0.0;
    
    int tau = (final_state==1) ? n1 : n1+1; // n1 includes the final state if applicable
  
    for (tt in 1:tau) {
      alpha = alpha * P;
  
      real prob_tt;
      if (tt <= n1) {
        prob_tt = sum(alpha[2:(n_states-1)]);
      } else if (final_state == 0) { // && tt == n1+1
        prob_tt = alpha[1];
      } else {// final_state == 2       && tt == n1+1
        prob_tt = alpha[n_states];
      }
      if (prob_tt <= 0) return negative_infinity();
      log_prob += log(prob_tt);
      
      // Zero out probabilities for (possibly latent) states that are
      // incompatible with the state observed at tt
      if (tt <= n1) {
        alpha[1]        = 0.0;
        alpha[n_states] = 0.0;
      } else if (final_state == 0) {
        alpha[2:n_states] = rep_row_vector(0.0, n_states-1);
      } else {
        alpha[1:(n_states-1)] = rep_row_vector(0.0, n_states-1);
      }
      // Normalize to use alpha as initial probs for next tt
      alpha = alpha / prob_tt;
    }
  
    return log_prob;
  }
}

data {
  real<lower=0> dt; // length of time between each pair of consecutive visits
  
  int<lower=0> n00; // total number of (0,0) visit intervals
  int<lower=0> n02; // total number of (0,2) visit intervals
  
  int<lower=1> n_types;                              // number of unique sequence types
  array[n_types] int<lower=1> n1_by_type;            // number of 1's in that sequence
  array[n_types] int<lower=0, upper=2> fs_by_type;   // final state in that sequence
  
  int<lower=1> n_profiles;                           // number of unique sequence profiles
  array[n_profiles] int<lower=1> n_subjects_profile; // number of subjects with that profile
  array[n_profiles] int<lower=1> profile_len;        // number of distinct sequence types in the profile
  
  array[sum(profile_len)] int<lower=1> profile_key_flat;   // all keys for all profiles
  array[sum(profile_len)] int<lower=1> profile_count_flat; // all counts for all profiles
  array[n_profiles] int<lower=1> profile_start;            // maps keys and counts to profiles
  
  int<lower=1> max_M0; // maximum value of M0 considered for approximating NB
  int<lower=1> max_M2; // likewise for M2
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
  log_lambda1 ~ normal(0, 1); // intentionally vague
  
  log_mean_s10 ~ normal(-1.15, 0.5); // log(mu0/lambda0)
  log_mu0      ~ normal(-0.62, 0.25);
  log_phi0     ~ normal(0, 0.25);
  
  log_mean_s12 ~ normal(1.14, 0.5);  // log(mu2/lambda2)
  log_mu2      ~ normal(1.61, 0.5);
  log_phi2     ~ normal(3, 0.5);
  
  // Log-likelihood
  // (0, 1, ..., 1, 2) sequences and (0, 2) intervals
  // Pre-compute distributions of M0 and M2
  vector[max_M0] log_M0_probs;
  vector[max_M2] log_M2_probs;
  for (M0 in 1:max_M0) log_M0_probs[M0] = neg_binomial_lpmf(M0-1 | phi0, phi0/mu0);
  for (M2 in 1:max_M2) log_M2_probs[M2] = neg_binomial_lpmf(M2-1 | phi2, phi2/mu2);
  
  // Condiitonal (on M) log-likelihood contributions for ...
  // ... (0,0) intervals
  matrix[max_M0, max_M2] loglik_c_00;
  // ... (0,2) intervals
  matrix[max_M0, max_M2] loglik_c_02;
  // ... (0,1,...,1,x) sequences
  array[n_types] matrix[max_M0, max_M2] loglik_c_012;
  
  // Pre-compute log-likelihood for each sequence type
  for (M0 in 1:max_M0) {
    for (M2 in 1:max_M2) {
      // Skip pairs with negligible probability, for computational efficiency
      if (log_M0_probs[M0] + log_M2_probs[M2] < -15) {
        loglik_c_00[M0, M2] = negative_infinity();
        loglik_c_02[M0, M2] = negative_infinity();
        for (ii in 1:n_types) loglik_c_012[ii][M0, M2] = negative_infinity();
      } else {
        int n_states = M0*M2 + 2;
        matrix[n_states, n_states] Q = Q_fn(M0, M2, lambda0, lambda1, lambda2);
        matrix[n_states, n_states] P = matrix_exp(Q * dt);
        
        loglik_c_00[M0, M2] = log(P[1, 1]);
        loglik_c_02[M0, M2] = log(P[1, n_states]);
        for (ii in 1:n_types) {
          loglik_c_012[ii][M0, M2] = seq_log_prob(P, n1_by_type[ii], fs_by_type[ii]);
        }
      }
    }
  }
  
  // Marginalize over M for ...
  // ... (0,0) intervals
  target += n00 * log_sum_exp(loglik_c_00 + 
                              rep_matrix(log_M0_probs,  max_M2) +
                              rep_matrix(log_M2_probs', max_M0));
  // ... (0,2) intervals
  target += n02 * log_sum_exp(loglik_c_02 + 
                              rep_matrix(log_M0_probs,  max_M2) +
                              rep_matrix(log_M2_probs', max_M0));
  // ... (0,1,...,1,x) sequences
  for (p in 1:n_profiles) {
    matrix[max_M0, max_M2] loglik_p = rep_matrix(0.0, max_M0, max_M2);
    for (r in profile_start[p]:(profile_start[p] + profile_len[p] - 1)) {
      int key   = profile_key_flat[r];
      int count = profile_count_flat[r];
      loglik_p += count * loglik_c_012[key];
    }
    target += n_subjects_profile[p] * log_sum_exp(loglik_p + 
                                                  rep_matrix(log_M0_probs,  max_M2) +
                                                  rep_matrix(log_M2_probs', max_M0));
  }
}


