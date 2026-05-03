
functions {
  // Build the transition intensity matrix Q
  matrix Q_fn_012(int M, real lambda1, real lambda2) {
    int n_states = M+2;  // state 0, M state 1 phases, state 2
    matrix[n_states, n_states] Q = rep_matrix(0.0, n_states, n_states);
  
    // State 0 (row 1) -> first phase (col 2)
    Q[1, 2] = lambda1;
    // Phase transitions (rows 2 to M+1)
    for (m in 1:(M-1)) {
      Q[m+1, m+2] = lambda2;  // advance to next phase
    }
    // Final phase -> state 2
    Q[M+1, n_states] = lambda2;
  
    // Diagonal: negative row sums
    for (i in 1:n_states) {
      Q[i, i] = -sum(Q[i, ]);
    }
    return Q;
  }
  
  // Log-probability of a sequence of states preceded by state 0,
  // where the transition probability matrix governing the transitions
  // between consecutive states is given by P.
  real seq_log_prob(matrix P, int n1) {
    int n_states = rows(P);
    
    // Calculate using a forward algorithm
    row_vector[n_states] alpha = rep_row_vector(0.0, n_states);
    alpha[1] = 1.0;  // sequences always start in state 0
  
    real log_prob = 0.0;
  
    for (tt in 1:(n1+1)) {
      alpha = alpha * P;
  
      real prob_tt;
      if (tt == n1+1) { // -> 2
        prob_tt = alpha[n_states];
      } else {          // -> 1
        prob_tt = sum(alpha[2:(n_states-1)]);
      }
      log_prob += log(prob_tt);
      
      // Zero out probabilities for (possibly latent) states that are
      // incompatible with the state observed at tt
      if (tt == n1+1) { // -> 2
        alpha[1:(n_states-1)] = rep_row_vector(0.0, n_states-1);
      } else {          // -> 1
        alpha[1]        = 0.0;
        alpha[n_states] = 0.0;
      }
      // Normalize to use alpha as initial probs for next tt
      alpha = alpha / prob_tt;
    }
  
    return log_prob;
  }
}

data {
  real<lower=0> dt; // length of time between each pair of consecutive visits
  
  int<lower=1> N;  // number of subjects
  
  int<lower=0> N0;  // number of subjects with at least one visit interval ending in state 0
  array[N0] real<lower=0> s0; // total observed sojourn time in state 0, per subject
  int<lower=0> n02; // number of subjects observed in states 0 and 2 at consecutive visits
  
  int<lower=1> nn1;  // number of unique numbers of visit intervals ending in state 1
  array[nn1] int<lower=1> n1_unique; // array of unique numbers of visit intervals ending in state 1
  array[nn1] int<lower=1> n1_count;  // number individuals with the corresponding n1_unique
  
  int<lower=1> max_M; // maximum value of M2 considered for approximating NB
}

parameters {
  // Exponential log-rate parameter for 0->1
  real log_lambda1;
  
  // Log-mean sojourn time for 1->2 when M2=mu2, i.e. log(mu2/lambda2)
  real log_mean_s12;
  // Negative Binomial log-mean and log-dispersion parameters for no. of state 1 phases M2
  real log_mu2;
  real log_phi2;
}

transformed parameters {
  real lambda1 = exp(log_lambda1);
  real lambda2 = exp(log_mu2 - log_mean_s12);
  real mu2  = exp(log_mu2);
  real phi2 = exp(log_phi2);
}

model {
  // Priors
  log_lambda1  ~ normal(0, 1); // intentionally vague
  log_mean_s12 ~ normal(1.14, 0.5); // log(mu2/lambda2)
  log_mu2      ~ normal(1.61, 0.5);
  log_phi2     ~ normal(3, 0.5);
  
  // Log-likelihood
  // (0, ..., 0) sequences
  target += exponential_lccdf(s0 | lambda1);
  
  // (0, 1, ..., 1, 2) sequences and (0, 2) intervals
  vector[max_M] log_M_probs;
  for (M in 1:max_M) {
    log_M_probs[M] = neg_binomial_lpmf(M-1 | phi2, phi2/mu2);
  }
  
  int n_states;
  
  vector[max_M] loglik_m_02;       // marginal log-likelihood for (0,2) intervals
  matrix[nn1, max_M] loglik_c_012; // conditional (on M) log-likelihood for (0,1,...,2) sequences
  
  for (M in 1:max_M) {
    n_states = M+2;
    matrix[n_states, n_states] Q = Q_fn_012(M, lambda1, lambda2);
    matrix[n_states, n_states] P = matrix_exp(Q * dt);
    
    loglik_m_02[M] = log(P[1, n_states]);
    
    for (ii in 1:nn1) {
      loglik_c_012[ii, M] = seq_log_prob(P, n1_unique[ii]);
    }
  }
  
  target += n02 * log_sum_exp(loglik_m_02 + log_M_probs);
  
  // Marginalize conditional terms over M
  matrix[nn1, max_M] loglik_m_012 = loglik_c_012 + rep_matrix(log_M_probs', nn1);
  
  for (ii in 1:nn1) {
    target += n1_count[ii] * log_sum_exp(loglik_m_012[ii]);
  }
}


