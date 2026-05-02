
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
  
  // CURRENTLY UNUSED: Probability of from->to transition with sojourn time s.
  real transition_log_prob(matrix Q, real s, int from, int to) {
    int n_states = rows(Q);
    int n1       = n_states - 2;  // number of phase grid states
  
    matrix[n_states, n_states] P = matrix_exp(Q*s);
  
    int from_idx = (from==0) ? 1 : 2;     // from=1 starts at the first phase
  
    if (to==1) {
      return log(sum(P[from_idx, 2:(n1+1)]));  // sum over all phase grid columns
    } else {
      int to_idx = (to==0) ? 1 : n_states;
      return log(P[from_idx, to_idx]);
    }
  }
  
  // Log-probability of a sequence of states preceded by state 0, where the
  // transition probability matrices governing the transitions between consecutive
  // states are given by P_unique.
  real seq_log_prob(array[] int states, array[] int delta_idx, array[] matrix P_unique) {
    int n_obs    = size(states);  // excludes the known initial state 0
    int n_states = rows(P_unique[1]);
    
    // Calculate using a forward algorithm
    row_vector[n_states] alpha = rep_row_vector(0.0, n_states);
    alpha[1] = 1.0;  // sequences always start in state 0
  
    real log_prob = 0.0;
  
    for (tt in 1:n_obs) {
      alpha = alpha * P_unique[delta_idx[tt]];
  
      real prob_tt;
      if (states[tt] == 0) {
        prob_tt = alpha[1];
      } else if (states[tt] == 1) {
        prob_tt = sum(alpha[2:(n_states-1)]);
      } else {
        prob_tt = alpha[n_states];
      }
      log_prob += log(prob_tt);
      
      // Zero out probabilities for (possibly latent) states that are
      // incompatible with the state observed at tt
      if (states[tt] == 0) {
        alpha[2:n_states] = rep_row_vector(0.0, n_states-1);
      } else if (states[tt] == 1) {
        alpha[1]        = 0.0;
        alpha[n_states] = 0.0;
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
  // (0, 0)
  int<lower=0> n00;         // number of intervals with 0 observed at both visits
  vector<lower=0>[n00] s00; // interval lengths

  // (0, 1, ..., 1, 0)
  // Sequences are stored flattened; start010[k] gives the first index in
  // t010_flat for sequence k, and J010[k] gives its length.
  int<lower=0> n010;
  int<lower=0> K010;                           // number of unique sequences
  array[n010] int<lower=1, upper=K010> idx010; // maps individuals to unique sequences
  array[K010] int<lower=2> J010;               // sequence lengths excluding leading 0 (J >= 2)
  array[sum(J010)] real<lower=0> t010_flat;    // cumulative times from sequence start
  array[K010] int<lower=1> start010;           // start indices into flat arrays

  // (0, 1, ..., 1) at end of follow-up
  int<lower=0> n011;
  int<lower=0> K011;
  array[n011] int<lower=1, upper=K011> idx011;
  array[K011] int<lower=1> J011;               // sequence lengths excluding leading 0 (J >= 1)
  array[sum(J011)] real<lower=0> t011_flat;
  array[K011] int<lower=1> start011;
  
  // (0, 1, ..., 1, 2) sequences
  int<lower=0> n012;
  int<lower=0> K012;
  array[n012] int<lower=1, upper=K012> idx012;
  array[K012] int<lower=2> J012;               // J >= 2: at least one 1 precedes the 2
  array[sum(J012)] real<lower=0> t012_flat;
  array[K012] int<lower=1> start012;
  
  // Unique interval lengths across all sequences, for de-duplicating matrix_exp calls.
  int<lower=1> nd; // number of unique interval lengths ("deltas")
  array[nd] real<lower=0> unique_deltas;
  array[sum(J010)] int<lower=1, upper=nd> delta010_idx; // maps flat positions to unique_deltas
  array[sum(J011)] int<lower=1, upper=nd> delta011_idx;
  array[sum(J012)] int<lower=1, upper=nd> delta012_idx;
  
  int<lower=1> max_M0; // maximum value of M0 considered for approximating NB
  int<lower=1> max_M2; // likewise for M2
}

transformed data {
  int max_M = max(max_M0, max_M2);
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
  log_mu2      ~ normal(0, 1);
  log_phi2     ~ normal(3, 0.5);
  
  // Log-likelihood for (0, 0)
  target += exponential_lccdf(s00 | lambda1);
  
  // Pre-compute joint distribution of (M0,M2)
  vector[max_M] log_M0_probs;
  vector[max_M] log_M2_probs;
  for (M in 1:max_M) {
    log_M0_probs[M] = neg_binomial_lpmf(M-1 | phi0, phi0/mu0);
    log_M2_probs[M] = neg_binomial_lpmf(M-1 | phi2, phi2/mu2);
  }
  
  // Pre-compute log-probabilities of each unique sequence, marginalising over
  // (M0, M2). The outer loop is over (M0, M2) so that the nd P matrices are
  // computed once per (M0, M2) and looked up for every sequence, rather than
  // re-computing matrix_exp for repeated interval lengths.
  // Use braces to enforce local scope and limit memory allocation needs.
  {
    vector[K010] seq010_lp = rep_vector(negative_infinity(), K010);
    vector[K011] seq011_lp = rep_vector(negative_infinity(), K011);
    vector[K012] seq012_lp = rep_vector(negative_infinity(), K012);
    
    for (M0 in 1:max_M0) {
      for (M2 in 1:max_M2) {
        int n_states = M0*M2 + 2;
        matrix[n_states, n_states] Q = Q_fn(M0, M2, lambda0, lambda1, lambda2);
        real lp_weight = log_M0_probs[M0] + log_M2_probs[M2];
        
        // Compute P = matrix_exp(Q * delta) for each unique interval length
        array[nd] matrix[n_states, n_states] P_unique;
        for (d in 1:nd) {
          P_unique[d] = matrix_exp(Q * unique_deltas[d]);
        }
        
        // (0, 1, ..., 1, 0): states passed to seq_log_prob are {1, ..., 1, 0}
        for (k in 1:K010) {
          int J  = J010[k];
          int st = start010[k];
          
          array[J] int states;
          for (j in 1:(J-1)) states[j] = 1;
          states[J] = 0;
          
          seq010_lp[k] = log_sum_exp(seq010_lp[k],
            seq_log_prob(states, delta010_idx[st:(st+J-1)], P_unique) + lp_weight);
        }
        
        // (0, 1, ..., 1): states passed to seq_log_prob are {1, ..., 1}
        for (k in 1:K011) {
          int J  = J011[k];
          int st = start011[k];
          
          array[J] int states;
          for (j in 1:J) states[j] = 1;
          
          seq011_lp[k] = log_sum_exp(seq011_lp[k],
            seq_log_prob(states, delta011_idx[st:(st+J-1)], P_unique) + lp_weight);
        }
        
        // (0, 1, ..., 1, 2): states passed to seq_log_prob are {1, ..., 1, 2}
        for (k in 1:K012) {
          int J  = J012[k];
          int st = start012[k];
          
          array[J] int states;
          for (j in 1:(J-1)) states[j] = 1;
          states[J] = 2;
          
          seq012_lp[k] = log_sum_exp(seq012_lp[k],
            seq_log_prob(states, delta012_idx[st:(st+J-1)], P_unique) + lp_weight);
        }
      }
    }
    
    // Log-likelihood for all non-(0,0) sequences, weighted by number of individuals
    // with each unique sequence via the idx arrays
    target += sum(seq010_lp[idx010]);
    target += sum(seq011_lp[idx011]);
    target += sum(seq012_lp[idx012]);
  }
}

generated quantities {
  real prior_log_lambda1  = normal_rng(0, 1);

  real prior_log_mean_s10 = normal_rng(-1, 1);
  real prior_log_mu0      = normal_rng(-0.5, 1);
  real prior_log_phi0     = normal_rng(0, 0.75);

  real prior_log_mean_s12 = normal_rng(1, 0.5);
  real prior_log_mu2      = normal_rng(0, 1);
  real prior_log_phi2     = normal_rng(3, 0.5);
}

