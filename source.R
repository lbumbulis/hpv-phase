
library(dplyr)
library(tidyr)

lambda0 <- 1.7
mu0 <- 0.54
phi0 <- 1

lambda1 <- 0.15

lambda2 <- 1.6
mu2 <- 5   # ideally 8, but reduced for computation reasons
phi2 <- 20

stan_params12 <- data.frame(
  param = c("lambda2", "log_mu2", "log_phi2"),
  value = c(lambda2, log(mu2), log(phi2))
)

stan_params012 <- rbind(
  data.frame(param="lambda1", value=lambda1), stan_params12
)

stan_params012_recur <- data.frame(
  param = paste0("log_", c("lambda1", "mean_s10", "mean_s12", "mu0", "phi0", "mu2", "phi2")),
  value = log(c(lambda1, mu0/lambda0, mu2/lambda2, mu0, phi0, mu2, phi2))
)

# Note on data generation:
# * In the competing gammas representation we generate a new M for each
#   transition of each individual.
#    - Unlike the approach where there is a fixed M common to all transitions
#      of a given type (1->0 or 1->2) across individuals within each dataset,
#      this results in a true gamma mixture distribution (and an identical
#      underlying Weibull approximation) in each dataset.
# * However, with panel data we need to use the Markov phase-type representation
#   of the multistate process, which requires a fixed M at least within each
#   pair of (i, transition type) to force a common dimension on the transition
#   probability matrices being multiplied together.
#    - It may be easier to use a common M across all individuals, but in addition
#      to being conceptually unappealing (approximating different Weibulls on
#      each simulation replicate), this makes it difficult to specify accurate
#      priors for the Negative Binomial hyperparameters on each replicate, as
#      the set of plausible value for these parameters depends on M.
#    - Thus, when we have panel data, we will fix M for each pair of
#      (i, transition type) but allow it to vary across individuals such that
#      it acts as a random effect.

###############################################################################
# Competing gammas representation
###############################################################################
# Generate sojourn time in state 0 before -> 1
s01.fn <- function(N) { rexp(N, rate=lambda1) }

# Generate sojourn time in state 1 before -> 0
s10.fn <- function(N, M=NA) {
  if (anyNA(M)) {
    M <- rnbinom(N, size=phi0, mu=mu0) + 1
  }
  rgamma(N, shape=M, rate=lambda0)
}

# Generate sojourn time in state 1 before -> 2
s12.fn <- function(N, M=NA) {
  if (anyNA(M)) {
    M <- rnbinom(N, size=phi2, mu=mu2) + 1
  }
  rgamma(N, shape=M, rate=lambda2)
}

# Convert progressive process data into natural history dataframe
dfify <- function(s01, s12) {
  n <- length(s01) # should be same as length(s12)
  data.frame(
    id     = rep(1:n, each=2),
    from   = c(rbind(0, s01)),
    to     = c(rbind(s01, s01+s12)),
    from_z = rep(0:1, times=n),
    to_z   = rep(1:2, times=n)
  )
}

# Code in here adapted from Fangya's hpv_precanc.R
###########################################################
# Simulate time to HPV appearance
appear <- function(strt, strt_z){
  n <- length(strt)
  s <- s01.fn(n)
  t_trans <- ifelse(strt_z==2, Inf, strt+s)
  z_trans <- ifelse(strt_z==2, 2, 1)
  return(data.frame(id=seq(n), from=strt, to=t_trans, from_z=strt_z, to_z=z_trans))
}
# Simulate clearance (resolution) of HPV infection
resolve <- function(strt, strt_z, M){
  n <- length(strt)
  clr <- s10.fn(n, M[1])
  prg <- s12.fn(n, M[2])
  t_trans <- ifelse(strt_z==2, Inf, strt + pmin(clr, prg))
  z_trans <- ifelse(strt_z==2, 2, ifelse(clr<prg, 0, 2))
  return(data.frame(id=seq(n), from=strt, to=t_trans, from_z=strt_z, to_z=z_trans))
}

# Create natural history (NH)
generate_NH <- function(n, fix_M=T) {
  if (fix_M) {
    M0 <- rnbinom(1, size=phi0, mu=mu0) + 1
    M2 <- rnbinom(1, size=phi2, mu=mu2) + 1
    M <- c(M0, M2)
  } else {
    M <- c(NA, NA)
  }
  
  out <- data.frame()
  out <- rbind(out, appear(rep(0,n), rep(0,n)))   # 1st appearance
  out <- rbind(out, resolve(out$to, out$to_z, M)) # 1st infection resolved
  
  # Repeat until individuals are too old to be included in studies
  max_iter <- 24
  for (i in 0:max_iter){
    out <- rbind(out, appear(out$to[((2*i+1)*n+1):((2*i+2)*n)], out$to_z[((2*i+1)*n+1):((2*i+2)*n)]))
    out <- rbind(out, resolve(out$to[((2*i+2)*n+1):((2*i+3)*n)], out$to_z[((2*i+2)*n+1):((2*i+3)*n)], M))
  }
  # # Check age and state distribution at end of simulated NH
  # sapply(out[(((max_iter+2)*2-1)*n+1):((max_iter+2)*2*n),], summary)
  
  # out$j <- rep(1:(2*(max_iter+2)), each=n) # integer time index
  
  out <- out[which(out$from_z != 2),]
  out <- out[order(out$id, out$to),]
  
  return(out)
}
###########################################################

# Administratively censor natural history data NH at age tau
# (for modelling based on exact transition times)
admincens <- function(NH, tau) {
  cens.rows <- NH %>%
    group_by(id) %>%
    filter((from < tau) & (to > tau)) %>%
    mutate(
      to = tau,
      to_z = from_z
    ) %>%
    ungroup()
  
  NH <- NH[which(NH$to < tau),]
  out <- as.data.frame(bind_rows(NH, cens.rows) %>%
                         arrange(id, to))
  return(out)
}

my.lag <- function(x, t, x0=NA) {
  x.prev <- c(x0, x[1:(length(x)-1)])
  x.prev[which(t==0)] <- x0
  x.prev
}

# Create panel data from natural history data NH, with visits every
# visit.freq years and administrative censoring at age tau.
panelize <- function(NH, visit.freq, tau) {
  # Add dummy state 2 observations so post-failure states are captured by panelization
  fail.rows <- NH %>%
    group_by(id) %>%
    slice_tail(n=1) %>%
    filter((from_z != 2) & (to_z == 2)) %>%
    mutate(
      from = to,
      to = ceiling(to/visit.freq) * visit.freq + visit.freq/2,
      from_z = to_z,
    ) %>%
    ungroup()
  full.NH <- as.data.frame(bind_rows(NH, fail.rows) %>%
                             arrange(id, to))
  out <- as.data.frame(
    full.NH %>%
      group_by(id) %>%
      reframe(t = seq(0, tau, by=visit.freq)) %>%
      inner_join(full.NH, by="id", relationship="many-to-many") %>%
      filter(t >= from, t < to) %>%
      mutate(to = pmin(to, tau)) %>%
      group_by(id) %>%
      mutate(j = row_number()-1) %>%
      ungroup() %>%
      select(id, j, t, to_z=from_z)
  )
  
  out$from_z <- my.lag(out$to_z, out$j, x0=-1)
  out$t_prev <- my.lag(out$t, out$j)
  out$dt <- out$t - out$t_prev
  
  # Re-order columns
  out <- out[which(!is.na(out$dt)), c(1:2, 6, 3, 7, 5, 4)]
  rownames(out) <- NULL
  
  return(out)
}

# Prep data for Stan under a progressive 0->1->2 process; assumes equal
# times between visits, and assumes state 2 is observed for all subjects
process.panel.012 <- function(pdat) {
  # Aggregate lengths of total observed sojourns in state 0
  nvisits <- pdat %>%
    group_by(id, to_z) %>%
    summarize(n=n(), .groups="drop") %>%
    complete(id, to_z=0:2, fill=list(n=0))
  
  dt <- pdat$t[2] - pdat$t[1]
  
  s0   <- nvisits$n[which(nvisits$to_z==0)] * dt
  idx0 <- which(s0 > 0)
  
  n1.df <- as.data.frame(table(nvisits$n[which(nvisits$to_z==1)]))
  n1.df <- n1.df[which(n1.df$Var1 != 0),]
  
  return(list(
    dt = dt,
    N  = nrow(nvisits) / 3,
    
    N0 = length(idx0),
    s0 = s0[idx0],
    
    nn1       = nrow(n1.df),
    n1_unique = n1.df$Var1,
    n1_count  = n1.df$Freq
  ))
}

# Prep data for Stan under the process with recurrent states 0 and 1
process.panel <- function(pdat) {
  # Label sequence starts; a new 01x sequence begins when from_z=0 and to_z=1
  pdat <- pdat %>%
    arrange(id, t_prev) %>%
    group_by(id) %>%
    mutate(
      seq_start = from_z == 0 & to_z == 1,
      seq_id    = cumsum(seq_start)
    ) %>%
    ungroup()
  
  # Each sequence is the initial (0,1) row plus all following (1,*) rows, up to and
  # including the first row where to_z changes back to 0 or follow-up ends with to_z=1.
  sequences <- pdat %>%
    filter(seq_start | from_z==1) %>%
    group_by(id, seq_id) %>%
    summarise(
      dt_seq      = list(dt),
      t_seq       = list(cumsum(dt)), # cumulative times from sequence start
      final_state = last(to_z),
      .groups     = "drop"
    ) %>%
    mutate(key = sapply(dt_seq, paste, collapse="_"))
  
  # Helper to build flattened arrays from a set of unique sequences
  flatten_seqs <- function(unique_seqs, all_seqs) {
    J_vec     <- sapply(unique_seqs$t_seq, length)
    start_vec <- c(1, cumsum(J_vec[-length(J_vec)])+1)
    t_flat    <- unlist(unique_seqs$t_seq)
    idx       <- all_seqs %>%
      left_join(unique_seqs %>% select(key, k), by="key") %>%
      pull(k)
    list(
      n       = nrow(all_seqs),
      K       = nrow(unique_seqs),
      idx     = idx,
      J       = J_vec,
      t_flat  = t_flat,
      start   = start_vec
    )
  }
  
  # (0, 1, ..., 1, 0) sequences
  seqs_010   <- sequences %>% filter(final_state==0)
  unique_010 <- seqs_010 %>%
    distinct(key, .keep_all = TRUE) %>%
    mutate(k = row_number())
  f010 <- flatten_seqs(unique_010, seqs_010)
  
  # (0, 1, ..., 1) sequences
  seqs_011   <- sequences %>% filter(final_state==1)
  unique_011 <- seqs_011 %>%
    distinct(key, .keep_all = TRUE) %>%
    mutate(k = row_number())
  f011 <- flatten_seqs(unique_011, seqs_011)
  
  # (0, 1, ..., 1, 2) sequences
  seqs_012   <- sequences %>% filter(final_state==2)
  unique_012 <- seqs_012 %>%
    distinct(key, .keep_all = TRUE) %>%
    mutate(k = row_number())
  f012 <- flatten_seqs(unique_012, seqs_012)
  
  # Compute interval lengths from cumulative times for all positions in each flat array
  deltas010 <- numeric(sum(f010$J))
  for (k in seq_len(f010$K)) {
    st      <- f010$start[k]
    J       <- f010$J[k]
    times_k <- f010$t_flat[st:(st + J - 1)]
    deltas010[st:(st + J - 1)] <- diff(c(0, times_k))
  }
  
  deltas011 <- numeric(sum(f011$J))
  for (k in seq_len(f011$K)) {
    st      <- f011$start[k]
    J       <- f011$J[k]
    times_k <- f011$t_flat[st:(st + J - 1)]
    deltas011[st:(st + J - 1)] <- diff(c(0, times_k))
  }
  
  deltas012 <- numeric(sum(f012$J))
  for (k in seq_len(f012$K)) {
    st      <- f012$start[k]
    J       <- f012$J[k]
    times_k <- f012$t_flat[st:(st + J - 1)]
    deltas012[st:(st + J - 1)] <- diff(c(0, times_k))
  }
  
  unique_deltas <- sort(unique(c(deltas010, deltas011, deltas012)))
  
  return(list(
    n010       = f010$n,
    K010       = f010$K,
    idx010     = f010$idx,
    J010       = f010$J,
    t010_flat  = f010$t_flat,
    start010   = f010$start,
    
    n011       = f011$n,
    K011       = f011$K,
    idx011     = f011$idx,
    J011       = f011$J,
    t011_flat  = f011$t_flat,
    start011   = f011$start,
    
    n012       = f012$n,
    K012       = f012$K,
    idx012     = f012$idx,
    J012       = f012$J,
    t012_flat  = f012$t_flat,
    start012   = f012$start,
    
    nd            = length(unique_deltas),
    unique_deltas = unique_deltas,
    delta010_idx  = match(deltas010, unique_deltas),
    delta011_idx  = match(deltas011, unique_deltas),
    delta012_idx  = match(deltas012, unique_deltas)
  ))
}

