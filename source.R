
library(dplyr)
library(tidyr)

lambda0 <- 1.7
mu0 <- 0.54
phi0 <- 1

lambda1 <- 0.15

lambda2 <- 1.6
mu2 <- 5   # ideally 8, but reduced for computation reasons
phi2 <- 20

stan_params012 <- data.frame(
  param = paste0("log_", c("lambda1", "mean_s12", "mu2", "phi2")),
  value = log(c(lambda1, mu2/lambda2, mu2, phi2))
)

stan_params <- data.frame(
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

# Generate sojourn time in state 0 before -> 1
s01.fn <- function(N) { rexp(N, rate=lambda1) }

# Generate sojourn time in state 1 before -> 0
s10.fn <- function(N, M=NULL) {
  if (is.null(M)) {
    M <- rnbinom(N, size=phi0, mu=mu0) + 1
  }
  rgamma(N, shape=M, rate=lambda0)
}

# Generate sojourn time in state 1 before -> 2
s12.fn <- function(N, M=NULL) {
  if (is.null(M)) {
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
  clr <- s10.fn(n, M[[1]])
  prg <- s12.fn(n, M[[2]])
  t_trans <- ifelse(strt_z==2, Inf, strt + pmin(clr, prg))
  z_trans <- ifelse(strt_z==2, 2, ifelse(clr<prg, 0, 2))
  return(data.frame(id=seq(n), from=strt, to=t_trans, from_z=strt_z, to_z=z_trans))
}

# Create natural history (NH)
generate_NH <- function(n, M=list(NULL, NULL)) {
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
process.panel.012 <- function(pdat, visit.freq) {
  # Aggregate lengths of total observed sojourns in state 0
  nvisits <- pdat %>%
    group_by(id, to_z) %>%
    summarize(n=n(), .groups="drop") %>%
    complete(id, to_z=0:2, fill=list(n=0))
  
  s0   <- nvisits$n[which(nvisits$to_z==0)] * visit.freq
  idx0 <- which(s0 > 0)
  
  # Subjects with n1=0 but n2>0: observed jumping from state 0 to state 2
  n02 <- sum(nvisits$n[which(nvisits$to_z==1)] == 0 & 
               nvisits$n[which(nvisits$to_z==2)] > 0)
  
  # Count the number of visits in state 1 for each individual
  n1.df <- as.data.frame(table(nvisits$n[which(nvisits$to_z==1)]))
  n1.df$Var1 <- as.integer(as.character(n1.df$Var1)) # convert from factor
  n1.df <- n1.df[which(n1.df$Var1 != 0),]
  
  return(list(
    dt = visit.freq,
    N  = nrow(nvisits) / 3,
    
    N0  = length(idx0),
    s0  = s0[idx0],
    n02 = n02,
    
    nn1       = nrow(n1.df),
    n1_unique = n1.df$Var1,
    n1_count  = n1.df$Freq
  ))
}

# Prep data for Stan under the process with recurrent states 0 and 1
process.panel <- function(pdat, visit.freq) {
  # Per-subject (0,0) and (0,2) interval counts
  subject_interval_counts <- pdat %>%
    filter(from_z==0, to_z %in% c(0,2)) %>%
    count(id, to_z) %>%
    pivot_wider(names_from=to_z, values_from=n, names_prefix="n0", values_fill=0)
  
  # Identify all (0,1,...,1,x) sequences
  pdat <- pdat %>%
    arrange(id, t) %>%
    group_by(id) %>%
    mutate(seq_start=from_z==0 & to_z==1, seq_id=cumsum(seq_start)) %>%
    ungroup()
  
  sequences <- pdat %>%
    filter(seq_start | from_z==1) %>%   # exclude (0,0) and (0,2) rows
    group_by(id, seq_id) %>%
    summarise(n1=sum(to_z==1), final_state=last(to_z), .groups="drop")
  
  # Key describing the properties of each type of sequence
  seq_key <- aggregate(seq_id ~ n1 + final_state, FUN=length, data=sequences)
  seq_key <- cbind(key=1:nrow(seq_key), seq_key)
  names(seq_key)[ncol(seq_key)] <- "count"
  
  # Merge the key labels back into the dataframe of sequences
  sequences <- merge(sequences, seq_key[,1:3], sort=F)
  # Count the number of each sequence type for each individual
  seq_agg <- aggregate(seq_id ~ id + key, FUN=length, data=sequences)
  names(seq_agg)[ncol(seq_agg)] <- "count"
  
  # Build per-subject profile strings, including (0,0) and (0,2) intervals
  seq_summary <- seq_agg %>%
    arrange(id, key) %>%   # sort by key so order of sequences doesn't matter
    group_by(id) %>%
    summarise(seq_str=paste(paste(key, count, sep="x"), collapse="_"), .groups="drop")
  
  subject_profiles <- data.frame(id=unique(pdat$id)) %>%
    left_join(subject_interval_counts, by="id") %>%
    replace_na(list(n00=0, n02=0)) %>%
    left_join(seq_summary, by="id") %>%
    replace_na(list(seq_str="")) %>%
    mutate(profile=paste(seq_str, n00, n02, sep="|")) %>%
    group_by(profile) %>%
    mutate(profile_id=cur_group_id(), n_subjects=n()) %>%
    ungroup()
  
  # One representative per profile
  representative_ids <- subject_profiles %>%
    group_by(profile_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(id, profile_id, n_subjects, n00, n02)
  
  # Flatten profile-level sequence aggregates into ragged arrays for Stan
  flat <- representative_ids %>%
    left_join(
      seq_agg %>%
        inner_join(representative_ids %>% select(id, profile_id), by="id") %>%
        arrange(profile_id, key) %>%
        group_by(profile_id) %>%
        summarise(keys=list(key), counts=list(count), n_seqtypes=n(), .groups="drop"),
      by="profile_id"
    ) %>%
    replace_na(list(n_seqtypes=0))
  
  # Replace NULL list entries (profiles with no sequences) with empty integer vectors
  flat$keys   <- lapply(flat$keys,   function(x) if (is.null(x)) integer(0) else x)
  flat$counts <- lapply(flat$counts, function(x) if (is.null(x)) integer(0) else x)
  
  return(list(
    dt = visit.freq,
    
    n_types    = nrow(seq_key),
    n1_by_type = seq_key$n1,
    fs_by_type = seq_key$final_state,
    
    n_profiles         = nrow(flat),
    n_subjects_profile = flat$n_subjects,
    profile_n00        = flat$n00,
    profile_n02        = flat$n02,
    profile_len        = flat$n_seqtypes,
    profile_start      = c(1, cumsum(flat$n_seqtypes[-nrow(flat)]) + 1),
    profile_key_flat   = unlist(flat$keys),
    profile_count_flat = unlist(flat$counts)
  ))
}

