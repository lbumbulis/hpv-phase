
library(dplyr)

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
#      this results in a true gamma mixture distribution in each dataset.
# * However, with panel data we need to use the Markov phase-type representation
#   of the multistate process, which requires a fixed M at least within each
#   pair of (i, transition type) to force a common dimension on the transition
#   probability matrices being multiplied together.
# * For now when we have panel data we will use a common M across all individuals.
#   If we fix M for each pair of (i, transition type) but allow it to vary across
#   individuals then it acts as a random effect; we may consider this later.

###############################################################################
# Competing gammas representation
###############################################################################
# Generate sojourn time in state 0 before -> 1
s01.fn <- function(N) { rexp(N, rate=lambda1) }

# Generate sojourn time in state 1 before -> 0
s10.fn <- function(N, M=NA) {
  if (is.na(M)) {
    M <- rnbinom(N, size=phi0, mu=mu0) + 1
  }
  rgamma(N, shape=M, rate=lambda0)
}

# Generate sojourn time in state 1 before -> 2
s12.fn <- function(N, M=NA) {
  if (is.na(M)) {
    M <- rnbinom(N, size=phi2, mu=mu2) + 1
  }
  rgamma(N, shape=M, rate=lambda2)
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

