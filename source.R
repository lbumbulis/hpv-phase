
library(dplyr)

lambda0 <- 1.7
mu0 <- 0.54
phi0 <- 1

lambda1 <- 0.15

lambda2 <- 1.6
mu2 <- 8
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

tau <- 90 # administrative censoring age

# Generate sojourn time in state 0 before -> 1
s01.fn <- function(N) { rexp(N, rate=lambda1) }

# Generate sojourn time in state 1 before -> 2
s10.fn <- function(N, M_type="individual") {
  if (M_type=="common") {
    M <- rnbinom(1, size=phi0, mu=mu0) + 1
  } else if (M_type=="individual") {
    M <- rnbinom(N, size=phi0, mu=mu0) + 1
  }
  rgamma(N, shape=M, rate=lambda0)
}

# Generate sojourn time in state 1 before -> 2
s12.fn <- function(N, M_type="individual") {
  if (M_type=="common") {
    M <- rnbinom(1, size=phi2, mu=mu2) + 1
  } else if (M_type=="individual") {
    M <- rnbinom(N, size=phi2, mu=mu2) + 1
  }
  rgamma(N, shape=M, rate=lambda2)
}

###############################################################################
# Code in here adapted from Fangya's hpv_precanc.R
###############################################################################
# Simulate time to HPV appearance
appear <- function(strt, strt_z){
  n <- length(strt)
  s <- s01.fn(n)
  t_trans <- ifelse(strt_z==2, Inf, strt+s)
  z_trans <- ifelse(strt_z==2, 2, 1)
  return(data.frame(id=seq(n), from=strt, to=t_trans, from_z=strt_z, to_z=z_trans))
}
# Simulate clearance (resolution) of HPV infection
resolve <- function(strt, strt_z){
  n <- length(strt)
  clr <- s10.fn(n)
  prg <- s12.fn(n)
  t_trans <- ifelse(strt_z==2, Inf, strt + pmin(clr, prg))
  z_trans <- ifelse(strt_z==2, 2, ifelse(clr<prg, 0, 2))
  return(data.frame(id=seq(n), from=strt, to=t_trans, from_z=strt_z, to_z=z_trans))
}

# Create natural history (NH)
generate_NH <- function(n) {
  out <- data.frame()
  out <- rbind(out, appear(0, rep(0,n)))       # 1st appearance
  out <- rbind(out, resolve(out$to, out$to_z)) # 1st infection resolved
  
  # Repeat until individuals are too old to be included in studies
  max_iter <- 24
  for (i in 0:max_iter){
    out <- rbind(out, appear(out$to[((2*i+1)*n+1):((2*i+2)*n)], out$to_z[((2*i+1)*n+1):((2*i+2)*n)]))
    out <- rbind(out, resolve(out$to[((2*i+2)*n+1):((2*i+3)*n)], out$to_z[((2*i+2)*n+1):((2*i+3)*n)]))
  }
  # # Check age and state distribution at end of simulated NH
  # sapply(out[(((max_iter+2)*2-1)*n+1):((max_iter+2)*2*n),], summary)
  
  # out$j <- rep(1:(2*(max_iter+2)), each=n) # integer time index
  
  out <- out[which(out$from_z != 2),]
  out <- out[order(out$id, out$to),]
  
  return(out)
}
###############################################################################

# Administratively censor natural history data NH at age tau
# (for modelling based on exact transition times)
admincens <- function(NH) {
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

