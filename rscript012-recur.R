
args <- commandArgs(trailingOnly=TRUE)
iter <- as.numeric(args[1])

library(cmdstanr)

source("source.R")

sim.seeds <- readRDS("sim_seeds_nsim1000.rds")

.Random.seed <- sim.seeds[[iter]]

N <- 1000
dat <- admincens(generate_NH(N))
dat$s <- dat$to - dat$from

idx01 <- which(dat$from_z==0 & dat$to_z==1)
idx10 <- which(dat$from_z==1 & dat$to_z==0)
idx12 <- which(dat$to_z==2)

idx00 <- which(dat$from_z==0 & dat$to_z==0) # admin censored
idx11 <- which(dat$from_z==1 & dat$to_z==1)

stan_data <- list(
  n01 = length(idx01),
  n10 = length(idx10),
  n12 = length(idx12),
  n00 = length(idx00),
  n11 = length(idx11),
  
  s01 = dat$s[idx01],
  s10 = dat$s[idx10],
  s12 = dat$s[idx12],
  s00 = dat$s[idx00],
  s11 = dat$s[idx11],
  
  max_M = 100
)

# set_cmdstan_path("C:/Users/lbumb/.cmdstan/cmdstan-2.38.0")
set_cmdstan_path("/home/lsbumbul/.cmdstan/cmdstan-2.38.0")

init_fn <- function() {
  list(
    log_lambda0 = runif(1, min=0,  max=1.5),
    log_mu0     = runif(1, min=-4, max=2),
    log_phi0    = runif(1, min=-1, max=2),
    
    log_lambda1 = runif(1, min=-3, max=0),
    
    log_lambda2 = runif(1, min=0,  max=1.5),
    log_mu2     = runif(1, min=0,  max=4), # could tighten to (0.5, 2.5)
    log_phi2    = runif(1, min=1,  max=4)  # could tighten to (2, 3.5)
  )
}

stan_model <- cmdstan_model("model012-recur.stan")

model_fit <- stan_model$sample(
  stan_data,
  init = init_fn,
  iter_sampling = 1500,
  iter_warmup = 500,
  parallel_chains = 4,
  refresh = 5
)
model_fit$save_object(file=paste0("./results/model012-recur_iter", iter, ".rds"))

