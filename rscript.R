
args <- commandArgs(trailingOnly=TRUE)
iter <- as.numeric(args[1])

library(cmdstanr)

source("source.R")

sim.seeds <- readRDS("sim_seeds_nsim1000.rds")

.Random.seed <- sim.seeds[[iter]]

pdat <- panelize(
  generate_NH(n=2000),
  visit.freq=1, tau=90
)

idx00 <- which(pdat$from_z==0 & pdat$to_z==0)

stan_data <- list(
  n00 = length(idx00),
  s00 = pdat$dt[idx00],
  
  max_M0 = 5,
  max_M2 = 12
)
stan_data <- c(stan_data, process.panel(pdat))

init_fn <- function() {
  list(
    log_lambda1  = runif(1, min=-3, max=0),
    
    log_mean_s10 = runif(1, min=-4, max=-1),
    log_mu0      = runif(1, min=-2, max=1),
    log_phi0     = runif(1, min=-1, max=1),
    
    log_mean_s12 = runif(1, min=1,  max=3),
    log_mu2      = runif(1, min=0,  max=2),
    log_phi2     = runif(1, min=1,  max=3.3)
  )
}

# set_cmdstan_path("C:/Users/lbumb/.cmdstan/cmdstan-2.38.0")
set_cmdstan_path("/home/lsbumbul/.cmdstan/cmdstan-2.38.0")

stan_model <- cmdstan_model("model.stan")

model_fit <- stan_model$sample(
  stan_data,
  init = init_fn,
  iter_sampling = 1000,
  iter_warmup = 1000,
  parallel_chains = 4,
  refresh = 5,
  save_warmup = TRUE
)
model_fit$save_object(file=paste0("./results/model_iter", iter, ".rds"))



