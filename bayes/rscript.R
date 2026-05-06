
args <- commandArgs(trailingOnly=TRUE)
iter <- as.numeric(args[1])

library(cmdstanr)

source("source.R")

lambda0 <- 1.7
mu0 <- 0.54
phi0 <- 1

lambda1 <- 0.15

lambda2 <- 1.6
mu2 <- 5
phi2 <- 20

sim.seeds <- readRDS("sim_seeds_nsim1000.rds")

.Random.seed <- sim.seeds[[iter]]

N <- 1000
visit.freq <- 3

M0 <- rnbinom(N, size=phi0, mu=mu0) + 1
M2 <- rnbinom(N, size=phi2, mu=mu2) + 1

pdat <- panelize(generate_NH(N, list(M0, M2)), visit.freq, tau=90)

# qnbinom(0.99, size=phi0, mu=mu0) + 1 # used to determine good value for max_M0
# qnbinom(0.99, size=phi2, mu=mu2) + 1 # ... max_M2

stan_data <- c(process.panel(pdat, visit.freq), list(max_M0=5, max_M2=15))

init_fn <- function() {
  list(
    log_lambda1  = runif(1, min=-3,   max=0),
    
    log_mean_s10 = runif(1, min=-2,   max=-1),
    log_mu0      = runif(1, min=-1,   max=0),
    log_phi0     = runif(1, min=-0.5, max=0.5),
    
    log_mean_s12 = runif(1, min=1,    max=2),
    log_mu2      = runif(1, min=1,    max=2),
    log_phi2     = runif(1, min=2.5,  max=3.3)
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
  adapt_delta = 0.99,
  chains = 2,
  parallel_chains = 2,
  refresh = 5,
  save_warmup = TRUE
)
model_fit$save_object(file=paste0("./results/model_iter", iter, ".rds"))



