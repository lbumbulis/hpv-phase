
args <- commandArgs(trailingOnly=TRUE)
iter <- as.numeric(args[1])

library(cmdstanr)

source("source.R")

sim.seeds <- readRDS("sim_seeds_nsim1000.rds")

.Random.seed <- sim.seeds[[iter]]

N <- 1000
s12 <- s12.fn(N)

stan_data <- list(
  N = N,
  s12 = s12,
  max_M = 100
)

# set_cmdstan_path("C:/Users/lbumb/.cmdstan/cmdstan-2.38.0")
set_cmdstan_path("/home/lsbumbul/.cmdstan/cmdstan-2.38.0")

init_fn <- function() {
  list(
    lambda2  = runif(1, min=0, max=3),
    log_mu2  = runif(1, min=0, max=4),
    log_phi2 = runif(1, min=1, max=4)
  )
}

stan_model <- cmdstan_model("model12.stan")
model_fit <- stan_model$sample(
  stan_data,
  init = init_fn,
  iter_sampling = 1500,
  iter_warmup = 500,
  parallel_chains = 4,
  refresh = 10
)

model_fit$save_object(file=paste0("./results/model12_iter", iter, ".rds"))

