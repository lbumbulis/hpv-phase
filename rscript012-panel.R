
args <- commandArgs(trailingOnly=TRUE)
iter <- as.numeric(args[1])

library(cmdstanr)

source("source.R")

sim.seeds <- readRDS("sim_seeds_nsim1000.rds")

.Random.seed <- sim.seeds[[iter]]

N <- 1000
M2 <- rnbinom(N, size=phi2, mu=mu2) + 1

s01 <- s01.fn(N)
s12 <- s12.fn(N, M=M2)

dat <- dfify(s01, s12)
pdat <- panelize(dat, visit.freq=3, tau=200)

stan_data <- c(process.panel.012(pdat), list(max_M=25))

init_fn <- function() {
  list(
    log_lambda1  = runif(1, min=-3,  max=0),
    log_mean_s12 = runif(1, min=0,   max=2),
    log_mu2      = runif(1, min=0.5, max=2),
    log_phi2     = runif(1, min=2,   max=3.3)
  )
}

# set_cmdstan_path("C:/Users/lbumb/.cmdstan/cmdstan-2.38.0")
set_cmdstan_path("/home/lsbumbul/.cmdstan/cmdstan-2.38.0")



stan_model <- cmdstan_model("model012-panel.stan")
model_fit <- stan_model$sample(
  stan_data,
  init = init_fn,
  iter_sampling = 1500,
  iter_warmup = 500,
  parallel_chains = 4,
  refresh = 10
)

model_fit$save_object(file=paste0("./results/model012-panel_iter", iter, ".rds"))

