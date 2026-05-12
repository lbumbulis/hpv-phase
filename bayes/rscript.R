
args <- commandArgs(trailingOnly=TRUE)
iter <- as.numeric(args[1])

library(cmdstanr)

source("source.R")

lambda1 <- 0.15

# 1->0 parameters
lambda0 <- 1.6
mu0 <- 0.5
max.M0 <- 2

# 1->2 parameters
lambda2 <- 1
mu2 <- 6.5
max.M2 <- 6

stan_params <- data.frame(
  param = paste0("log_", c("lambda1", "lambda0", "lambda2", "mu0", "mu2")),
  value = log(c(lambda1, lambda0, lambda2, mu0, mu2))
)

sim.seeds <- readRDS("sim_seeds_nsim1000.rds")

.Random.seed <- sim.seeds[[iter]]

N <- 1000
visit.freq <- 3
se <- 0.25

M0 <- rM(N, mu=mu0, max.M=max.M0) # truncated Poisson distribution
M2 <- rM(N, mu=mu2, max.M=max.M2)

pdat <- panelize(generate_NH(N, theta=list(
  lambda1=lambda1, lambda0=lambda0, mu0=mu0, max.M0=max.M0, lambda2=lambda2, mu2=mu2, max.M2=max.M2
), M=list(M0, M2)), visit.freq, tau=45)

stan_data <- c(process.panel(pdat, visit.freq), list(max_M0=max.M0, max_M2=max.M2, se=se))

init_fn <- function() {
  list(
    log_lambda1 = runif(1, min=-2.5, max=-1),
    
    log_lambda0 = runif(1, min=-0.5, max=1.5),
    log_mu0     = runif(1, min=-1.5, max=0),
    
    log_lambda2 = runif(1, min=-1,   max=1),
    log_mu2     = runif(1, min=0.5,  max=2)
  )
}

# set_cmdstan_path("C:/Users/lbumb/.cmdstan/cmdstan-2.38.0")
set_cmdstan_path("/home/lsbumbul/.cmdstan/cmdstan-2.38.0")

stan_model <- cmdstan_model("model.stan")

model_fit <- stan_model$sample(
  stan_data,
  init            = init_fn,
  iter_sampling   = 1000,
  iter_warmup     = 1000,
  # adapt_delta     = 0.99,
  chains          = 4,
  parallel_chains = 4,
  refresh         = 5,
  save_warmup     = TRUE #,
  # metric          = "dense_e"  # may help handle correlation between parameters
)
# With N=1000, visit.freq=3, tau=45, metric="diag_e", default adapt_delta (0.8):
# - Takes about 7min with 0.25 SE on priors, and estimation is good
# - Takes about 8min with 0.5  SE on priors, and estimation is a bit off
#   (e.g. hat(log_lambda2)=0.225) with 1% divergences
# - Takes about 11min with 1   SE on priors, and estimation is noticeably worse
#   (e.g. hat(log_lambda2)=0.359) with 20% divergences
# Moving to adapt_delta=0.9:
# - Takes about 9min with 0.5 SE on priors and produces 1% divergences;
#   estimates are barely different from those obtained under default adapt_delta
# Moving to adapt_delta=0.99:
# - Takes about 23min with 0.5 SE on priors and produces 1% divergences;
#   helps a bit, but not much (e.g. hat(log_lambda2)=0.215)
# Going back to default adapt_delta (0.8) but now switching to metric="dense_e":
# _ Takes about 3min with 0.5 SE on priors and produces 1% divergences;
#   estimates are worse than with metric="diag_e"

model_fit$save_object(file=paste0("./results/model_se", se, "_iter", iter, ".rds"))



