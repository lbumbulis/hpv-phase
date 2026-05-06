
args <- commandArgs(trailingOnly=TRUE)
iter <- as.numeric(args[1])

library(cmdstanr)

source("source.R")

sim.seeds <- readRDS("sim_seeds_nsim1000.rds")

.Random.seed <- sim.seeds[[iter]]

N <- 1000
dat <- admincens(generate_NH(N), tau=90)
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
  
  max_M = 50,
  
  grainsize = 1
)

# set_cmdstan_path("C:/Users/lbumb/.cmdstan/cmdstan-2.38.0")
set_cmdstan_path("/home/lsbumbul/.cmdstan/cmdstan-2.38.0")

init_fn <- function() {
  list(
    log_lambda1  = runif(1, min=-3, max=0),
    
    log_mean_s10 = runif(1, min=-4, max=-1),
    log_mu0      = runif(1, min=-4, max=2),
    log_phi0     = runif(1, min=-1, max=2),
    
    log_mean_s12 = runif(1, min=1,  max=3),
    log_mu2      = runif(1, min=0,  max=4),
    log_phi2     = runif(1, min=1,  max=3.3)
  )
}

stan_model <- cmdstan_model("model012-recur.stan", cpp_options=list(stan_threads=T))

model_fit <- stan_model$sample(
  stan_data,
  init = init_fn,
  iter_sampling = 1000,
  iter_warmup = 500,
  parallel_chains = 4,
  threads_per_chain = 4,
  refresh = 5
)
model_fit$save_object(file=paste0("./results/model012-recur_iter", iter, ".rds"))

# test_fit <- stan_model$sample(
#   stan_data,
#   init = init_fn,
#   iter_sampling = 1000,
#   iter_warmup = 500,
#   chains = 2,
#   parallel_chains = 2,
#   refresh = 1,
#   save_warmup = T
# )
# # For (lambda, mu) param, this takes ~ 50min without "dense_e"; with "dense_e" it seems even slower.
# 
# test_fit$draws() |> posterior::as_draws_df() |>
#   select(starts_with("log_")) |>
#   pairs()
# test_fit$draws() |> posterior::as_draws_df() |>
#   select(starts_with("log_")) |>
#   cor()
# # log_lambda0 is highly (positively) correlated with log_mu0 (~0.96),
# # and likewise for log_lambda2 and log_mu2.
# # - More phases => need to move through the phases at a faster rate.


