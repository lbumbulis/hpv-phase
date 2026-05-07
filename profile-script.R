
args <- commandArgs(trailingOnly=TRUE)
p2   <- as.numeric(args[1]) # 0.5, 0.7, 0.9
profile.param <-   args[2]  # "lambda2", "phi2", "mu2"

source("source.R")
source("likelihood.R")

phi2  <- 1
mu2   <- 3
max.M <- 25  # chosen using qnbinom(0.999, size=1, mu=3)+1
tau   <- 10

theta <- set.lambda(theta=list(phi2=phi2, mu2=mu2), max.M, tau, prop0=0.5, p2=p2)
# Note: This gives lambda1 and lambda2 that are much lower than in the data
# application since the data application has recurrence (1->0) to compete with
# the 1->2 transition, so progressive rates can be higher in that context without
# inflating the end-of-study failure probability.

set.seed(1)
system.time(param.grid <- sim.loglik012.m(theta, max.M, tau, profile.param))

write.csv(param.grid, file=paste0(
  "param_grid_phi", phi2, "_mu", mu2, "_tau", tau, "_p2-", p2, "_profile-", profile.param, ".csv"
), row.names=F)
