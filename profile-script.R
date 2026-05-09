
args   <- commandArgs(trailingOnly=TRUE)
pfail  <- as.numeric(args[1]) # 0.5, 0.7, 0.9
mu2    <- as.numeric(args[2]) # 0.5, 6.5
# p2     <- as.numeric(args[2]) # 0.17, 0.8
max.M2 <- as.numeric(args[3]) # 2, 6

source("source.R")
source("profile-source.R")

tau <- 4

theta <- set.lambda(theta=list(mu2=mu2, max.M2=max.M2), tau, prop0=0.5, pfail=pfail)
# theta <- set.lambda(theta=list(p2=p2, max.M2=max.M2), tau, prop0=0.5, pfail=pfail)

set.seed(1)
system.time(param.grid <- sim.loglik012.m(theta, tau))

write.csv(param.grid, file=paste0(
  "param_grid_mu", mu2,
  # "param_grid_p", p2,
  "_maxM", max.M2, "_tau", tau, "_pfail", pfail, ".csv"
), row.names=F)
