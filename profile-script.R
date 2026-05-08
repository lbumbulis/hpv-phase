
args   <- commandArgs(trailingOnly=TRUE)
pfail  <- as.numeric(args[1]) # 0.5, 0.7, 0.9
p2     <- as.numeric(args[2]) # 0.17, 0.8
max.M2 <- as.numeric(args[3]) # 3, 7

source("source.R")
source("likelihood.R")

tau <- 4

theta <- set.lambda(theta=list(p2=p2, max.M2=max.M2), tau, prop0=0.5, pfail=pfail)
# Note: This gives lambda1 and lambda2 that are much lower than in the data
# application since the data application has recurrence (1->0) to compete with
# the 1->2 transition, so progressive rates can be higher in that context without
# inflating the end-of-study failure probability.

set.seed(1)
system.time(param.grid <- sim.loglik012.m(theta, tau))

write.csv(param.grid, file=paste0(
  "param_grid_p", p2, "_maxM", max.M2, "_tau", tau, "_pfail", pfail, ".csv"
), row.names=F)
