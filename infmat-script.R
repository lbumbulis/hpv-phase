
args   <- commandArgs(trailingOnly=TRUE)
mu2    <- as.numeric(args[1]) # 2, 4, 6
max.M2 <- as.numeric(args[2]) # 5, 7
prop0  <- as.numeric(args[3]) # 0.3, 0.5, 0.7
pfail  <- as.numeric(args[4]) # 0.5, 0.7, 0.9

source("source.R")
source("infmat-source.R")

dt.vals <- c(1, 2, 4, 8, 10, 20)
res <- as.data.frame(data.table::rbindlist(lapply(dt.vals, function(dt) {
  temp  <- ESI.fn(mu2, max.M2, prop0, pfail, dt, tau=80)
  var.i <- diag(solve(temp$I))
  out   <- c(dt=dt, var.i)
  names(out)[2:4] <- c("log.lambda1", "log.lambda2", "log.mu2")
  return(as.data.frame(t(out)))
})))

write.csv(res, file=paste0(
  "var_mu", mu2, "_maxM", max.M2, "_prop", prop0, "_pfail", pfail, ".csv"
), row.names=F)
