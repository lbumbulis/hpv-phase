
loglik012.c.l2.term <- function(dat, l2, M) {
  term2 <- dat$d1 * (
      (1-dat$d2) * log(pgamma(dat$du, shape=M, rate=l2, lower.tail=F)) +
      dat$d2 *     log(dgamma(dat$du, shape=M, rate=l2))
  )
  term2[which(dat$du==0)] <- 0 # remove NAs from f(0) in rows where neither 0->1 nor 1->2 was observed
  return(term2)
}

# Log-likelihood for 0->1->2 process, marginal over M
loglik012.m.fn <- function(dat, theta, max.M) {
  l1 <- theta$lambda1
  l2 <- theta$lambda2
  
  # Calculate the lambda1 terms of the log-likelihood
  l1.terms <- -dat$u1*l1 + dat$d1*log(l1)
  # Calculate the lambda2 terms conditional on each possible M
  loglik.c <- l1.terms + sapply(1:max.M, function(M) { loglik012.c.l2.term(dat, l2, M) })
  # Marginalize over M
  log.M.probs <- dnbinom(0:(max.M-1), size=theta$phi2, mu=theta$mu2, log=T)
  sum(matrixStats::rowLogSumExps(sweep(loglik.c, 2, log.M.probs, "+")))
}

sim.loglik012.m <- function(theta, max.M, tau, profile.param, N=10^5) {
  s0  <- s01.fn(N, theta$lambda1)
  s1  <- s12.fn(N, theta)
  
  dat <- data.frame(
    u1 = pmin(s0,    tau),
    u2 = pmin(s0+s1, tau),
    d1 = as.numeric(s0    <= tau),
    d2 = as.numeric(s0+s1 <= tau)
  )
  dat$du <- dat$u2 - dat$u1
  
  log.l2.values   <- round(log(theta$lambda2), 2) - 2 + seq(from=0.2,  by=0.2, length.out=20)
  log.phi2.values <- round(log(theta$phi2), 1)    - 2 + seq(from=0.2,  by=0.2,  length.out=20)
  log.mu2.values  <- round(log(theta$mu2), 1)     - 2 + seq(from=0.2,  by=0.2,  length.out=20)
  
  if (profile.param == "lambda2") {
    param.grid <- expand.grid(log.phi2=log.phi2.values,  log.mu2=log.mu2.values)
  } else if (profile.param == "phi2") {
    param.grid <- expand.grid(log.lambda2=log.l2.values, log.mu2=log.mu2.values)
  } else if (profile.param == "mu2") {
    param.grid <- expand.grid(log.lambda2=log.l2.values, log.phi2=log.phi2.values)
  }
  
  log.start <- c(0,0)
  profile.param.names <- c("lambda1", profile.param)
  names(log.start) <- profile.param.names
  
  # Function to profile out lambda1 and profile.param
  profile.fn <- function(param.list) {
    # If not profiling out lambda2, we can pre-compute the (l2,M)-dependent terms before optim
    if (profile.param != "lambda2") {
      l2 <- param.list$lambda2
      l2.terms <- sapply(1:max.M, function(M) { loglik012.c.l2.term(dat, l2, M) })  # N x max.M matrix
      
      profile.optim.fn <- function(log.profile.vec) {
        theta <- c(as.list(exp(log.profile.vec)), param.list)
        log.M.probs <- dnbinom(0:(max.M-1), size=theta$phi2, mu=theta$mu2, log=T)
        
        l1 <- theta$lambda1
        l1.terms <- -dat$u1*l1 + dat$d1*log(l1) # length N vector
        loglik.c <- l1.terms + l2.terms         # N x max.M matrix
        
        sum(matrixStats::rowLogSumExps(sweep(loglik.c, 2, log.M.probs, "+")))
      }
    } else {
      profile.optim.fn <- function(log.profile.vec) {
        loglik012.m.fn(dat, c(as.list(exp(log.profile.vec)), param.list), max.M)
      }
    }
    res <- optim(log.start, profile.optim.fn, method="BFGS", control=list(fnscale=-1))
    return(exp(res$par))
  } # takes < 1min for profile.param = "phi2" or "mu2"
  
  # Calculate the profile log-likelihood
  param.grid$loglik <- sapply(1:nrow(param.grid), function(rr) {
    print(paste0(Sys.time(), ": row ", rr))
    
    theta <- exp(unlist(param.grid[rr,]))
    names(theta) <- gsub("log.", "", names(theta), fixed=T)
    param.list <- as.list(theta)
    
    profile.hat <- profile.fn(param.list)
    names(profile.hat) <- profile.param.names
    
    loglik012.m.fn(dat, c(as.list(profile.hat), param.list), max.M)
  })
  return(param.grid)
}

