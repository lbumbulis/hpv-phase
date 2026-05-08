
loglik012.c.l2.term <- function(dat, l2, M) {
  term2 <- dat$d1 * (
      (1-dat$d2) * log(pgamma(dat$du, shape=M, rate=l2, lower.tail=F)) +
      dat$d2 *     log(dgamma(dat$du, shape=M, rate=l2))
  )
  term2[which(dat$du==0)] <- 0 # remove NAs from f(0) in rows where neither 0->1 nor 1->2 was observed
  return(term2)
}

# Log-likelihood for 0->1->2 process, marginal over M
loglik012.m.fn <- function(dat, theta) {
  l1 <- theta$lambda1
  l2 <- theta$lambda2
  max.M2 <- theta$max.M2
  
  # Calculate the lambda1 terms of the log-likelihood
  l1.terms <- -dat$u1*l1 + dat$d1*log(l1)
  # Calculate the lambda2 terms conditional on each possible M
  loglik.c <- l1.terms + sapply(1:max.M2, function(M) { loglik012.c.l2.term(dat, l2, M) })
  # Marginalize over M
  log.M2.pmf <- log(M.pmf.fn(theta$p2, max.M2))
  sum(matrixStats::rowLogSumExps(sweep(loglik.c, 2, log.M2.pmf, "+")))
}

sim.loglik012.m <- function(theta.true, tau, N=10^5) {
  s0  <- s01.fn(N, theta.true$lambda1)
  s1  <- s12.fn(N, theta.true)
  
  dat <- data.frame(
    u1 = pmin(s0,    tau),
    u2 = pmin(s0+s1, tau),
    d1 = as.numeric(s0    <= tau),
    d2 = as.numeric(s0+s1 <= tau)
  )
  dat$du <- dat$u2 - dat$u1
  
  log.l2.values <- round(log(theta.true$lambda2), 2) - 0.2 + seq(from=0.02, by=0.02, length.out=20)
  logit.p2.values <- round(logit(theta.true$p2), 2)  - 0.5 + seq(from=0.05, by=0.05, length.out=20)
  param.grid <- expand.grid(log.lambda2=log.l2.values, logit.p2=logit.p2.values)
  
  log.l1.start <- c(lambda1=0)
  
  # Function to profile out lambda1 and profile.param
  profile.fn <- function(param.list) {
    # Pre-compute the l1-independent terms before optim
    l2 <- param.list$lambda2
    p2 <- param.list$p2
    max.M2 <- param.list$max.M2
    log.M2.pmf <- log(M.pmf.fn(p2, max.M2))
    
    l2.terms <- sapply(1:max.M2, function(M) { loglik012.c.l2.term(dat, l2, M) })  # N x max.M2 matrix
    
    profile.optim.fn <- function(log.l1) {
      l1 <- exp(log.l1)
      l1.terms <- -dat$u1*l1 + dat$d1*log(l1) # length N vector
      loglik.c <- l1.terms + l2.terms         # N x max.M2 matrix
      
      sum(matrixStats::rowLogSumExps(sweep(loglik.c, 2, log.M2.pmf, "+")))
    }
    res <- optim(log.l1.start, profile.optim.fn, method="BFGS", control=list(fnscale=-1))
    return(expify(res$par))
  } # takes < 1min for profile.param="p2"
  
  # Calculate the profile log-likelihood
  param.grid$loglik <- sapply(1:nrow(param.grid), function(rr) {
    print(paste0(Sys.time(), ": row ", rr))
    
    theta <- unlist(param.grid[rr,])
    names(theta) <- gsub("logit.", "", names(theta), fixed=T)
    names(theta) <- gsub("log.", "", names(theta), fixed=T)
    theta <- expify(theta)
    
    param.list <- c(as.list(theta), list(max.M2=theta.true$max.M2))
    
    loglik012.m.fn(dat, c(list(lambda1=profile.fn(param.list)), param.list))
  })
  return(param.grid)
}

