
library(expm)
library(matrixStats)
library(numDeriv)

# Q matrix for given M: state 0, then (M+1) phases of state 1, then state 2.
# Total n.states = M+3.
Q.fn <- function(M, l1, l2) {
  n <- M+3
  Q <- matrix(0, n, n)
  Q[1, 2] <- l1                             # state 0 -> first phase
  for (m in 1:(M+1)) { Q[m+1, m+2] <- l2 }  # phase -> next phase (or state 2)
  diag(Q) <- -rowSums(Q)
  Q
}

# ---------------------------------------------------------------------------
# Forward algorithm
# ---------------------------------------------------------------------------

# Log-probability of observed sequence given M, lambda1, lambda2.
# Computed via forward algorithm.
# seq.obs: integer vector of length n.visits, each in {0, 1, 2};
#          first element corresponds to first visit (not the known t=0 start).
seq.logprob.fn <- function(seq.obs, M, l1, l2, dt) {
  n.states <- M+3
  P <- expm(Q.fn(M, l1, l2) * dt)
  
  fwd    <- rep(0, n.states)
  fwd[1] <- 1.0   # known to start in state 0
  log.prob <- 0.0
  
  for (obs in seq.obs) {
    fwd <- fwd %*% P
    
    p.obs <- switch(obs + 1,
                    fwd[1],                    # obs == 0
                    sum(fwd[2:(n.states-1)]),  # obs == 1 (sum over phases)
                    fwd[n.states]              # obs == 2
    )
    if (p.obs <= 0) { return(-Inf) } # guard against underflow
    log.prob <- log.prob + log(p.obs)
    
    # Zero out entries for states that are incompatible with the observed state
    if (obs == 0) {
      fwd[-1]        <- 0
    } else if (obs == 1) {
      fwd[1]         <- 0
      fwd[n.states]  <- 0
    } else {
      fwd[-n.states] <- 0
    }
    # Normalize to use as starting probabilities for next time point
    fwd <- fwd / p.obs
  }
  log.prob
}

# Log-probability of a sequence, marginalized over M ~ TruncPois(mu, max.M).
# theta: parameter vector (log.lambda1, log.lambda2, log.mu2)
seq.logprob.marginal <- function(theta, seq.obs, max.M2, dt) {
  l1  <- exp(theta[1])
  l2  <- exp(theta[2])
  mu2 <- exp(theta[3])
  
  log.M.pmf <- log(M.pmf.fn(mu2, max.M2))
  log.seq.M <- sapply(0:max.M2, function(M) { # log-likelihood of sequence given M
    seq.logprob.fn(seq.obs, M, l1, l2, dt)
  })
  logSumExp(log.seq.M + log.M.pmf)
}

# ---------------------------------------------------------------------------
# Enumerate valid sequences
# ---------------------------------------------------------------------------

# All non-decreasing sequences of length n.visits from {0, 1, 2}.
# Pattern: a zeros, b ones, c twos with a + b + c = n.visits.
gen.sequences <- function(n.visits) {
  # No. of sequences is (n.visits choose 2) since we choose the locations of
  # the last 0 observation and the last 1 observation in the sequence.
  seqs <- vector("list", choose(n.visits+2, 2))
  idx  <- 1L
  for (a in 0:n.visits) {
    for (b in 0:(n.visits-a)) {
      c.  <- n.visits - a - b
      seqs[[idx]] <- c(rep(0L, a), rep(1L, b), rep(2L, c.))
      idx <- idx + 1L
    }
  }
  seqs[1:(idx - 1L)]
}

# ---------------------------------------------------------------------------
# Expected score vector and information matrix
# ---------------------------------------------------------------------------

# Parameters are (log.lambda1, log.lambda2, log.mu); derivatives are taken
# on this log scale so that each parameter is unconstrained.

# n.visits is the number of observations after t=0

ESI.fn <- function(mu2, max.M2, prop0, pfail, dt, tau) {
  n.visits <- tau/dt
  theta.list  <- set.lambda(
    theta=list(mu2=mu2, max.M2=max.M2), tau=tau, prop0=prop0, pfail=pfail
  )
  theta       <- c(log(theta.list$lambda1), log(theta.list$lambda2), log(mu2))
  param.names <- c("log.lambda1", "log.lambda2", "log.mu2")
  n.params    <- length(theta)
  sequences   <- gen.sequences(n.visits)
  n.seqs      <- length(sequences)
  cat(sprintf("Evaluating %d sequences x %d M values...\n", n.seqs, max.M2+1))
  
  # Sequence probabilities at theta (weights for expectation)
  l.probs <- sapply(sequences, function(s) { seq.logprob.marginal(theta, s, max.M2, dt) })
  probs   <- exp(l.probs)
  # cat(sprintf("Sanity check - sum of sequence probabilities: %.8f\n", sum(probs)))
  
  # For each sequence, compute Hessian of log p(seq | theta) and accumulate
  # the weighted sum. numDeriv::hessian uses Richardson extrapolation by default.
  # Likewise for the score vector.
  S.vec <- numeric(n.params)
  I.mat <- matrix(0, n.params, n.params)
  for (i in seq_along(sequences)) {
    S <- grad(seq.logprob.marginal, theta, seq.obs=sequences[[i]], max.M=max.M2, dt=dt)
    H <- hessian(seq.logprob.marginal, theta, seq.obs=sequences[[i]], max.M=max.M2, dt=dt)
    S.vec <- S.vec + probs[i]*S
    I.mat <- I.mat - probs[i]*H
  }
  
  names(S.vec)    <- param.names
  dimnames(I.mat) <- list(param.names, param.names)
  list(S=S.vec, I=I.mat, seq.probs=setNames(probs, sapply(sequences, paste, collapse="")))
}


# ## Quick check
# # Timing
# system.time(foo3 <- ESI.fn(mu2=6.5, max.M2=6, prop0=0.5, pfail=0.5, dt=3, tau=30)) # < 10s ; 66 seqs
# system.time(foo2 <- ESI.fn(mu2=6.5, max.M2=6, prop0=0.5, pfail=0.5, dt=2, tau=30)) # < 15s ; 136 seqs
# system.time(foo1 <- ESI.fn(mu2=6.5, max.M2=6, prop0=0.5, pfail=0.5, dt=1, tau=30)) # < 1min; 496 seqs
# system.time(bar  <- ESI.fn(mu2=6.5, max.M2=6, prop0=0.5, pfail=0.5, dt=1, tau=80)) # < 7min; 3321 seqs
# # => time cost appears to be linear in the number of sequences, as expected
# 
# # Check E(S)=0
# isTRUE(all.equal(round(unname(foo3$S), 10), rep(0,3)))
# isTRUE(all.equal(round(unname(foo2$S), 10), rep(0,3)))
# isTRUE(all.equal(round(unname(foo1$S), 10), rep(0,3)))
# 
# # Check that SEs get smaller as dt decreases
# sqrt(diag(solve(foo3$I)))
# sqrt(diag(solve(foo2$I)))
# sqrt(diag(solve(foo1$I)))


