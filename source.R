
mu2 <- 8
phi2 <- 20

lambda1 <- 0.15
lambda2 <- 1.5

stan_params12 <- data.frame(
  param = c("lambda2", "log_mu2", "log_phi2"),
  value = c(lambda2, log(mu2), log(phi2))
)

stan_params012 <- rbind(
  data.frame(param="lambda1", value=lambda1), stan_params12
)

# Generate sojourn time in state 0 before -> 1
s01.fn <- function(N) { rexp(N, rate=lambda1) }

# Generate sojourn time in state 1 before -> 2
s12.fn <- function(N, M_type) {
  if (M_type=="common") {
    M <- rnbinom(1, size=phi2, mu=mu2)+1
  } else if (M_type=="individual") {
    M <- rnbinom(N, size=phi2, mu=mu2)+1
  }
  rgamma(N, shape=M, rate=lambda2)
}
