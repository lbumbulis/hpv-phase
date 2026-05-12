
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

source("source.R")

theme.component <- theme_minimal(base_size=12) +
  theme(
    plot.title    = element_text(face="bold"),
    strip.text    = element_text(face="bold", size=10),
    panel.border  = element_rect(fill=NA, linewidth=0.5),
    strip.background = element_rect(fill="grey95")
  )

rho0 <- 1/0.75
kappa0 <- 0.74
rho2 <- 1/6.25
kappa2 <- 2.43

###############################################################################
# Truncated Poisson distributions as a function of mu, lambda, and max.M
###############################################################################
plot.tpois <- function() {
  mu.vals    <- 0:8
  mu.vals[1] <- 0.1
  max.M.vals <- 2:7
  
  param.grid <- expand.grid(mu=mu.vals, max.M=max.M.vals)
  df <- as.data.frame(data.table::rbindlist(lapply(1:nrow(param.grid), function(r) {
    mu    <- param.grid$mu[r]
    max.M <- param.grid$max.M[r]
    data.frame(mu=mu, max.M=max.M, M=0:max.M, prob=M.pmf.fn(mu, max.M))
  })))
  
  df <- df |> mutate(
    mu_label = factor(paste0("mu==", mu), levels=paste0("mu==", mu.vals)),
    max_M_label = factor(paste0("M^'*'==", max.M), levels=paste0("M^'*'==", max.M.vals))
  )
  
  ggplot(df, aes(x=M, y=prob)) +
    geom_col(colour="black", fill="steelblue", alpha=0.3) +
    facet_grid(max_M_label ~ mu_label, labeller=label_parsed) +
    labs(y="PMF") +
    theme.component
}

plot.tpois()

###############################################################################
# Sojourn time distributions as a function of lambda and mu or p
###############################################################################
# Plot sojourn distribution marginal over M
plot.sojourn.M <- function(lambda.vals, param.vals, max.x, target,
                           M.dist, max.M=NULL, out.col="density") {
  x.vals <- seq(0.01, max.x, by=0.01)
  max.M.vals <- c(2,3,5,7,9)
  if (is.null(max.M)) {
    test.params <- expand.grid(x=x.vals, lambda=lambda.vals, param=param.vals, max.M=max.M.vals)
  } else {
    test.params <- expand.grid(x=x.vals, lambda=lambda.vals, param=param.vals)
  }
  
  df <- test.params
  df <- cbind(df, t(sapply(1:nrow(df), function(r) {
    if (is.null(max.M)) { max.M <- df$max.M[r] }
    M.pmf <- M.pmf.fn(df$param[r], max.M, M.dist)
    
    x <- df$x[r]
    dens <- dgamma(x, shape=1:(max.M+1), rate=df$lambda[r])
    surv <- pgamma(x, shape=1:(max.M+1), rate=df$lambda[r], lower.tail=F)
    
    t(M.pmf) %*% cbind(dens, surv, dens/surv)
  })))
  names(df)[ncol(df)-(2:0)] <- c("density", "survival", "hazard")
  df <- df |>
    mutate(
      d10 = dweibull(x, shape=kappa0, scale=1/rho0),
      d12 = dweibull(x, shape=kappa2, scale=1/rho2),
      lambda_label = factor(paste0("lambda==", lambda), levels=paste0("lambda==", lambda.vals))
    )
  
  col.idx <- which(names(df)==out.col)
  names(df)[col.idx] <- "y"
  target.idx <- which(names(df)==target)
  names(df)[target.idx] <- "d"
  
  if (out.col=="density") {
    y.fn <- "f"
  } else if (out.col=="hazard") {
    y.fn <- "h"
  } else if (out.col=="survival") {
    y.fn <- "1 - F"
  }
  
  if (target=="none") {
    target_component <- geom_blank()
  } else {
    target_component <- geom_line(aes(y=d), colour="grey70", linetype="dashed", linewidth=0.6)
  }
  
  if (M.dist=="tpois") {
    df <- df |> mutate(param_label = factor(paste0("mu==", param), levels=paste0("mu==", param.vals)))
    if (is.null(max.M)) {
      lab.component <- labs(
        title = bquote(bold("Erlang mixture" ~ .(out.col) * ": rate" ~ lambda *
                              ", shape = 1 + M;  M ~ tPois(" * mu * "," ~ M^"*" * ")")),
        x = "x",
        y = paste0(y.fn, "(x)")
      )
    } else {
      lab.component <- labs(
        title = bquote(bold("Erlang mixture" ~ .(out.col) * ": rate" ~ lambda *
                              ", shape = 1 + M;  M ~ tPois(" * mu * "," ~ M^"*"==.(max.M) * ")")),
        x = "x",
        y = paste0(y.fn, "(x)")
      )
    }
  } else if (M.dist=="binom") {
    df <- df |> mutate(param_label = factor(paste0("p==", param), levels=paste0("p==", param.vals)))
    lab.component <- labs(
      title = bquote(bold("Erlang" ~ .(out.col) * ": rate" ~ lambda *
                          ", shape = 1 + M;  M ~ Bin(" * n==.(max.M) * "," ~ p * ")")),
      x = "x",
      y = paste0(y.fn, "(x)")
    )
  }
  
  if (is.null(max.M)) {
    g <- ggplot(df, aes(x=x, y=y, colour=factor(max.M)))
    n <- length(max.M.vals)
    scale.component <- scale_colour_manual(
      name   = expression(M^"*"),
      values = rev(brewer.pal(n+3, "Blues")[(1:n)+3])
    )
  } else {
    g <- ggplot(df, aes(x=x, y=y, colour=factor(lambda), fill=factor(lambda)))
    scale.component <-
      scale_colour_brewer(palette="Dark2", guide="none") +
      scale_fill_brewer(palette="Dark2", guide="none")
  }
  
  g + target_component +
    geom_line(linewidth=0.8) +
    facet_grid(param_label ~ lambda_label, labeller=label_parsed) +
    scale.component +
    lab.component +
    theme.component +
    theme(plot.title = element_text(size=12))
}

### Try to emulate 1->0 density
lambda.vals <- c(1.6, 2, 2.5, 3)
mu.vals <- c(0.1, 0.5, 1)
plot.sojourn.M(lambda.vals, mu.vals, max.M=2, max.x=5, target="d10", M.dist="tpois")
# plot.sojourn.M(lambda.vals, mu.vals, max.M=3, max.x=5, target="d10", M.dist="tpois")
# plot.sojourn.M(lambda.vals, mu.vals, max.M=4, max.x=5, target="d10", M.dist="tpois")
# save as params10_tpois.png (Width=550, Height=450)

lambda.vals <- c(1.25, 1.5, 2)
p.vals <- c(0.1, 0.15, 0.2)
# plot.sojourn.M(lambda.vals, p.vals, max.M=1, max.x=5, target="d10", M.dist="binom")
plot.sojourn.M(lambda.vals, p.vals, max.M=2, max.x=5, target="d10", M.dist="binom")
# plot.sojourn.M(lambda.vals, p.vals, max.M=3, max.x=5, target="d10", M.dist="binom")
# plot.sojourn.M(lambda.vals, p.vals, max.M=4, max.x=5, target="d10", M.dist="binom")


### Try to emulate 1->2 density
lambda.vals <- c(0.5, 1, 1.5)
mu.vals <- c(5, 6, 6.5)
plot.sojourn.M(lambda.vals, mu.vals, max.M=6, max.x=12, target="d12", M.dist="tpois")
# plot.sojourn.M(lambda.vals, mu.vals, max.M=7, max.x=12, target="d12", M.dist="tpois")
# plot.sojourn.M(lambda.vals, mu.vals, max.M=9, max.x=12, target="d12", M.dist="tpois")
# save as params12_tpois.png (Width=600, Height=500)

lambda.vals <- c(0.8, 1, 1.5)
p.vals <- c(0.7, 0.75, 0.8, 0.9)
plot.sojourn.M(lambda.vals, p.vals, max.M=6, max.x=12, target="d12", M.dist="binom")
plot.sojourn.M(lambda.vals, p.vals, max.M=7, max.x=12, target="d12", M.dist="binom")

### Explore variety of shapes
plot.sojourn.M(lambda.vals, mu.vals, max.x=8, target="none", M.dist="tpois")
# save as erlang_tpois.png (Width=670, Height=550)


###############################################################################
# Expected sojourn time as a function of mu or p
###############################################################################
ES.fn <- function(lambda, param, max.M, M.dist, squared=F) {
  M.pmf <- M.pmf.fn(param, max.M, M.dist)
  
  if (squared) {
    integrand <- Vectorize(function(x) {
      M.pmf %*% dgamma(x, shape=1:(max.M+1), rate=lambda) * x^2
    })
  } else {
    integrand <- Vectorize(function(x) {
      M.pmf %*% dgamma(x, shape=1:(max.M+1), rate=lambda) * x
      # M.pmf %*% pgamma(x, shape=1:(max.M+1), rate=lambda, lower.tail=F) # same
    })
  }
  integrate(integrand, 0, Inf)$value
}

VarS.fn <- function(lambda, param, max.M, M.dist) {
  ES.fn(lambda, param, max.M, M.dist, squared=T) - ES.fn(lambda, param, max.M, M.dist, squared=F)^2
}

### Try to emulate 1->0 density
# Mean sojourn time
ES.fn(lambda=2,   param=0.1,  max.M=2, M.dist="tpois") # 0.550
ES.fn(lambda=1.6, param=0.5,  max.M=2, M.dist="tpois") # 0.913 *
ES.fn(lambda=1.5, param=0.1,  max.M=2, M.dist="binom") # 0.800
ES.fn(lambda=1.5, param=0.15, max.M=2, M.dist="binom") # 0.867
ES.fn(lambda=1.5, param=0.17, max.M=2, M.dist="binom") # 0.893 **
ES.fn(lambda=1.5, param=0.2,  max.M=2, M.dist="binom") # 0.933
1/rho0 * gamma(1 + 1/kappa0)                           # 0.903

# Variance of sojourn time
VarS.fn(lambda=2,   param=0.1,  max.M=2, M.dist="tpois") # 0.300
VarS.fn(lambda=1.6, param=0.5,  max.M=2, M.dist="tpois") # 0.728 *
VarS.fn(lambda=1.5, param=0.1,  max.M=2, M.dist="binom") # 0.613
VarS.fn(lambda=1.5, param=0.15, max.M=2, M.dist="binom") # 0.691
VarS.fn(lambda=1.5, param=0.17, max.M=2, M.dist="binom") # 0.721 **
VarS.fn(lambda=1.5, param=0.2,  max.M=2, M.dist="binom") # 0.764
1/(rho0^2) * ( gamma(1+2/kappa0) - gamma(1+1/kappa0)^2 )  # 1.538


### Try to emulate 1->2 density
# Mean sojourn time
ES.fn(lambda=1, param=6,   max.M=6, M.dist="tpois") # 5.410
ES.fn(lambda=1, param=6.5, max.M=6, M.dist="tpois") # 5.556 *
ES.fn(lambda=1, param=0.8, max.M=6, M.dist="binom") # 5.800 **
1/rho2 * gamma(1 + 1/kappa2)                        # 5.542

# Variance of sojourn time
VarS.fn(lambda=1, param=6,   max.M=6, M.dist="tpois")    # 7.294
VarS.fn(lambda=1, param=6.5, max.M=6, M.dist="tpois")    # 7.304 *
VarS.fn(lambda=1, param=0.8, max.M=6, M.dist="binom")    # 6.760 **
1/(rho2^2) * ( gamma(1+2/kappa2) - gamma(1+1/kappa2)^2 ) # 5.917


