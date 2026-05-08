
library(ggplot2)
library(dplyr)
library(tidyr)

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

###### Sojourn time distributions as a function of lambda and mu ##############

# Plot sojourn distribution marginal over M
plot.sojourn.M <- function(lambda.vals, mu.vals, max.M, max.x, target, out.col="density") {
  test.params <- expand.grid(x=seq(0.01, max.x, by=0.01), lambda=lambda.vals, mu=mu.vals)
  
  df <- test.params
  df <- cbind(df, t(sapply(1:nrow(df), function(r) {
    M.pmf <- dtpois.all(df$mu[r], max.M)
    
    x <- df$x[r]
    dens <- dgamma(x, shape=1:max.M, rate=df$lambda[r])
    surv <- pgamma(x, shape=1:max.M, rate=df$lambda[r], lower.tail=F)
    
    t(M.pmf) %*% cbind(dens, surv, dens/surv)
  })))
  names(df)[4:6] <- c("density", "survival", "hazard")
  df <- df |>
    mutate(
      d10        = dweibull(x, shape=kappa0, scale=1/rho0),
      d12        = dweibull(x, shape=kappa2, scale=1/rho2),
      lambda_label = factor(paste0("lambda==", lambda), levels=paste0("lambda==", lambda.vals)),
      mu_label   = factor(paste0("mu==", mu), levels=paste0("mu==", mu.vals))
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
  
  ggplot(df, aes(x=x, y=y, colour=factor(lambda), fill=factor(lambda))) +
    target_component +
    geom_line(linewidth=0.8) +
    facet_grid(mu_label ~ lambda_label, labeller=label_parsed) +
    scale_colour_brewer(palette="Dark2", guide="none") +
    scale_fill_brewer(palette="Dark2", guide="none") +
    labs(
      title = bquote(bold("Erlang" ~ .(out.col) * ": rate" ~ lambda *
                          ", shape M ~ tPois(" * mu * "," ~ tau==.(max.M-1) * ")+1")),
      x = "x",
      y = paste0(y.fn, "(x)")
    ) +
    theme.component
}

# Try to emulate 1->0 density
lambda.vals <- c(2, 3, 4, 5)
mu.vals <- c(0.1, 0.5, 1)
plot.sojourn.M(lambda.vals, mu.vals, max.M=3, max.x=5, "d10")
plot.sojourn.M(lambda.vals, mu.vals, max.M=4, max.x=5, "d10")
plot.sojourn.M(lambda.vals, mu.vals, max.M=5, max.x=5, "d10")
# save as params10.png (Width=500, Height=450)

# Try to emulate 1->2 density
lambda.vals <- c(0.5, 1, 1.5)
mu.vals <- c(4, 5, 6)
plot.sojourn.M(lambda.vals, mu.vals, max.M=7, max.x=12, "d12")
plot.sojourn.M(lambda.vals, mu.vals, max.M=8, max.x=12, "d12")
plot.sojourn.M(lambda.vals, mu.vals, max.M=10, max.x=12, "d12")
# save as params12.png (Width=600, Height=500)



###### Expected sojourn time as a function of mu and phi ######################
ES.fn <- function(lambda, mu, max.M, squared=F) {
  M.pmf <- dtpois.all(mu, max.M)
  
  if (squared) {
    integrand <- Vectorize(function(x) {
      M.pmf %*% dgamma(x, shape=1:max.M, rate=lambda) * x^2
    })
  } else {
    integrand <- Vectorize(function(x) {
      M.pmf %*% dgamma(x, shape=1:max.M, rate=lambda) * x
      # M.pmf %*% pgamma(x, shape=1:max.M, rate=lambda, lower.tail=F) # same
    })
  }
  integrate(integrand, 0, Inf)$value
}

VarS.fn <- function(lambda, mu, max.M) {
  ES.fn(lambda, mu, max.M, squared=T) - ES.fn(lambda, mu, max.M, squared=F)^2
}

# Try to emulate 1->0 mean sojourn time
ES.fn(lambda=2, mu=0.1, max.M=3) # 0.550
1/rho0 * gamma(1 + 1/kappa0)     # 0.903

# Try to emulate 1->2 mean sojourn time
ES.fn(lambda=1, mu=6, max.M=7) # 5.410
1/rho2 * gamma(1 + 1/kappa2)   # 5.542


# Try to emulate 1->0 variance of sojourn time
VarS.fn(lambda=2, mu=0.1, max.M=5)                       # 0.416
1/(rho0^2) * ( gamma(1+2/kappa0) - gamma(1+1/kappa0)^2 ) # 1.538

# Try to emulate 1->2 variance of sojourn time
VarS.fn(lambda=1, mu=6, max.M=7)                         # 7.294
1/(rho2^2) * ( gamma(1+2/kappa2) - gamma(1+1/kappa2)^2 ) # 5.917


