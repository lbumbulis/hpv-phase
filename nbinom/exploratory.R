
library(ggplot2)
library(dplyr)
library(tidyr)

theme.component <- theme_minimal(base_size=12) +
  theme(
    plot.title    = element_text(face="bold"),
    strip.text    = element_text(face="bold", size=10),
    panel.border  = element_rect(fill=NA, linewidth=0.5),
    strip.background = element_rect(fill="grey95")
  )

lambda1 <- 0.15
rho0 <- 1/0.75
kappa0 <- 0.74
rho2 <- 1/6.25
kappa2 <- 2.43

###### NB distributions as a function of the parameters #######################
size.vals <- c(1, 2, 4)
mu.vals <- c(1, 2, 4)
NB.test.params <- expand.grid(x=0:20, size=size.vals, mu=mu.vals)

NB.plot.df <- NB.test.params |>
  mutate(
    density    = dnbinom(x, size=size, mu=mu),
    size_label = factor(paste0("phi==", size), levels=paste0("phi==", size.vals)),
    mu_label   = factor(paste0("mu==", mu), levels=paste0("mu==", mu.vals)),
  )

ggplot(NB.plot.df, aes(x=x, y=density, colour=factor(size), fill=factor(size))) +
  geom_col(width=1, alpha=0.6) +
  facet_grid(mu_label ~ size_label, labeller=label_parsed) +
  scale_colour_brewer(palette="Dark2", guide="none") +
  scale_fill_brewer(palette="Dark2", guide="none") +
  labs(
    title    = expression(bold("Negative Binomial")),
    subtitle = expression("Rows:" ~ mu ~ "   |   Columns:" ~ phi ~ "(dispersion parameter)"),
    x        = "x",
    y        = "P(X = x)"
  ) +
  theme.component

# Larger phi => less right-skewed
# Larger mu  => larger mean and variance

# To start, we will use M-2 ~ NB(phi=1, mu=1).
# We add 2 to the NB to get M since we always want M >= 2.



###### Sojourn time distributions as a function of M ##########################
M.vals <- 2:5
l2.vals <- c(2, 3, 4)
gamma.test.params <- expand.grid(x=seq(0.01, 5, by=0.01), M=M.vals, l2=l2.vals)

gamma.plot.df <- gamma.test.params |>
  mutate(
    density  = dgamma(x, shape=M, rate=l2),
    survival = pgamma(x, shape=M, rate=l2, lower.tail=F),
    hazard   = density / survival,
    M_label  = factor(paste0("M==", M), levels=paste0("M==", M.vals)),
    l2_label = factor(paste0("lambda[2]==", l2), levels=paste0("lambda[2]==", l2.vals)),
  )

# Plot sojourn distribution conditional on M
plot.sojourn.C <- function(out.col) {
  df <- plot.df
  col.idx <- which(names(df)==out.col)
  names(df)[col.idx] <- "y"
  
  if (out.col=="density") {
    y.fn <- "f"
  } else if (out.col=="hazard") {
    y.fn <- "h"
  } else if (out.col=="survival") {
    y.fn <- "1 - F"
  }
  
  ggplot(df, aes(x=x, y=y, colour=factor(M), fill=factor(M))) +
    geom_line(linewidth=0.8) +
    facet_grid(l2_label ~ M_label, labeller=label_parsed) +
    scale_colour_brewer(palette="Dark2", guide="none") +
    scale_fill_brewer(palette="Dark2", guide="none") +
    labs(
      title    = bquote(bold("Gamma" ~ .(out.col))),
      subtitle = expression("Rows: " * lambda[2] * " (rate)   |   Columns: " * M * " (shape)"),
      x        = "x",
      y        = paste0(y.fn, "(x)")
    ) +
    theme.component
}

plot.sojourn.C("density")
plot.sojourn.C("survival")
plot.sojourn.C("hazard")

# To start, we will use lambda2=2



###### Sojourn time distributions as a function of mu and phi #################
M.vals <- 1:100

# Plot sojourn distribution marginal over M
plot.sojourn.M <- function(lambda, size.vals, mu.vals, max.x, target, out.col="density") {
  test.params <- expand.grid(x=seq(0.01, max.x, by=0.01), size=size.vals, mu=mu.vals)
  
  df <- test.params
  df <- cbind(df, t(sapply(1:nrow(df), function(r) {
    M.pmf <- dnbinom(M.vals-1, size=df$size[r], mu=df$mu[r])
    
    x <- df$x[r]
    dens <- dgamma(x, shape=M.vals, rate=lambda)
    surv <- pgamma(x, shape=M.vals, rate=lambda, lower.tail=F)
    
    t(M.pmf) %*% cbind(dens, surv, dens/surv)
  })))
  names(df)[4:6] <- c("density", "survival", "hazard")
  df <- df |>
    mutate(
      d10        = dweibull(x, shape=kappa0, scale=1/rho0),
      d12        = dweibull(x, shape=kappa2, scale=1/rho2),
      size_label = factor(paste0("phi==", size), levels=paste0("phi==", size.vals)),
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
  
  ggplot(df, aes(x=x, y=y, colour=factor(size), fill=factor(size))) +
    target_component +
    geom_line(linewidth=0.8) +
    facet_grid(mu_label ~ size_label, labeller=label_parsed) +
    scale_colour_brewer(palette="Dark2", guide="none") +
    scale_fill_brewer(palette="Dark2", guide="none") +
    labs(
      title    = bquote(bold("Erlang" ~ .(out.col) * ": rate" ~ lambda==.(lambda) ~
                               ", shape M ~ NB(" * mu * ", " * phi * ")+1")),
      x        = "x",
      y        = paste0(y.fn, "(x)")
    ) +
    theme.component
}

# Try to emulate 1->0 density
size.vals <- c(1, 5, 50)
mu.vals <- c(0.1, 0.54, 1)
plot.sojourn.M(lambda=1.7, size.vals, mu.vals, max.x=5, "d10")
# save as params10.png (Width=500, Height=450)
plot.sojourn.M(lambda=1.7, size.vals, mu.vals, max.x=5, "none", "hazard")
plot.sojourn.M(lambda=1.7, size.vals, mu.vals, max.x=5, "none", "survival")

# Try to emulate 1->2 density
size.vals <- c(1, 5, 20, 50)
mu.vals <- c(1, 2, 5, 8)
plot.sojourn.M(lambda=1.6, size.vals, mu.vals, max.x=10, "d12")
# save as params12.png (Width=600, Height=500)
plot.sojourn.M(lambda=1.6, size.vals, mu.vals, max.x=10, "none", "hazard")
plot.sojourn.M(lambda=1.6, size.vals, mu.vals, max.x=10, "none", "survival")



###### Expected sojourn time as a function of mu and phi ######################
M.vals <- 1:100

ES.fn <- function(lambda, phi, mu, squared=F) {
  M.pmf <- dnbinom(M.vals-1, size=phi, mu=mu)
  
  if (squared) {
    integrand <- Vectorize(function(x) {
      M.pmf %*% dgamma(x, shape=M.vals, rate=lambda) * x^2
    })
  } else {
    integrand <- Vectorize(function(x) {
      M.pmf %*% dgamma(x, shape=M.vals, rate=lambda) * x
      # M.pmf %*% pgamma(x, shape=M.vals, rate=lambda, lower.tail=F) # same
    })
  }
  integrate(integrand, 0, Inf)$value
}

VarS.fn <- function(lambda, phi, mu) {
  ES.fn(lambda, phi, mu, squared=T) - ES.fn(lambda, phi, mu, squared=F)^2
}

# Try to emulate 1->0 mean sojourn time
ES.fn(lambda=1.7, phi=5, mu=0.1)  # 0.647
ES.fn(lambda=1.7, phi=5, mu=0.5)  # 0.882
ES.fn(lambda=1.7, phi=5, mu=0.54) # 0.906
ES.fn(lambda=1.7, phi=1, mu=0.54) # 0.906
1/rho0 * gamma(1 + 1/kappa0)      # 0.903

# Try to emulate 1->2 mean sojourn time
ES.fn(lambda=1.5, phi=20, mu=8) # 6.000
ES.fn(lambda=1.6, phi=20, mu=8) # 5.625; hard to estimate?
ES.fn(lambda=1.6, phi=20, mu=5) # 3.750
ES.fn(lambda=1.6, phi=20, mu=6) # 4.375
ES.fn(lambda=1.6, phi=20, mu=7) # 5.000
1/rho2 * gamma(1 + 1/kappa2)    # 5.542

# Note: For fixed lambda and mu, phi does not affect the mean

# Try to emulate 1->0 variance of sojourn time
VarS.fn(lambda=1.7, phi=5, mu=0.1)                       # 0.416
VarS.fn(lambda=1.7, phi=5, mu=0.5)                       # 0.709
VarS.fn(lambda=1.7, phi=5, mu=0.54)                      # 0.740
VarS.fn(lambda=1.7, phi=1, mu=0.54)                      # 0.821
1/(rho0^2) * ( gamma(1+2/kappa0) - gamma(1+1/kappa0)^2 ) # 1.538

# Try to emulate 1->2 variance of sojourn time
VarS.fn(lambda=1.5, phi=20, mu=8)                        # 8.978
VarS.fn(lambda=1.6, phi=20, mu=8)                        # 7.891
VarS.fn(lambda=1.6, phi=50, mu=8)                        # 7.141; hard to estimate phi
VarS.fn(lambda=1.6, phi=50, mu=5)                        # 4.492
VarS.fn(lambda=1.6, phi=50, mu=6)                        # 5.359
VarS.fn(lambda=1.6, phi=50, mu=7)                        # 6.242
1/(rho2^2) * ( gamma(1+2/kappa2) - gamma(1+1/kappa2)^2 ) # 5.917


