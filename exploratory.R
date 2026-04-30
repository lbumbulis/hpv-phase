
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
plot.sojourn.M <- function(lambda2, size.vals, mu.vals, max.x, target, out.col="density") {
  test.params <- expand.grid(x=seq(0.01, max.x, by=0.01), size=size.vals, mu=mu.vals)
  
  df <- test.params
  df <- cbind(df, t(sapply(1:nrow(df), function(r) {
    dM <- dnbinom(M.vals-1, size=df$size[r], mu=df$mu[r])
    
    x <- df$x[r]
    dens <- dgamma(x, shape=M.vals, rate=lambda2)
    surv <- pgamma(x, shape=M.vals, rate=lambda2, lower.tail=F)
    
    t(dM) %*% cbind(dens, surv, dens/surv)
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
      title    = bquote(bold("Gamma" ~ .(out.col) * ": rate" ~ lambda[2]==.(lambda2) ~
                               ", shape M ~ NB(" * mu * ", " * phi * ")+1")),
      x        = "x",
      y        = paste0(y.fn, "(x)")
    ) +
    theme.component
}

# Try to emulate 1->0 density
size.vals <- c(1, 5, 50)
mu.vals <- c(0, 1, 2)
plot.sojourn.M(lambda2=1.7, size.vals, mu.vals, max.x=5, "d10")
plot.sojourn.M(lambda2=1.7, size.vals, mu.vals, max.x=5, "none", "hazard")
plot.sojourn.M(lambda2=1.7, size.vals, mu.vals, max.x=5, "none", "survival")

# Try to emulate 1->2 density
size.vals <- c(1, 5, 50)
mu.vals <- c(1, 2, 8)
plot.sojourn.M(lambda2=1.5, size.vals, mu.vals, max.x=10, "d12")
plot.sojourn.M(lambda2=1.5, size.vals, mu.vals, max.x=10, "none", "hazard")
plot.sojourn.M(lambda2=1.5, size.vals, mu.vals, max.x=10, "none", "survival")



