
library(ggplot2)
library(bayesplot)

##### Inspect results from a single model fit #################################
mcmc_trace(model_fit$draws(), pars=stan_params12$param)
model_fit$summary()

