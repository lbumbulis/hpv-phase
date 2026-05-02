
library(ggplot2)
library(bayesplot)

##### Inspect results from a single model fit #################################
mcmc_trace(model_fit$draws(), pars=stan_params12$param)
model_fit$summary()

##### Process results from many replicates ####################################
library(dplyr)
library(tidyr)
model_type <- "012-recur"

results_files <- list.files(paste0("./results/"), full.names=T)
target_files <- results_files[grep(paste0("model", model_type, "_"), results_files, fixed=T)]
model_fit_list <- lapply(target_files, readRDS)
saveRDS(model_fit_list, paste0("model", model_type, "_fits.rds"))

all_draws_list <- lapply(model_fit_list, function(m) {
  posterior::as_draws_df(subset(m$draws(), variable="log_", regex=T))
})
all_draws <- bind_rows(all_draws_list, .id="iter") %>%
  mutate(iter=as.numeric(iter))
saveRDS(all_draws, file=paste0("model", model_type, "_draws.rds"))

all_summary <- as.data.frame(
  all_draws %>%
    pivot_longer(cols=starts_with("log_"), names_to="param", values_to="value") %>% # TODO
    group_by(iter, param) %>%
    summarise(
      est=mean(value), q025=quantile(value, 0.025), q975=quantile(value, 0.975),
      .groups="drop"
    )
)
saveRDS(all_summary, file=paste0("model", model_type, "_summary.rds"))

##### Inspect results from many replicates ####################################
coverage <- function(df, truth) {
  nparam <- nrow(truth)
  
  data.frame(
    param = truth$param,
    prob = sapply(1:nparam, function(jj) {
      p <- truth$param[jj]
      v <- truth$value[jj]
      temp <- df[which(df$param==p),]
      
      mean(temp$q025 <= v & temp$q975 >= v)
    })
  )
}

plot_est <- function(est_se, truth) {
  nparam <- nrow(truth)
  
  if (nparam==3) {
    old.par <- par(mfrow=c(1,3), mar=c(2,4,2,1))
    layout(matrix(c(1,2,3,rep(4,3)), nrow=2, byrow=T), heights=c(1, 0.1))
  } else if (nparam==7) {
    old.par <- par(mfrow=c(2,4), mar=c(2,4,2,1))
    layout(matrix(c(1:7, rep(8,7)), nrow=2, byrow=T), heights=c(1, 0.1))
  }
  
  for (jj in 1:nparam) {
    p <- truth$param[jj]
    temp <- est_se[which(est_se$param==p),]
    
    hist(temp$est, main=truth$param[jj], xlab="")
    abline(v=mean(temp$est), lwd=2, lty=2)
    abline(v=truth$val[jj], col="red", lwd=2, lty=2)
  }
  par(mar=c(0,0,0,0)+0.1)
  plot.new()
  legend("center", legend=c("estimate", "truth"), col=c("black","red"),
         lty=2, lwd=2, bty="n", horiz=T)
  
  par(old.par)
}

model_type <- "012-recur"
truth <- stan_params012_recur

all_summary <- readRDS(paste0("./results/model", model_type, "_summary.rds"))

coverage(all_summary, truth)
plot_est(all_summary, truth)

