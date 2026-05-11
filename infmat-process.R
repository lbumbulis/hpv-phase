
library(tidyr)
library(dplyr)
library(ggrepel)

file.pattern <- paste0(
  "var_mu(\\d+\\.?\\d*)_maxM(\\d+\\.?\\d*)_prop(\\d+\\.?\\d*)_pfail(\\d+\\.?\\d*)\\.csv"
)

theme.component <- theme_minimal(base_size=12) +
  theme(
    plot.title    = element_text(face="bold"),
    strip.text    = element_text(face="bold", size=10),
    panel.border  = element_rect(fill=NA, linewidth=0.5),
    strip.background = element_rect(fill="grey95")
  )

results.files <- list.files("./results/", full.names=T)

plot.var <- function(max.M2, param) {
  target.files  <- results.files[which(
    grepl("var_", results.files, fixed=T) &
      grepl(paste0("maxM", max.M2), results.files, fixed=T)
  )]
  dat <- as.data.frame(data.table::rbindlist(lapply(target.files, function(f) {
    matches <- regexec(file.pattern, f)
    vals    <- regmatches(f, matches)[[1]]
    params  <- as.list(as.numeric(vals[2:5]))
    names(params) <- c("mu2", "max.M2", "prop0", "pfail")
    
    theta <- set.lambda(
      theta=list(mu2=params$mu2, max.M2=max.M2), tau=80, prop0=params$prop0, pfail=params$pfail
    )
    names(theta) <- paste0("true.", names(theta))
    cbind(data.frame(mu2=params$mu2, prop0=params$prop0, pfail=params$pfail),
          as.data.frame(theta), read.csv(f))
  })))
  
  plot.dat <- pivot_longer(dat, cols=starts_with("log."), names_to="param", values_to="var")
  plot.dat$truth <- with(plot.dat,
    log(ifelse(param=="log.lambda1", true.lambda1,
              ifelse(param=="log.lambda2", true.lambda2, true.mu2)))
  )
  
  plot.dat <- plot.dat |> mutate(
    rse1000      = abs(sqrt(var/1000) / truth),
    prop0_label = factor(paste0("pi[0]==", prop0), levels=paste0("pi[0]==", c(0.3, 0.5, 0.7))),
    pfail_label = factor(paste0("P(Z(tau)==2)==", pfail),
                         levels=paste0("P(Z(tau)==2)==", c(0.5, 0.7, 0.9)))
  )
  
  plot.dat <- plot.dat[grep(param, plot.dat$param, fixed=T),]
  if (param=="mu2") {
    ylims    <- c(0, 0.3)
    my.title <- bquote(log ~ mu[2] * "," ~ M[2]^"*" == .(max.M2))
    g        <- ggplot(plot.dat, aes(x=dt, y=rse1000, linetype=factor(mu2)))
    scale.colour.component <- geom_blank()
  } else {
    if (param=="lambda1") {
      ylims    <- c(0, 0.02)
      my.title <- bquote(log ~ lambda[1] * "," ~ M[2]^"*" == .(max.M2))
    } else if (param=="lambda2") {
      ylims    <- c(0, 0.07)
      my.title <- bquote(log ~ lambda[2] * "," ~ M[2]^"*" == .(max.M2))
    }
    g <- ggplot(plot.dat, aes(x=dt, y=rse1000, colour=truth, linetype=factor(mu2)))
    scale.colour.component <- scale_colour_continuous(name="true value")
  }
  
  g + geom_line(linewidth=0.8) +
    facet_grid(prop0_label ~ pfail_label, labeller=label_parsed) +
    scale_y_continuous(limits=ylims) +
    scale.colour.component +
    scale_linetype_manual(values=c("solid", "dashed", "dotted")) +
    guides(
      linetype = guide_legend(order=1),   # linetype legend above colour legend
      colour   = guide_colourbar(order=2)
    ) +
    labs(x="TIME BETWEEN VISITS", y="APPROX. RELATIVE SE (N = 1000)",
         linetype=expression(mu[2]), title=my.title) +
    theme.component +
    theme(plot.title = element_text(size=12))
}

plot.var(max.M2=5, param="mu2")
plot.var(max.M2=5, param="lambda1")
plot.var(max.M2=5, param="lambda2")

plot.var(max.M2=7, param="mu2")
plot.var(max.M2=7, param="lambda1")
plot.var(max.M2=7, param="lambda2")

