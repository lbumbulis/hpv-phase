
library(plotly)

plot.loglik <- function(mu2, max.M2, tau, pfail) {
  dat <- read.csv(paste0(
    "./results/param_grid_mu", mu2, "_maxM", max.M2, "_tau", tau, "_pfail", pfail, ".csv"
  ))
  
  plot.grid <- dat[order(dat[,1], dat[,2]),]
  
  xvals <- unique(plot.grid[,1])
  yvals <- unique(plot.grid[,2])
  grid.names <- names(plot.grid)
  
  theta <- set.lambda(theta=list(mu2=mu2, max.M2=max.M2), tau=tau, prop0=0.5, pfail=pfail)
  truth <- round(as.data.frame(as.list(unlist(theta)[
    which(names(theta) %in% gsub("log.", "", grid.names[1:2], fixed=T))
  ])), 2)
  truth <- merge(truth, plot.grid)
  truth <- truth[which(truth$loglik == max(truth$loglik)),]
  
  plotly.obj <- plot_ly(
    x = xvals,
    y = yvals,
    z = matrix(
      plot.grid$loglik,
      nrow = length(yvals),
      ncol = length(xvals),
      byrow = F
    ),
    type = "surface"
  )
  plotly.obj|>
    add_markers(
      data = truth,
      x = as.formula(paste0("~", grid.names[1])),
      y = as.formula(paste0("~", grid.names[2])),
      z = ~loglik,
      marker = list(color="red", size=5),
      name = ""
    ) %>% layout(
      xaxis = list(title=grid.names[1]),
      yaxis = list(title=grid.names[2])
    ) 
}

plot.loglik(mu2=0.5, max.M2=2, tau=4, pfail=0.5)
plot.loglik(mu2=0.5, max.M2=2, tau=4, pfail=0.7)
plot.loglik(mu2=0.5, max.M2=2, tau=4, pfail=0.9)

plot.loglik(mu2=6.5, max.M2=6, tau=4, pfail=0.5)
plot.loglik(mu2=6.5, max.M2=6, tau=4, pfail=0.7)
plot.loglik(mu2=6.5, max.M2=6, tau=4, pfail=0.9)

# There is a ridge, which is maybe slightly flatter for larger pfail.
# Patterns don't seem to vary much with (mu2, max.M2).


