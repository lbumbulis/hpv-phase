
args <- commandArgs(trailingOnly=TRUE)
iter <- as.numeric(args[1])

library(cmdstanr)

source("source.R")

sim.seeds <- readRDS("sim_seeds_nsim1000.rds")

.Random.seed <- sim.seeds[[iter]]

pdat <- panelize(
  generate_NH(n=2000),
  visit.freq=1, tau=90
)

idx00 <- which(pdat$from_z==0 & pdat$to_z==0)

###########################################################
# Data processing
###########################################################
# Label sequence starts; a new 01x sequence begins when from_z=0 and to_z=1
pdat <- pdat %>%
  arrange(id, t_prev) %>%
  group_by(id) %>%
  mutate(
    seq_start = from_z == 0 & to_z == 1,
    seq_id    = cumsum(seq_start)
  ) %>%
  ungroup()

# Each sequence is the initial (0,1) row plus all following (1,*) rows, up to and
# including the first row where to_z changes back to 0 or follow-up ends with to_z=1.
sequences <- pdat %>%
  filter(seq_start | from_z==1) %>%
  group_by(id, seq_id) %>%
  summarise(
    dt_seq      = list(dt),
    t_seq       = list(cumsum(dt)), # cumulative times from sequence start
    final_state = last(to_z),
    .groups     = "drop"
  ) %>%
  mutate(key = sapply(dt_seq, paste, collapse="_"))

# Helper to build flattened arrays from a set of unique sequences
flatten_seqs <- function(unique_seqs, all_seqs) {
  J_vec     <- sapply(unique_seqs$t_seq, length)
  start_vec <- c(1, cumsum(J_vec[-length(J_vec)])+1)
  t_flat    <- unlist(unique_seqs$t_seq)
  idx       <- all_seqs %>%
    left_join(unique_seqs %>% select(key, k), by="key") %>%
    pull(k)
  list(
    n       = nrow(all_seqs),
    K       = nrow(unique_seqs),
    idx     = idx,
    J       = J_vec,
    t_flat  = t_flat,
    start   = start_vec
  )
}

# (0, 1, ..., 1, 0) sequences
seqs_010   <- sequences %>% filter(final_state==0)
unique_010 <- seqs_010 %>%
  distinct(key, .keep_all = TRUE) %>%
  mutate(k = row_number())
f010 <- flatten_seqs(unique_010, seqs_010)

# (0, 1, ..., 1) sequences
seqs_011   <- sequences %>% filter(final_state==1)
unique_011 <- seqs_011 %>%
  distinct(key, .keep_all = TRUE) %>%
  mutate(k = row_number())
f011 <- flatten_seqs(unique_011, seqs_011)

# (0, 1, ..., 1, 2) sequences
seqs_012   <- sequences %>% filter(final_state==2)
unique_012 <- seqs_012 %>%
  distinct(key, .keep_all = TRUE) %>%
  mutate(k = row_number())
f012 <- flatten_seqs(unique_012, seqs_012)

# Compute interval lengths from cumulative times for all positions in each flat array
deltas010 <- numeric(sum(f010$J))
for (k in seq_len(f010$K)) {
  st      <- f010$start[k]
  J       <- f010$J[k]
  times_k <- f010$t_flat[st:(st + J - 1)]
  deltas010[st:(st + J - 1)] <- diff(c(0, times_k))
}

deltas011 <- numeric(sum(f011$J))
for (k in seq_len(f011$K)) {
  st      <- f011$start[k]
  J       <- f011$J[k]
  times_k <- f011$t_flat[st:(st + J - 1)]
  deltas011[st:(st + J - 1)] <- diff(c(0, times_k))
}

deltas012 <- numeric(sum(f012$J))
for (k in seq_len(f012$K)) {
  st      <- f012$start[k]
  J       <- f012$J[k]
  times_k <- f012$t_flat[st:(st + J - 1)]
  deltas012[st:(st + J - 1)] <- diff(c(0, times_k))
}

unique_deltas <- sort(unique(c(deltas010, deltas011, deltas012)))

###########################################################

stan_data <- list(
  n00 = length(idx00),
  s00 = pdat$dt[idx00],
  
  n010       = f010$n,
  K010       = f010$K,
  idx010     = f010$idx,
  J010       = f010$J,
  t010_flat  = f010$t_flat,
  start010   = f010$start,
  
  n011       = f011$n,
  K011       = f011$K,
  idx011     = f011$idx,
  J011       = f011$J,
  t011_flat  = f011$t_flat,
  start011   = f011$start,
  
  n012       = f012$n,
  K012       = f012$K,
  idx012     = f012$idx,
  J012       = f012$J,
  t012_flat  = f012$t_flat,
  start012   = f012$start,
  
  nd            = length(unique_deltas),
  unique_deltas = unique_deltas,
  delta010_idx  = match(deltas010, unique_deltas),
  delta011_idx  = match(deltas011, unique_deltas),
  delta012_idx  = match(deltas012, unique_deltas),
  
  max_M0 = 5,
  max_M2 = 12
)

init_fn <- function() {
  list(
    log_lambda1  = runif(1, min=-3, max=0),
    
    log_mean_s10 = runif(1, min=-4, max=-1),
    log_mu0      = runif(1, min=-2, max=1),
    log_phi0     = runif(1, min=-1, max=1),
    
    log_mean_s12 = runif(1, min=1,  max=3),
    log_mu2      = runif(1, min=0,  max=2),
    log_phi2     = runif(1, min=1,  max=3.3)
  )
}

# set_cmdstan_path("C:/Users/lbumb/.cmdstan/cmdstan-2.38.0")
set_cmdstan_path("/home/lsbumbul/.cmdstan/cmdstan-2.38.0")

stan_model <- cmdstan_model("model.stan")

model_fit <- stan_model$sample(
  stan_data,
  init = init_fn,
  iter_sampling = 1000,
  iter_warmup = 1000,
  parallel_chains = 4,
  refresh = 5,
  save_warmup = TRUE
)
model_fit$save_object(file=paste0("./results/model_iter", iter, ".rds"))



