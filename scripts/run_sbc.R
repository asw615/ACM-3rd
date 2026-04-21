suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(ggplot2)
})

# Small SBC check.
#
# This is a lightweight systematic calibration pass:
# - draw true parameters from the model priors
# - simulate a dataset from those parameters
# - refit the same model
# - compute the rank of the true parameter inside the posterior draws
#
# If inference is well calibrated, the ranks should look roughly uniform.

SEED <- 20260421L
N_SIM <- 20L
CHAINS <- 4L
ITER_WARMUP <- 500L
ITER_SAMPLING <- 500L
PARALLEL_CHAINS <- 4L
CHOICE_TOTAL <- 7L

dir.create(file.path("outputs", "sbc"), recursive = TRUE, showWarnings = FALSE)


# -----------------------------
# Task and simulation helpers
# -----------------------------

draw_beta_binomial_count <- function(size, alpha, beta) {
  probability <- rbeta(1, shape1 = alpha, shape2 = beta)
  rbinom(1, size = size, prob = probability)
}

build_evidence_grid <- function(total_direct = 7L, total_social = 7L, repeats_per_cell = 1L) {
  base_grid <- expand.grid(
    direct_count = 0:total_direct,
    social_count = 0:total_social
  )

  grid <- base_grid[rep(seq_len(nrow(base_grid)), each = repeats_per_cell), , drop = FALSE]
  grid$trial_id <- seq_len(nrow(grid))
  grid$total_direct <- total_direct
  grid$total_social <- total_social
  grid
}

simulate_from_pba <- function(task, p, seed, prior_alpha = 0.5, prior_beta = 0.5) {
  set.seed(seed)

  alpha_post <- prior_alpha +
    p * task$direct_count +
    (1 - p) * task$social_count

  beta_post <- prior_beta +
    p * (task$total_direct - task$direct_count) +
    (1 - p) * (task$total_social - task$social_count)

  choice <- vapply(
    seq_len(nrow(task)),
    function(i) draw_beta_binomial_count(CHOICE_TOTAL, alpha_post[i], beta_post[i]),
    integer(1)
  )

  data.frame(
    direct_count = task$direct_count,
    social_count = task$social_count,
    choice = choice
  )
}

simulate_from_wba <- function(task, rho, kappa, seed, prior_alpha = 0.5, prior_beta = 0.5) {
  set.seed(seed)

  weight_direct <- rho * kappa
  weight_social <- (1 - rho) * kappa

  alpha_post <- prior_alpha +
    weight_direct * task$direct_count +
    weight_social * task$social_count

  beta_post <- prior_beta +
    weight_direct * (task$total_direct - task$direct_count) +
    weight_social * (task$total_social - task$social_count)

  choice <- vapply(
    seq_len(nrow(task)),
    function(i) draw_beta_binomial_count(CHOICE_TOTAL, alpha_post[i], beta_post[i]),
    integer(1)
  )

  data.frame(
    direct_count = task$direct_count,
    social_count = task$social_count,
    choice = choice
  )
}

make_stan_data <- function(sim_data) {
  list(
    N = nrow(sim_data),
    choice = as.integer(sim_data$choice),
    direct_count = as.integer(sim_data$direct_count),
    social_count = as.integer(sim_data$social_count),
    choice_total = CHOICE_TOTAL
  )
}

compute_ebfmi <- function(fit) {
  sampler_array <- fit$sampler_diagnostics()
  n_chains <- dim(sampler_array)[2]
  ebfmi_by_chain <- numeric(n_chains)

  for (chain_id in seq_len(n_chains)) {
    energy <- sampler_array[, chain_id, "energy__"]
    variance_energy <- stats::var(energy)

    if (is.na(variance_energy) || variance_energy == 0) {
      ebfmi_by_chain[chain_id] <- NA_real_
    } else {
      ebfmi_by_chain[chain_id] <- mean(diff(energy)^2) / variance_energy
    }
  }

  ebfmi_by_chain
}

compute_rank <- function(true_value, posterior_draws) {
  sum(posterior_draws < true_value)
}


# -----------------------------
# Model setup
# -----------------------------

set.seed(SEED)
task <- build_evidence_grid(repeats_per_cell = 1L)

models <- list(
  pba = cmdstan_model(file.path("stan", "pba.stan")),
  wba = cmdstan_model(file.path("stan", "wba_reparam.stan"))
)


# -----------------------------
# PBA SBC
# -----------------------------

pba_rows <- vector("list", N_SIM)

for (sim_id in seq_len(N_SIM)) {
  true_p <- rbeta(1, 2, 2)
  sim_data <- simulate_from_pba(task, p = true_p, seed = SEED + sim_id)
  stan_data <- make_stan_data(sim_data)

  fit <- models$pba$sample(
    data = stan_data,
    seed = SEED + 1000L + sim_id,
    chains = CHAINS,
    parallel_chains = PARALLEL_CHAINS,
    iter_warmup = ITER_WARMUP,
    iter_sampling = ITER_SAMPLING,
    refresh = 0,
    adapt_delta = 0.9
  )

  draws <- as_draws_matrix(fit$draws("p"))[, 1]
  sampler_array <- fit$sampler_diagnostics()
  param_summary <- fit$summary("p")

  pba_rows[[sim_id]] <- data.frame(
    model = "PBA",
    simulation = sim_id,
    parameter = "p",
    true_value = true_p,
    rank = compute_rank(true_p, draws),
    n_draws = length(draws),
    divergences = sum(sampler_array[, , "divergent__"]),
    max_rhat = param_summary$rhat,
    min_ess_bulk = param_summary$ess_bulk,
    min_ebfmi = min(compute_ebfmi(fit), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}


# -----------------------------
# WBA SBC
# -----------------------------

wba_rows <- vector("list", N_SIM)

for (sim_id in seq_len(N_SIM)) {
  true_rho <- rbeta(1, 2, 2)
  true_kappa <- rlnorm(1, log(2), 0.5)
  sim_data <- simulate_from_wba(task, rho = true_rho, kappa = true_kappa, seed = SEED + 2000L + sim_id)
  stan_data <- make_stan_data(sim_data)

  fit <- models$wba$sample(
    data = stan_data,
    seed = SEED + 3000L + sim_id,
    chains = CHAINS,
    parallel_chains = PARALLEL_CHAINS,
    iter_warmup = ITER_WARMUP,
    iter_sampling = ITER_SAMPLING,
    refresh = 0,
    adapt_delta = 0.9
  )

  draws <- as_draws_matrix(fit$draws(c("rho", "kappa")))
  sampler_array <- fit$sampler_diagnostics()
  param_summary <- fit$summary(c("rho", "kappa"))

  wba_rows[[sim_id]] <- data.frame(
    model = "WBA",
    simulation = sim_id,
    parameter = c("rho", "kappa"),
    true_value = c(true_rho, true_kappa),
    rank = c(
      compute_rank(true_rho, draws[, "rho"]),
      compute_rank(true_kappa, draws[, "kappa"])
    ),
    n_draws = nrow(draws),
    divergences = sum(sampler_array[, , "divergent__"]),
    max_rhat = rep(max(param_summary$rhat, na.rm = TRUE), 2),
    min_ess_bulk = rep(min(param_summary$ess_bulk, na.rm = TRUE), 2),
    min_ebfmi = rep(min(compute_ebfmi(fit), na.rm = TRUE), 2),
    stringsAsFactors = FALSE
  )
}

sbc_results <- do.call(rbind, c(pba_rows, wba_rows))


# -----------------------------
# Rank plots
# -----------------------------

rank_plot <- ggplot(sbc_results, aes(x = rank)) +
  geom_histogram(bins = 10, fill = "gray80", color = "gray30") +
  facet_grid(parameter ~ model, scales = "free_x") +
  labs(
    title = "Small SBC Rank Histograms",
    x = "Rank of true parameter within posterior draws",
    y = "Count"
  ) +
  theme_classic()

ggsave(
  file.path("outputs", "sbc", "sbc_rank_histograms.pdf"),
  plot = rank_plot,
  width = 10,
  height = 8
)


# -----------------------------
# Status summary
# -----------------------------

status_lines <- c(
  paste("seed:", SEED),
  paste("n_simulations_per_model:", N_SIM),
  paste("total_divergences:", sum(sbc_results$divergences)),
  paste("max_rhat:", max(sbc_results$max_rhat, na.rm = TRUE)),
  paste("min_ess_bulk:", min(sbc_results$min_ess_bulk, na.rm = TRUE)),
  paste("min_ebfmi:", min(sbc_results$min_ebfmi, na.rm = TRUE))
)

writeLines(status_lines, con = file.path("outputs", "sbc", "sbc_status.txt"))

message("Saved SBC plots to outputs/sbc/")
