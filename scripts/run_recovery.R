suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(loo)
})

# First-pass recovery script.
#
# This script:
# 1. Simulates datasets from a small grid of true PBA and WBA parameters.
# 2. Fits both Stan models to every simulated dataset.
# 3. Compares models with LOO.
# 4. Saves parameter recovery summaries.
# 5. Saves sampler diagnostics and funnel-relevant plots.
# 6. Flags whether a broader stress-test should be run next.

SEED <- 20260416L
RECOVERY_MODE <- Sys.getenv("RECOVERY_MODE", unset = "first_pass")
CHAINS <- 4L
ITER_WARMUP <- 500L
ITER_SAMPLING <- 500L
PARALLEL_CHAINS <- 4L
CHOICE_TOTAL <- 7L
DEFAULT_MAX_TREEDEPTH <- 10L
OUTPUT_DIR <- if (RECOVERY_MODE == "stress") file.path("outputs", "recovery_stress") else file.path("outputs", "recovery")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)


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
    trial_id = task$trial_id,
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
    trial_id = task$trial_id,
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


# -----------------------------
# Recovery design
# -----------------------------

if (RECOVERY_MODE == "stress") {
  pba_grid <- data.frame(
    true_model = "pba",
    setting_id = c("pba_p010", "pba_p025", "pba_p050", "pba_p075", "pba_p090"),
    p = c(0.10, 0.25, 0.50, 0.75, 0.90),
    rho = NA_real_,
    kappa = NA_real_,
    stringsAsFactors = FALSE
  )

  wba_grid <- expand.grid(
    rho = c(0.10, 0.25, 0.50, 0.75, 0.90),
    kappa = c(0.5, 2.0, 6.0)
  )

  wba_grid$true_model <- "wba"
  wba_grid$setting_id <- paste0(
    "wba_r",
    sprintf("%03d", as.integer(wba_grid$rho * 100)),
    "_k",
    gsub("\\.", "", format(wba_grid$kappa, trim = TRUE))
  )
  wba_grid$p <- NA_real_
  wba_grid <- wba_grid[, c("true_model", "setting_id", "p", "rho", "kappa")]
} else {
  # Option 1 chosen by the group:
  # - moderate parameter grid
  # - one simulated dataset per setting for the first pass
  pba_grid <- data.frame(
    true_model = "pba",
    setting_id = c("pba_p025", "pba_p050", "pba_p075"),
    p = c(0.25, 0.50, 0.75),
    rho = NA_real_,
    kappa = NA_real_,
    stringsAsFactors = FALSE
  )

  wba_grid <- data.frame(
    true_model = "wba",
    setting_id = c("wba_r025_k2", "wba_r050_k2", "wba_r075_k2"),
    p = NA_real_,
    rho = c(0.25, 0.50, 0.75),
    kappa = c(2.0, 2.0, 2.0),
    stringsAsFactors = FALSE
  )
}

recovery_grid <- rbind(pba_grid, wba_grid)


# -----------------------------
# Model compilation
# -----------------------------

models <- list(
  pba = cmdstan_model(file.path("stan", "pba.stan")),
  wba = cmdstan_model(file.path("stan", "wba_reparam.stan"))
)


# -----------------------------
# Diagnostics helpers
# -----------------------------

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

extract_draws_with_diagnostics <- function(fit, variables) {
  draws_df <- posterior::as_draws_df(fit$draws(variables = variables))
  sampler_array <- fit$sampler_diagnostics()
  n_iter <- dim(sampler_array)[1]
  n_chains <- dim(sampler_array)[2]

  diag_df <- expand.grid(
    .iteration = seq_len(n_iter),
    .chain = seq_len(n_chains)
  )

  diag_df <- diag_df[order(diag_df$.chain, diag_df$.iteration), , drop = FALSE]
  diag_df$divergent__ <- as.vector(sampler_array[, , "divergent__"])
  diag_df$treedepth__ <- as.vector(sampler_array[, , "treedepth__"])
  diag_df$energy__ <- as.vector(sampler_array[, , "energy__"])

  merge(draws_df, diag_df, by = c(".chain", ".iteration"), sort = FALSE)
}

save_diagnostic_plot <- function(fit, fitted_model, label, output_file) {
  if (fitted_model == "wba") {
    draws_df <- extract_draws_with_diagnostics(fit, c("rho", "kappa"))
    chain_colors <- c("black", "royalblue", "darkgreen", "purple")
    point_colors <- ifelse(draws_df$divergent__ == 1, "red", rgb(0, 0, 0, 0.25))

    pdf(output_file, width = 10, height = 8)
    par(mfrow = c(2, 2))

    plot(
      seq_len(nrow(draws_df)),
      draws_df$rho,
      col = chain_colors[draws_df$.chain],
      pch = 16,
      cex = 0.4,
      xlab = "Draw Index",
      ylab = "rho",
      main = paste(label, "- trace: rho")
    )

    plot(
      seq_len(nrow(draws_df)),
      draws_df$kappa,
      col = chain_colors[draws_df$.chain],
      pch = 16,
      cex = 0.4,
      xlab = "Draw Index",
      ylab = "kappa",
      main = paste(label, "- trace: kappa")
    )

    plot(
      draws_df$rho,
      draws_df$kappa,
      col = point_colors,
      pch = 16,
      cex = 0.6,
      xlab = "rho",
      ylab = "kappa",
      main = paste(label, "- rho vs kappa")
    )

    plot(
      log(draws_df$kappa),
      draws_df$rho,
      col = point_colors,
      pch = 16,
      cex = 0.6,
      xlab = "log(kappa)",
      ylab = "rho",
      main = paste(label, "- funnel check")
    )

    dev.off()
  }

  if (fitted_model == "pba") {
    draws_df <- extract_draws_with_diagnostics(fit, "p")
    chain_colors <- c("black", "royalblue", "darkgreen", "purple")

    pdf(output_file, width = 10, height = 8)
    par(mfrow = c(2, 2))

    plot(
      seq_len(nrow(draws_df)),
      draws_df$p,
      col = chain_colors[draws_df$.chain],
      pch = 16,
      cex = 0.4,
      xlab = "Draw Index",
      ylab = "p",
      main = paste(label, "- trace: p")
    )

    hist(
      draws_df$p,
      breaks = 30,
      col = "gray80",
      border = "white",
      xlab = "p",
      main = paste(label, "- posterior histogram")
    )

    plot(
      stats::density(draws_df$p),
      lwd = 2,
      xlab = "p",
      main = paste(label, "- posterior density")
    )

    plot(
      seq_len(nrow(draws_df)),
      draws_df$energy__,
      col = chain_colors[draws_df$.chain],
      pch = 16,
      cex = 0.4,
      xlab = "Draw Index",
      ylab = "energy__",
      main = paste(label, "- energy trace")
    )

    dev.off()
  }
}

fit_to_summary_row <- function(fit, fitted_model, true_model, setting_id, true_p, true_rho, true_kappa) {
  if (fitted_model == "pba") {
    param_summary <- fit$summary(variables = "p")
    param_names <- "p"
  } else {
    param_summary <- fit$summary(variables = c("rho", "kappa"))
    param_names <- c("rho", "kappa")
  }

  sampler_array <- fit$sampler_diagnostics()
  divergences <- sum(sampler_array[, , "divergent__"])
  treedepth_hits <- sum(sampler_array[, , "treedepth__"] >= DEFAULT_MAX_TREEDEPTH)
  ebfmi_values <- compute_ebfmi(fit)

  data.frame(
    true_model = true_model,
    fitted_model = fitted_model,
    setting_id = setting_id,
    true_p = true_p,
    true_rho = true_rho,
    true_kappa = true_kappa,
    divergences = divergences,
    treedepth_hits = treedepth_hits,
    max_rhat = max(param_summary$rhat, na.rm = TRUE),
    min_ess_bulk = min(param_summary$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(param_summary$ess_tail, na.rm = TRUE),
    min_ebfmi = min(ebfmi_values, na.rm = TRUE),
    unstable = divergences > 0 ||
      max(param_summary$rhat, na.rm = TRUE) > 1.01 ||
      min(param_summary$ess_bulk, na.rm = TRUE) < 200 ||
      min(ebfmi_values, na.rm = TRUE) < 0.3,
    stringsAsFactors = FALSE
  )
}

fit_to_parameter_recovery <- function(fit, fitted_model, true_model, setting_id, true_p, true_rho, true_kappa) {
  if (fitted_model == "pba") {
    param_summary <- fit$summary(variables = "p")

    return(data.frame(
      true_model = true_model,
      fitted_model = fitted_model,
      setting_id = setting_id,
      parameter = "p",
      true_value = true_p,
      mean = param_summary$mean,
      median = param_summary$median,
      q5 = param_summary$q5,
      q95 = param_summary$q95,
      recovered_in_90 = true_p >= param_summary$q5 && true_p <= param_summary$q95,
      stringsAsFactors = FALSE
    ))
  }

  param_summary <- fit$summary(variables = c("rho", "kappa"))
  true_lookup <- c(rho = true_rho, kappa = true_kappa)

  data.frame(
    true_model = true_model,
    fitted_model = fitted_model,
    setting_id = setting_id,
    parameter = param_summary$variable,
    true_value = unname(true_lookup[param_summary$variable]),
    mean = param_summary$mean,
    median = param_summary$median,
    q5 = param_summary$q5,
    q95 = param_summary$q95,
    recovered_in_90 = unname(true_lookup[param_summary$variable]) >= param_summary$q5 &
      unname(true_lookup[param_summary$variable]) <= param_summary$q95,
    stringsAsFactors = FALSE
  )
}


# -----------------------------
# Recovery run
# -----------------------------

task <- build_evidence_grid(repeats_per_cell = 1L)

loo_rows <- list()
diagnostic_rows <- list()
parameter_rows <- list()
dataset_rows <- list()

for (row_id in seq_len(nrow(recovery_grid))) {
  spec <- recovery_grid[row_id, ]
  sim_seed <- SEED + row_id

  if (spec$true_model == "pba") {
    sim_data <- simulate_from_pba(task, p = spec$p, seed = sim_seed)
  } else {
    sim_data <- simulate_from_wba(task, rho = spec$rho, kappa = spec$kappa, seed = sim_seed)
  }

  sim_data$true_model <- spec$true_model
  sim_data$setting_id <- spec$setting_id
  sim_data$first_rating <- sim_data$direct_count + 1L
  sim_data$group_rating <- sim_data$social_count + 1L
  sim_data$second_rating <- sim_data$choice + 1L
  dataset_rows[[row_id]] <- sim_data

  stan_data <- make_stan_data(sim_data)

  fits <- list(
    pba = models$pba$sample(
      data = stan_data,
      seed = sim_seed + 100L,
      chains = CHAINS,
      parallel_chains = PARALLEL_CHAINS,
      iter_warmup = ITER_WARMUP,
      iter_sampling = ITER_SAMPLING,
      refresh = 0,
      adapt_delta = 0.9
    ),
    wba = models$wba$sample(
      data = stan_data,
      seed = sim_seed + 200L,
      chains = CHAINS,
      parallel_chains = PARALLEL_CHAINS,
      iter_warmup = ITER_WARMUP,
      iter_sampling = ITER_SAMPLING,
      refresh = 0,
      adapt_delta = 0.9
    )
  )

  fit_labels <- c(
    pba = paste(spec$setting_id, "fit_pba", sep = "_"),
    wba = paste(spec$setting_id, "fit_wba", sep = "_")
  )

  for (fit_name in names(fits)) {
    fit <- fits[[fit_name]]

    diagnostic_rows[[length(diagnostic_rows) + 1L]] <- fit_to_summary_row(
      fit = fit,
      fitted_model = fit_name,
      true_model = spec$true_model,
      setting_id = spec$setting_id,
      true_p = spec$p,
      true_rho = spec$rho,
      true_kappa = spec$kappa
    )

    if (fit_name == spec$true_model) {
      parameter_rows[[length(parameter_rows) + 1L]] <- fit_to_parameter_recovery(
        fit = fit,
        fitted_model = fit_name,
        true_model = spec$true_model,
        setting_id = spec$setting_id,
        true_p = spec$p,
        true_rho = spec$rho,
        true_kappa = spec$kappa
      )
    }

    save_diagnostic_plot(
      fit = fit,
      fitted_model = fit_name,
      label = fit_labels[[fit_name]],
      output_file = file.path(OUTPUT_DIR, "plots", paste0(fit_labels[[fit_name]], ".pdf"))
    )
  }

  log_lik_pba <- fits$pba$draws(variables = "log_lik", format = "matrix")
  log_lik_wba <- fits$wba$draws(variables = "log_lik", format = "matrix")

  loo_pba <- loo::loo(log_lik_pba)
  loo_wba <- loo::loo(log_lik_wba)

  elpd_pba <- loo_pba$estimates["elpd_loo", "Estimate"]
  elpd_wba <- loo_wba$estimates["elpd_loo", "Estimate"]
  se_pba <- loo_pba$estimates["elpd_loo", "SE"]
  se_wba <- loo_wba$estimates["elpd_loo", "SE"]
  winner <- if (elpd_pba > elpd_wba) "pba" else "wba"

  loo_rows[[length(loo_rows) + 1L]] <- data.frame(
    true_model = spec$true_model,
    setting_id = spec$setting_id,
    true_p = spec$p,
    true_rho = spec$rho,
    true_kappa = spec$kappa,
    elpd_pba = elpd_pba,
    se_pba = se_pba,
    elpd_wba = elpd_wba,
    se_wba = se_wba,
    elpd_difference_wba_minus_pba = elpd_wba - elpd_pba,
    winner = winner,
    recovered_true_model = winner == spec$true_model,
    stringsAsFactors = FALSE
  )
}

recovery_datasets <- do.call(rbind, dataset_rows)
loo_summary <- do.call(rbind, loo_rows)
diagnostic_summary <- do.call(rbind, diagnostic_rows)
parameter_recovery <- do.call(rbind, parameter_rows)

warning_signs <- any(diagnostic_summary$unstable) || any(!loo_summary$recovered_true_model)

status_lines <- c(
  paste("seed:", SEED),
  paste("warning_signs:", warning_signs),
  paste("recovered_all_true_models:", all(loo_summary$recovered_true_model)),
  paste("any_unstable_fit:", any(diagnostic_summary$unstable))
)

writeLines(status_lines, con = file.path(OUTPUT_DIR, "recovery_status.txt"))

message("Recovery mode: ", RECOVERY_MODE)
message("Recovery complete.")
message("Warning signs present: ", warning_signs)
