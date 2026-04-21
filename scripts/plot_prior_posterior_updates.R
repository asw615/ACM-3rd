suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Prior/posterior update and predictive-check script.
#
# This follows the course logic:
# 1. Fit the models to empirical data
# 2. Compare parameter priors vs posteriors
# 3. Compare prior predictive, posterior predictive, and observed ratings

SEED <- 20260421L
CHAINS <- 4L
ITER_WARMUP <- 1000L
ITER_SAMPLING <- 1000L
PARALLEL_CHAINS <- 4L
CHOICE_TOTAL <- 7L

dir.create(file.path("outputs", "model_checks"), recursive = TRUE, showWarnings = FALSE)


# -----------------------------
# Data preparation
# -----------------------------

prepare_empirical_data <- function(path = file.path("data", "cogsci_clean.csv")) {
  dat <- read.csv(path, stringsAsFactors = FALSE)

  dat <- dat[!is.na(dat$GroupRating) & dat$GroupRating >= 1 & dat$GroupRating <= 8, , drop = FALSE]
  dat <- dat[!is.na(dat$FirstRating) & !is.na(dat$SecondRating), , drop = FALSE]

  dat$direct_count <- as.integer(dat$FirstRating - 1L)
  dat$social_count <- as.integer(dat$GroupRating - 1L)
  dat$choice <- as.integer(dat$SecondRating - 1L)

  dat
}

make_stan_data <- function(dat) {
  list(
    N = nrow(dat),
    choice = as.integer(dat$choice),
    direct_count = as.integer(dat$direct_count),
    social_count = as.integer(dat$social_count),
    choice_total = CHOICE_TOTAL
  )
}


# -----------------------------
# Plot helpers
# -----------------------------

plot_density_update <- function(posterior_values, prior_values, parameter_name, output_file) {
  plot_df <- bind_rows(
    data.frame(value = posterior_values, distribution = "Posterior"),
    data.frame(value = prior_values, distribution = "Prior")
  )

  p <- ggplot(plot_df, aes(x = value, color = distribution, fill = distribution)) +
    geom_density(alpha = 0.2, linewidth = 1) +
    scale_color_manual(values = c("Posterior" = "blue", "Prior" = "red")) +
    scale_fill_manual(values = c("Posterior" = "blue", "Prior" = "red")) +
    labs(
      title = paste("Prior-Posterior Update:", parameter_name),
      x = parameter_name,
      y = "Density"
    ) +
    theme_classic()

  ggsave(output_file, plot = p, width = 8, height = 5)
}

extract_predictive_long <- function(fit, variable_name, observed_choice) {
  pred_df <- as.data.frame(fit$draws(variable_name, format = "draws_matrix"))
  names(pred_df) <- paste0("obs_", seq_len(ncol(pred_df)))

  pred_long <- pred_df |>
    mutate(draw_id = seq_len(nrow(pred_df))) |>
    pivot_longer(
      cols = starts_with("obs_"),
      names_to = "observation",
      values_to = "choice"
    ) |>
    mutate(choice = choice + 1L)

  observed_df <- data.frame(
    choice = observed_choice + 1L,
    distribution = "Observed"
  )

  pred_summary <- pred_long |>
    group_by(draw_id) |>
    count(choice, name = "n") |>
    group_by(choice) |>
    summarise(mean_n = mean(n), .groups = "drop")

  pred_summary
}

plot_predictive_check <- function(fit, observed_choice, variable_name, title_text, output_file) {
  pred_draws <- as.data.frame(fit$draws(variable_name, format = "draws_matrix"))
  names(pred_draws) <- paste0("obs_", seq_len(ncol(pred_draws)))

  pred_long <- pred_draws |>
    mutate(draw_id = seq_len(nrow(pred_draws))) |>
    pivot_longer(
      cols = starts_with("obs_"),
      names_to = "observation",
      values_to = "choice"
    ) |>
    mutate(choice = choice + 1L)

  pred_summary <- pred_long |>
    group_by(draw_id, choice) |>
    summarise(n = n(), .groups = "drop") |>
    group_by(choice) |>
    summarise(
      mean_n = mean(n),
      q10 = quantile(n, 0.10),
      q90 = quantile(n, 0.90),
      .groups = "drop"
    )

  observed_summary <- data.frame(choice = observed_choice + 1L) |>
    count(choice, name = "observed_n")

  plot_df <- left_join(pred_summary, observed_summary, by = "choice")
  plot_df$observed_n[is.na(plot_df$observed_n)] <- 0

  p <- ggplot(plot_df, aes(x = factor(choice))) +
    geom_col(aes(y = mean_n), fill = "gray80", color = "gray50") +
    geom_errorbar(aes(ymin = q10, ymax = q90), width = 0.2, color = "gray30") +
    geom_point(aes(y = observed_n), color = "blue", size = 2) +
    labs(
      title = title_text,
      x = "Second Rating",
      y = "Count"
    ) +
    theme_classic()

  ggsave(output_file, plot = p, width = 8, height = 5)
}


# -----------------------------
# Fit models
# -----------------------------

empirical_data <- prepare_empirical_data()
stan_data <- make_stan_data(empirical_data)

models <- list(
  pba = cmdstan_model(file.path("stan", "pba.stan")),
  wba = cmdstan_model(file.path("stan", "wba_reparam.stan"))
)

fits <- list(
  pba = models$pba$sample(
    data = stan_data,
    seed = SEED,
    chains = CHAINS,
    parallel_chains = PARALLEL_CHAINS,
    iter_warmup = ITER_WARMUP,
    iter_sampling = ITER_SAMPLING,
    refresh = 0,
    adapt_delta = 0.9
  ),
  wba = models$wba$sample(
    data = stan_data,
    seed = SEED + 1L,
    chains = CHAINS,
    parallel_chains = PARALLEL_CHAINS,
    iter_warmup = ITER_WARMUP,
    iter_sampling = ITER_SAMPLING,
    refresh = 0,
    adapt_delta = 0.9
  )
)


# -----------------------------
# Prior-posterior parameter updates
# -----------------------------

pba_draws <- as_draws_df(fits$pba$draws(c("p", "p_prior")))
plot_density_update(
  posterior_values = pba_draws$p,
  prior_values = pba_draws$p_prior,
  parameter_name = "p",
  output_file = file.path("outputs", "model_checks", "pba_prior_posterior_update_p.pdf")
)

wba_draws <- as_draws_df(fits$wba$draws(c("rho", "kappa", "rho_prior", "kappa_prior")))
plot_density_update(
  posterior_values = wba_draws$rho,
  prior_values = wba_draws$rho_prior,
  parameter_name = "rho",
  output_file = file.path("outputs", "model_checks", "wba_prior_posterior_update_rho.pdf")
)

plot_density_update(
  posterior_values = wba_draws$kappa,
  prior_values = wba_draws$kappa_prior,
  parameter_name = "kappa",
  output_file = file.path("outputs", "model_checks", "wba_prior_posterior_update_kappa.pdf")
)


# -----------------------------
# Prior and posterior predictive checks
# -----------------------------

plot_predictive_check(
  fit = fits$pba,
  observed_choice = empirical_data$choice,
  variable_name = "prior_pred",
  title_text = "PBA Prior Predictive Check",
  output_file = file.path("outputs", "model_checks", "pba_prior_predictive.pdf")
)

plot_predictive_check(
  fit = fits$pba,
  observed_choice = empirical_data$choice,
  variable_name = "posterior_pred",
  title_text = "PBA Posterior Predictive Check",
  output_file = file.path("outputs", "model_checks", "pba_posterior_predictive.pdf")
)

plot_predictive_check(
  fit = fits$wba,
  observed_choice = empirical_data$choice,
  variable_name = "prior_pred",
  title_text = "WBA Prior Predictive Check",
  output_file = file.path("outputs", "model_checks", "wba_prior_predictive.pdf")
)

plot_predictive_check(
  fit = fits$wba,
  observed_choice = empirical_data$choice,
  variable_name = "posterior_pred",
  title_text = "WBA Posterior Predictive Check",
  output_file = file.path("outputs", "model_checks", "wba_posterior_predictive.pdf")
)

message("Saved prior/posterior update plots to outputs/model_checks/")
