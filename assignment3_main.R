set.seed(2026)

suppressPackageStartupMessages({
  if (requireNamespace("pacman", quietly = TRUE)) {
    pacman::p_load(cmdstanr, posterior, loo, ggplot2, dplyr, readr, tibble, purrr, tidyr)
  }
})

# 0. SETTINGS

empirical_path <- "data/cogsci_clean.csv"
run_sampling <- TRUE
run_recovery <- TRUE

dir.create("plots", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# 1. DATA PREPARATION

prepare_cogsci_data <- function(path = "data/cogsci_clean.csv") {
  dat <- read.csv(path, stringsAsFactors = FALSE)

  required_cols <- c(
    "ID",
    "FaceID",
    "FirstRating",
    "GroupRating",
    "SecondRating",
    "Feedback",
    "Change"
  )

  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Some rows have GroupRating == 0 and missing feedback.
  # Those rows do not contain a usable social signal on the 1-8 scale,
  # so they are excluded before the Bayesian recoding.
  dat <- dat[!is.na(dat$GroupRating) & dat$GroupRating >= 1 & dat$GroupRating <= 8, , drop = FALSE]

  # Recode 1-8 ratings to 0-7 pseudo-counts so the models can follow
  # the beta-binomial setup from the course notes.
  dat$participant_id <- dat$ID
  dat$face_id <- dat$FaceID
  dat$direct_count <- as.integer(dat$FirstRating - 1L)
  dat$social_count <- as.integer(dat$GroupRating - 1L)
  dat$second_count <- as.integer(dat$SecondRating - 1L)
  dat$total_direct <- 7L
  dat$total_social <- 7L
  dat$choice_total <- 7L
  dat$evidence_gap <- as.integer(dat$GroupRating - dat$FirstRating)
  dat$rating_change <- as.integer(dat$SecondRating - dat$FirstRating)

  dat
}

make_stan_data <- function(dat) {
  list(
    N = nrow(dat),
    choice = as.integer(dat$second_count),
    blue1 = as.integer(dat$direct_count),
    blue2 = as.integer(dat$social_count),
    total1 = rep(7L, nrow(dat)),
    total2 = rep(7L, nrow(dat)),
    choice_total = 7L
  )
}

summarise_cogsci_data <- function(dat) {
  participant_level <- aggregate(
    cbind(FirstRating, GroupRating, SecondRating, rating_change) ~ participant_id,
    data = dat,
    FUN = mean
  )

  trial_level <- aggregate(
    cbind(FirstRating, GroupRating, SecondRating, rating_change) ~ face_id,
    data = dat,
    FUN = mean
  )

  list(
    participant_level = participant_level,
    trial_level = trial_level
  )
}

# 2. SIMULATION HELPERS

.draw_beta_binomial_count <- function(size, alpha, beta) {
  prob <- stats::rbeta(1, alpha, beta)
  stats::rbinom(1, size = size, prob = prob)
}

build_evidence_grid <- function(total_direct = 7L, total_social = 7L) {
  expand.grid(
    direct_count = 0:total_direct,
    social_count = 0:total_social
  ) |>
    transform(
      total_direct = total_direct,
      total_social = total_social
    )
}

simulate_pba_trials <- function(evidence_grid, p, seed = NULL, prior_alpha = 0.5, prior_beta = 0.5) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  weight_direct <- p
  weight_social <- 1 - p

  alpha_post <- prior_alpha +
    weight_direct * evidence_grid$direct_count +
    weight_social * evidence_grid$social_count

  beta_post <- prior_beta +
    weight_direct * (evidence_grid$total_direct - evidence_grid$direct_count) +
    weight_social * (evidence_grid$total_social - evidence_grid$social_count)

  simulated_count <- vapply(
    seq_len(nrow(evidence_grid)),
    function(i) .draw_beta_binomial_count(7L, alpha_post[i], beta_post[i]),
    integer(1)
  )

  cbind(
    evidence_grid,
    weight_direct = weight_direct,
    weight_social = weight_social,
    alpha_post = alpha_post,
    beta_post = beta_post,
    simulated_count = simulated_count,
    simulated_rating = simulated_count + 1L
  )
}

simulate_wba_reparam_trials <- function(evidence_grid, rho, kappa, seed = NULL, prior_alpha = 0.5, prior_beta = 0.5) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  weight_direct <- rho * kappa
  weight_social <- (1 - rho) * kappa

  alpha_post <- prior_alpha +
    weight_direct * evidence_grid$direct_count +
    weight_social * evidence_grid$social_count

  beta_post <- prior_beta +
    weight_direct * (evidence_grid$total_direct - evidence_grid$direct_count) +
    weight_social * (evidence_grid$total_social - evidence_grid$social_count)

  simulated_count <- vapply(
    seq_len(nrow(evidence_grid)),
    function(i) .draw_beta_binomial_count(7L, alpha_post[i], beta_post[i]),
    integer(1)
  )

  cbind(
    evidence_grid,
    rho = rho,
    kappa = kappa,
    weight_direct = weight_direct,
    weight_social = weight_social,
    alpha_post = alpha_post,
    beta_post = beta_post,
    simulated_count = simulated_count,
    simulated_rating = simulated_count + 1L
  )
}

# 3. MODEL FITTING

compile_assignment3_models <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("cmdstanr is required to compile the Stan models.")
  }

  list(
    pba = cmdstanr::cmdstan_model("pba.stan"),
    wba_reparam = cmdstanr::cmdstan_model("wba_reparam.stan")
  )
}

fit_empirical_models <- function(stan_data, seed = 2026) {
  models <- compile_assignment3_models()

  list(
    pba = models$pba$sample(
      data = stan_data,
      seed = seed,
      chains = 4,
      iter_warmup = 1000,
      iter_sampling = 1000,
      refresh = 250
    ),
    wba_reparam = models$wba_reparam$sample(
      data = stan_data,
      seed = seed,
      chains = 4,
      iter_warmup = 1000,
      iter_sampling = 1000,
      refresh = 250
    )
  )
}

# 4. RUN THE BASIC PIPELINE

empirical_data <- prepare_cogsci_data(empirical_path)
empirical_summary <- summarise_cogsci_data(empirical_data)
stan_data <- make_stan_data(empirical_data)

write.csv(empirical_summary$participant_level, "results/participant_summary.csv", row.names = FALSE)
write.csv(empirical_summary$trial_level, "results/trial_summary.csv", row.names = FALSE)

message("Prepared empirical data from: ", empirical_path)
message("Rows: ", nrow(empirical_data))
message("Participants: ", length(unique(empirical_data$participant_id)))
message("Faces: ", length(unique(empirical_data$face_id)))

if (run_recovery) {
  recovery_grid <- build_evidence_grid()

  recovery_stub <- list(
    pba_example = simulate_pba_trials(recovery_grid, p = 0.65, seed = 2026),
    wba_example = simulate_wba_reparam_trials(recovery_grid, rho = 0.75, kappa = 2.0, seed = 2026)
  )

  saveRDS(recovery_stub, "results/recovery_scaffold.rds")
  message("Saved recovery scaffold to results/recovery_scaffold.rds")
}

if (run_sampling) {
  fits <- fit_empirical_models(stan_data)
  saveRDS(fits, "results/empirical_fits.rds")
  message("Saved empirical fits to results/empirical_fits.rds")
}
