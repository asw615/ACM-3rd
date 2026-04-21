suppressPackageStartupMessages({
  library(ggplot2)
})

# Scenario comparison script.
#
# This script does three things:
# 1. Formalizes the 3 scenarios from workbook.Rmd
# 2. Simulates both WBA and PBA on the same task for each scenario
# 3. Saves side-by-side plots for:
#    - SecondRating
#    - Change
#    - Feedback -> Change
#
# Important assumption:
# The workbook names the scenarios conceptually as:
# - Self-Focused
# - Balanced
# - Socially-Influenced
#
# Here they are mapped onto the current model parameterizations:
# - PBA uses p
# - WBA uses rho with a fixed kappa = 2
#
# This keeps the scenario labels aligned with the current Stan models.

SEED <- 20260421L
dir.create(file.path("outputs", "scenarios"), recursive = TRUE, showWarnings = FALSE)


# -----------------------------
# Task and simulation helpers
# -----------------------------

draw_beta_binomial_count <- function(size, alpha, beta) {
  probability <- rbeta(1, shape1 = alpha, shape2 = beta)
  rbinom(1, size = size, prob = probability)
}

build_evidence_grid <- function(total_direct = 7L, total_social = 7L, repeats_per_cell = 4L) {
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

  vapply(
    seq_len(nrow(task)),
    function(i) draw_beta_binomial_count(7L, alpha_post[i], beta_post[i]),
    integer(1)
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

  vapply(
    seq_len(nrow(task)),
    function(i) draw_beta_binomial_count(7L, alpha_post[i], beta_post[i]),
    integer(1)
  )
}


# -----------------------------
# Scenario definitions
# -----------------------------

scenario_grid <- data.frame(
  scenario = c("Self-Focused", "Balanced", "Socially-Influenced"),
  pba_p = c(0.75, 0.50, 0.25),
  wba_rho = c(0.75, 0.50, 0.25),
  wba_kappa = c(2.0, 2.0, 2.0),
  stringsAsFactors = FALSE
)


# -----------------------------
# Run scenarios
# -----------------------------

task <- build_evidence_grid(repeats_per_cell = 4L)
scenario_results <- vector("list", nrow(scenario_grid) * 2L)
row_index <- 1L

for (i in seq_len(nrow(scenario_grid))) {
  scenario_name <- scenario_grid$scenario[i]

  pba_choice <- simulate_from_pba(
    task = task,
    p = scenario_grid$pba_p[i],
    seed = SEED + i
  )

  wba_choice <- simulate_from_wba(
    task = task,
    rho = scenario_grid$wba_rho[i],
    kappa = scenario_grid$wba_kappa[i],
    seed = SEED + 100L + i
  )

  pba_df <- data.frame(
    scenario = scenario_name,
    model = "PBA",
    first_rating = task$direct_count + 1L,
    group_rating = task$social_count + 1L,
    second_rating = pba_choice + 1L
  )

  wba_df <- data.frame(
    scenario = scenario_name,
    model = "WBA",
    first_rating = task$direct_count + 1L,
    group_rating = task$social_count + 1L,
    second_rating = wba_choice + 1L
  )

  pba_df$feedback <- pba_df$group_rating - pba_df$first_rating
  pba_df$change <- pba_df$second_rating - pba_df$first_rating

  wba_df$feedback <- wba_df$group_rating - wba_df$first_rating
  wba_df$change <- wba_df$second_rating - wba_df$first_rating

  scenario_results[[row_index]] <- pba_df
  scenario_results[[row_index + 1L]] <- wba_df
  row_index <- row_index + 2L
}

scenario_results <- do.call(rbind, scenario_results)


# -----------------------------
# Plot helpers
# -----------------------------

plot_feedback_change <- aggregate(
  change ~ scenario + model + feedback,
  data = scenario_results,
  FUN = mean
)

plot_second <- ggplot(
  scenario_results,
  aes(x = second_rating, fill = model, color = model)
) +
  geom_density(alpha = 0.25) +
  facet_wrap(~scenario, ncol = 1) +
  scale_x_continuous(breaks = 1:8, limits = c(1, 8)) +
  labs(
    title = "Second Rating by Scenario",
    x = "Second Rating",
    y = "Density"
  ) +
  theme_minimal()

plot_change <- ggplot(
  scenario_results,
  aes(x = change, fill = model, color = model)
) +
  geom_density(alpha = 0.25) +
  facet_wrap(~scenario, ncol = 1) +
  labs(
    title = "Change From First Rating by Scenario",
    x = "Change",
    y = "Density"
  ) +
  theme_minimal()

plot_feedback <- ggplot(
  plot_feedback_change,
  aes(x = feedback, y = change, color = model)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~scenario, ncol = 1) +
  labs(
    title = "Average Change as a Function of Feedback",
    x = "Feedback",
    y = "Average Change"
  ) +
  theme_minimal()


# -----------------------------
# Save plots
# -----------------------------

pdf(file.path("outputs", "scenarios", "scenario_comparison_plots.pdf"), width = 10, height = 15)
print(plot_second)
print(plot_change)
print(plot_feedback)
dev.off()

message("Saved scenario plots to outputs/scenarios/scenario_comparison_plots.pdf")
