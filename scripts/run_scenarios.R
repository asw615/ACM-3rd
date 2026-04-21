suppressPackageStartupMessages({
  library(ggplot2)
})

# Scenario comparison script aligned with Daniel's workbook.
#
# What this script does:
# 1. Builds one random 100-trial task.
# 2. Lets PBA and WBA play that same task.
# 3. Uses Daniel's literal scenario weights.
# 4. Saves side-by-side scenario plots.
#
# Important note:
# Daniel's PBA scenarios use literal own-weight values 1.5, 1.0, 0.5.
# In the PBA update rule, social weight is (1 - own_weighting).
# That means the self-focused PBA scenario uses a negative social weight.
# This is unusual, but it is implemented here on purpose to match his workbook.

SEED <- 20260421L
N_TRIALS <- 100L
ALPHA_PRIOR <- 4
BETA_PRIOR <- 3

dir.create(file.path("outputs", "scenarios"), recursive = TRUE, showWarnings = FALSE)

draw_beta_binomial_rating <- function(alpha_post, beta_post) {
  # If a scenario produces an invalid beta distribution, return NA.
  # This matches Daniel's practical workflow, where invalid rows are dropped later.
  if (alpha_post <= 0 || beta_post <= 0) {
    return(NA_integer_)
  }

  probability <- rbeta(1, shape1 = alpha_post, shape2 = beta_post)
  1L + rbinom(1, size = 7L, prob = probability)
}

build_random_task <- function(n_trials = N_TRIALS, seed = SEED) {
  set.seed(seed)

  data.frame(
    trial_id = seq_len(n_trials),
    first_rating = round(runif(n_trials, min = 1, max = 8)),
    group_rating = round(runif(n_trials, min = 1, max = 8))
  )
}

simulate_pba <- function(task, own_weighting, seed) {
  set.seed(seed)

  second_rating <- vapply(
    seq_len(nrow(task)),
    function(i) {
      alpha_post <- ALPHA_PRIOR +
        own_weighting * (task$first_rating[i] - 1) +
        (1 - own_weighting) * (task$group_rating[i] - 1)

      beta_post <- BETA_PRIOR +
        own_weighting * (8 - task$first_rating[i]) +
        (1 - own_weighting) * (8 - task$group_rating[i])

      draw_beta_binomial_rating(alpha_post, beta_post)
    },
    integer(1)
  )

  second_rating
}

simulate_wba <- function(task, own_weighting, external_weighting, seed) {
  set.seed(seed)

  second_rating <- vapply(
    seq_len(nrow(task)),
    function(i) {
      alpha_post <- ALPHA_PRIOR +
        own_weighting * (task$first_rating[i] - 1) +
        external_weighting * (task$group_rating[i] - 1)

      beta_post <- BETA_PRIOR +
        own_weighting * (8 - task$first_rating[i]) +
        external_weighting * (8 - task$group_rating[i])

      draw_beta_binomial_rating(alpha_post, beta_post)
    },
    integer(1)
  )

  second_rating
}

scenario_grid <- data.frame(
  scenario = c("Self-Focused", "Balanced", "Socially-Influenced"),
  pba_own_weighting = c(1.5, 1.0, 0.5),
  wba_own_weighting = c(1.5, 1.0, 0.5),
  wba_external_weighting = c(0.5, 1.0, 2.0),
  stringsAsFactors = FALSE
)

task <- build_random_task()
scenario_results <- vector("list", nrow(scenario_grid) * 2L)
row_index <- 1L

for (i in seq_len(nrow(scenario_grid))) {
  scenario_name <- scenario_grid$scenario[i]

  pba_second <- simulate_pba(
    task = task,
    own_weighting = scenario_grid$pba_own_weighting[i],
    seed = SEED + i
  )

  wba_second <- simulate_wba(
    task = task,
    own_weighting = scenario_grid$wba_own_weighting[i],
    external_weighting = scenario_grid$wba_external_weighting[i],
    seed = SEED + 100L + i
  )

  pba_df <- data.frame(
    scenario = scenario_name,
    model = "PBA",
    trial_id = task$trial_id,
    first_rating = task$first_rating,
    group_rating = task$group_rating,
    second_rating = pba_second
  )

  wba_df <- data.frame(
    scenario = scenario_name,
    model = "WBA",
    trial_id = task$trial_id,
    first_rating = task$first_rating,
    group_rating = task$group_rating,
    second_rating = wba_second
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
scenario_results <- scenario_results[!is.na(scenario_results$second_rating), , drop = FALSE]

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

pdf(file.path("outputs", "scenarios", "scenario_comparison_plots.pdf"), width = 10, height = 15)
print(plot_second)
print(plot_change)
print(plot_feedback)
dev.off()

message("Saved scenario plots to outputs/scenarios/scenario_comparison_plots.pdf")
