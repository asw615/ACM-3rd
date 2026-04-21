# Minimal simulation script aligned with the shared branch scaffold.
#
# 1. Keeps the later fitting representation on the 0:7 pseudo-count scale.
# 2. Simulates one random 100-trial task, as in Daniel's workbook.
# 3. Models:
#    - Model 1: reparameterized WBA (rho, kappa)
#    - Model 2: PBA (p)
# 4. Exports files:
#    - dataset_1_task.csv
#    - dataset_2_model_choices.csv
SEED <- 123
N_TRIALS <- 100L
dir.create(file.path("outputs", "simulated"), recursive = TRUE, showWarnings = FALSE)
# -----------------------------
# Empirical-data helper
# -----------------------------

# This helper is for later fitting to the real data. Rules:
  # rows with GroupRating outside 1:8 are removed before fitting.
  # That means GroupRating == 0 / no-feedback rows are excluded.
prepare_for_fitting <- function(dat) {
  dat <- dat[!is.na(dat$GroupRating) & dat$GroupRating >= 1 & dat$GroupRating <= 8, , drop = FALSE]

  dat$direct_count <- as.integer(dat$FirstRating - 1L)
  dat$social_count <- as.integer(dat$GroupRating - 1L)
  dat$second_count <- as.integer(dat$SecondRating - 1L)

  dat
}

# -----------------------------
# Beta-binomial helper
# -----------------------------

# Draw one discrete response on the 0:7 scale:
# 1. sample a probability from Beta(alpha, beta)
# 2. sample a count from Binomial(7, probability)
draw_beta_binomial_count <- function(size, alpha, beta) {
  probability <- rbeta(1, shape1 = alpha, shape2 = beta)
  rbinom(1, size = size, prob = probability)
}


# -----------------------------
# Build the simulated task
# -----------------------------
# Daniel's workbook uses one random task with 100 trials.
# We sample observed ratings on the 1:8 scale, then convert them to 0:7
# internally so the Stan-ready representation stays consistent.
build_random_task <- function(n_trials = N_TRIALS, seed = SEED) {
  set.seed(seed)

  task <- data.frame(
    trial_id = seq_len(n_trials),
    first_rating = round(runif(n_trials, min = 1, max = 8)),
    group_rating = round(runif(n_trials, min = 1, max = 8))
  )

  task$direct_count <- task$first_rating - 1L
  task$social_count <- task$group_rating - 1L
  task$total_direct <- 7L
  task$total_social <- 7L

  task
}

# -----------------------------
# Model 1: reparameterized WBA
# -----------------------------

# WBA parameterization:
  # - rho controls how total weight is split between self and group
  # - kappa controls how much total evidence weight is used overall

# implied weights:
  # - direct weight = rho * kappa
  # - social weight = (1 - rho) * kappa
play_wba_reparam <- function(task,
                             rho = 0.75, kappa = 2.0,
                             seed = SEED + 1L,
                             prior_alpha = 0.5, prior_beta = 0.5) {
  set.seed(seed)

  weight_direct <- rho * kappa
  weight_social <- (1 - rho) * kappa

  alpha_post <- prior_alpha +
    weight_direct * task$direct_count +
    weight_social * task$social_count

  beta_post <- prior_beta +
    weight_direct * (task$total_direct - task$direct_count) +
    weight_social * (task$total_social - task$social_count)

  second_count <- vapply(
    seq_len(nrow(task)),
    function(i) draw_beta_binomial_count(7L, alpha_post[i], beta_post[i]),
    integer(1)
  )
  second_count
}
# -----------------------------
# Model 2: PBA
# -----------------------------

# PBA parameterization from the other branch:
  # - p controls how much relative weight goes to direct evidence
  # - (1 - p) goes to social evidence
play_pba <- function(task,
                     p = 0.65,
                     seed = SEED + 2L,
                     prior_alpha = 0.5,
                     prior_beta = 0.5) {
  set.seed(seed)

  alpha_post <- prior_alpha +
    p * task$direct_count +
    (1 - p) * task$social_count

  beta_post <- prior_beta +
    p * (task$total_direct - task$direct_count) +
    (1 - p) * (task$total_social - task$social_count)

  second_count <- vapply(
    seq_len(nrow(task)),
    function(i) draw_beta_binomial_count(7L, alpha_post[i], beta_post[i]),
    integer(1)
  )

  second_count
}

# -----------------------------
# RUN THE SIMULATION
# -----------------------------

set.seed(SEED)

task_counts <- build_random_task()

second_count_model_1 <- play_wba_reparam(task_counts)
second_count_model_2 <- play_pba(task_counts)

# Dataset 1:
# export the fixed task in the observed 1:8 rating scale
dataset_1_task <- data.frame(
  trial_id = task_counts$trial_id,
  first_rating = task_counts$first_rating,
  group_rating = task_counts$group_rating
)

# Match the empirical dataset definition:
# Feedback = GroupRating - FirstRating
dataset_1_task$feedback <- dataset_1_task$group_rating - dataset_1_task$first_rating

# Dataset 2:
# same task, plus each model's simulated second rating
dataset_2_model_choices <- data.frame(
  trial_id = task_counts$trial_id,
  first_rating = task_counts$first_rating,
  group_rating = task_counts$group_rating,
  second_rating_model_1 = second_count_model_1 + 1L,
  second_rating_model_2 = second_count_model_2 + 1L
)

dataset_2_model_choices$feedback <- dataset_2_model_choices$group_rating - dataset_2_model_choices$first_rating
dataset_2_model_choices$change_model_1 <- dataset_2_model_choices$second_rating_model_1 - dataset_2_model_choices$first_rating
dataset_2_model_choices$change_model_2 <- dataset_2_model_choices$second_rating_model_2 - dataset_2_model_choices$first_rating

write.csv(
  dataset_1_task,
  file.path("outputs", "simulated", "dataset_1_task.csv"),
  row.names = FALSE
)

write.csv(
  dataset_2_model_choices,
  file.path("outputs", "simulated", "dataset_2_model_choices.csv"),
  row.names = FALSE
)

# SIMPLE PLOTTING
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)
par(mfrow = c(3, 1))

plot_choices <- aggregate(
  cbind(second_rating_model_1, second_rating_model_2) ~ group_rating,
  data = dataset_2_model_choices,
  FUN = mean
)

plot(
  plot_choices$group_rating,
  plot_choices$group_rating,
  type = "b",
  pch = 16,
  lty = 2,
  col = "gray40",
  xlim = c(1, 8),
  ylim = c(1, 8),
  xlab = "Group Rating",
  ylab = "Average Rating",
  main = "Model Choices vs Group Rating"
)

lines(
  plot_choices$group_rating,
  plot_choices$second_rating_model_1,
  type = "b",
  pch = 19,
  col = "steelblue",
  lwd = 2
)

lines(
  plot_choices$group_rating,
  plot_choices$second_rating_model_2,
  type = "b",
  pch = 17,
  col = "firebrick",
  lwd = 2
)

legend(
  "topleft",
  legend = c("Group Rating", "Model 1: WBA", "Model 2: PBA"),
  col = c("gray40", "steelblue", "firebrick"),
  lty = c(2, 1, 1),
  pch = c(16, 19, 17),
  bty = "n"
)

plot_change <- aggregate(
  cbind(change_model_1, change_model_2) ~ feedback,
  data = dataset_2_model_choices,
  FUN = mean
)

plot(
  plot_change$feedback,
  plot_change$feedback,
  type = "b",
  pch = 16,
  lty = 2,
  col = "gray40",
  xlim = range(dataset_2_model_choices$feedback),
  ylim = range(
    c(
      dataset_2_model_choices$feedback,
      dataset_2_model_choices$change_model_1,
      dataset_2_model_choices$change_model_2
    )
  ),
  xlab = "Feedback",
  ylab = "Average Change From First Rating",
  main = "Model Change vs Feedback"
)

lines(
  plot_change$feedback,
  plot_change$change_model_1,
  type = "b",
  pch = 19,
  col = "steelblue",
  lwd = 2
)

lines(
  plot_change$feedback,
  plot_change$change_model_2,
  type = "b",
  pch = 17,
  col = "firebrick",
  lwd = 2
)

legend(
  "topleft",
  legend = c("Feedback", "Model 1: WBA", "Model 2: PBA"),
  col = c("gray40", "steelblue", "firebrick"),
  lty = c(2, 1, 1),
  pch = c(16, 19, 17),
  bty = "n"
)

# Feedback is discrete, so the density plot is more useful for the simulated
# change distributions than for feedback itself.
dens_1 <- density(dataset_2_model_choices$change_model_1)
dens_2 <- density(dataset_2_model_choices$change_model_2)

plot(
  dens_1,
  col = "steelblue",
  lwd = 2,
  main = "Density of Simulated Change",
  xlab = "Change From First Rating",
  ylim = range(c(dens_1$y, dens_2$y))
)

lines(
  dens_2,
  col = "firebrick",
  lwd = 2
)

legend(
  "topright",
  legend = c("Model 1: WBA", "Model 2: PBA"),
  col = c("steelblue", "firebrick"),
  lwd = 2,
  bty = "n"
)
