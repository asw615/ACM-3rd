# Minimal simulation script aligned with the branch scaffold.
#
# 1. Scaling internally as 0:7
# 2. Builds an evidence grid task.
# 3. Models:
#    - Model 1: reparameterized WBA (rho, kappa)
#    - Model 2: PBA (p)
# 4. Exports files:
#    - dataset_1_task.csv
#    - dataset_2_model_choices.csv
SEED <- 123
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
  # Keep the task on the 0:7 scale internally, 
  #to avoid confusion with the 1:8 scale used in the user-facing CSV files.
build_evidence_grid <- function(total_direct = 7L,
                                total_social = 7L,
                                repeats_per_cell = 1L) {
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

task_counts <- build_evidence_grid(
  total_direct = 7L,
  total_social = 7L,
  repeats_per_cell = 1L
)

second_count_model_1 <- play_wba_reparam(task_counts)
second_count_model_2 <- play_pba(task_counts)

# Dataset 1:
  # export the fixed task in the observed 1:8 rating scale
dataset_1_task <- data.frame(
  trial_id = task_counts$trial_id,
  first_rating = task_counts$direct_count + 1L,
  group_rating = task_counts$social_count + 1L
)

# Dataset 2:
  # same task, plus each model's simulated second rating
dataset_2_model_choices <- data.frame(
  trial_id = task_counts$trial_id,
  first_rating = task_counts$direct_count + 1L,
  group_rating = task_counts$social_count + 1L,
  second_rating_model_1 = second_count_model_1 + 1L,
  second_rating_model_2 = second_count_model_2 + 1L
)

write.csv(
  dataset_1_task,
  file.path("outputs", "simulated", "dataset_1_task.csv"),row.names = FALSE
)

write.csv(
  dataset_2_model_choices,
  file.path("outputs", "simulated", "dataset_2_model_choices.csv"),row.names = FALSE
)

# SIMPLE PLOTTING
plot_data <- aggregate(
  cbind(second_rating_model_1, second_rating_model_2) ~ group_rating,
  data = dataset_2_model_choices,
  FUN = mean
)

plot(
  plot_data$group_rating,
  plot_data$group_rating,
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
  plot_data$group_rating,
  plot_data$second_rating_model_1,
  type = "b",
  pch = 19,
  col = "steelblue",
  lwd = 2
)

lines(
  plot_data$group_rating,
  plot_data$second_rating_model_2,
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
