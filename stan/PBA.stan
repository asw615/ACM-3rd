// PROPORTIONAL BAYESIAN AGENT (PBA)

data {
  int<lower=1> N; //
  array[N] int<lower=1, upper=8> FirstRating; //
  array[N] int<lower=1, upper=8> GroupRating; //
  array[N] int<lower=1, upper=8> SecondRating; //
}

parameters {
  real<lower=0, upper=1> p;
}

model {
  // PRIOR
  target += beta_lpdf(p | 2, 2);

  // POSTERIOR LOOP
  for (i in 1:N) {
    real alpha_post = 0.5
      + p * (FirstRating[i] - 1)
      + (1 - p) * (GroupRating[i] - 1);
    real beta_post = 0.5
      + p * (8 - FirstRating[i])
      + (1 - p) * (8 - GroupRating[i]);
    target += beta_binomial_lpmf(SecondRating[i] - 1 | 7, alpha_post, beta_post);
  }
}

generated quantities {
  real lprior;
  real p_prior = beta_rng(2, 2);
  vector[N] log_lik;
  array[N] int prior_pred;
  array[N] int posterior_pred;

  lprior = beta_lpdf(p | 2, 2);

  for (i in 1:N) {
    real alpha_post = 0.5
      + p * (FirstRating[i] - 1)
      + (1 - p) * (GroupRating[i] - 1);
    real beta_post = 0.5
      + p * (8 - FirstRating[i])
      + (1 - p) * (8 - GroupRating[i]);

    real alpha_prior_pred = 0.5
      + p_prior * (FirstRating[i] - 1)
      + (1 - p_prior) * (GroupRating[i] - 1);
    real beta_prior_pred = 0.5
      + p_prior * (8 - FirstRating[i])
      + (1 - p_prior) * (8 - GroupRating[i]);

    log_lik[i] = beta_binomial_lpmf(SecondRating[i] - 1 | 7, alpha_post, beta_post);
    prior_pred[i] = 1 + beta_binomial_rng(7, alpha_prior_pred, beta_prior_pred);
    posterior_pred[i] = 1 + beta_binomial_rng(7, alpha_post, beta_post);
  }
}
