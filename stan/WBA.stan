// WEIGHTED BAYESIAN AGENT (WBA)

data {
  int<lower=1> N; //
  array[N] int<lower=1, upper=8> FirstRating; //
  array[N] int<lower=1, upper=8> GroupRating; //
  array[N] int<lower=1, upper=8> SecondRating; //
}

parameters {
  real<lower=0, upper=1> rho;
  real<lower=0> kappa;
}

transformed parameters {
  real<lower=0> own_weighting = rho * kappa;
  real<lower=0> external_weighting = (1 - rho) * kappa;
}

model {
  // PRIOR
  target += beta_lpdf(rho | 2, 2);
  target += lognormal_lpdf(kappa | log(2), 0.5);

  // POSTERIOR LOOP
  for (i in 1:N) {
    real alpha_post = 0.5
      + own_weighting * (FirstRating[i] - 1)
      + external_weighting * (GroupRating[i] - 1);
    real beta_post = 0.5
      + own_weighting * (8 - FirstRating[i])
      + external_weighting * (8 - GroupRating[i]);
    target += beta_binomial_lpmf(SecondRating[i] - 1 | 7, alpha_post, beta_post);
  }
}

generated quantities {
  real lprior;
  real rho_prior = beta_rng(2, 2);
  real kappa_prior = lognormal_rng(log(2), 0.5);
  real own_weighting_prior = rho_prior * kappa_prior;
  real external_weighting_prior = (1 - rho_prior) * kappa_prior;

  vector[N] log_lik;
  array[N] int prior_pred;
  array[N] int posterior_pred;

  lprior = beta_lpdf(rho | 2, 2) + lognormal_lpdf(kappa | log(2), 0.5);

  for (i in 1:N) {
    real alpha_post = 0.5
      + own_weighting * (FirstRating[i] - 1)
      + external_weighting * (GroupRating[i] - 1);
    real beta_post = 0.5
      + own_weighting * (8 - FirstRating[i])
      + external_weighting * (8 - GroupRating[i]);

    real alpha_prior_pred = 0.5
      + own_weighting_prior * (FirstRating[i] - 1)
      + external_weighting_prior * (GroupRating[i] - 1);
    real beta_prior_pred = 0.5
      + own_weighting_prior * (8 - FirstRating[i])
      + external_weighting_prior * (8 - GroupRating[i]);

    log_lik[i] = beta_binomial_lpmf(SecondRating[i] - 1 | 7, alpha_post, beta_post);
    prior_pred[i] = 1 + beta_binomial_rng(7, alpha_prior_pred, beta_prior_pred);
    posterior_pred[i] = 1 + beta_binomial_rng(7, alpha_post, beta_post);
  }
}
