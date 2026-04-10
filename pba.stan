data {
  int<lower=1> N;                           // number of observations
  array[N] int<lower=0, upper=7> choice;   // observed second rating, recoded from 1-8 to 0-7
  array[N] int<lower=0, upper=7> blue1;    // direct evidence from the first rating, recoded to 0-7
  array[N] int<lower=0, upper=7> blue2;    // social evidence from the group rating, recoded to 0-7
  array[N] int<lower=1> total1;            // total pseudo-count scale for direct evidence, fixed at 7
  array[N] int<lower=1> total2;            // total pseudo-count scale for social evidence, fixed at 7
  int<lower=1> choice_total;               // total number of rating steps in the beta-binomial observation model
}

parameters {
  real<lower=0, upper=1> p;                // relative weight placed on direct evidence
}

model {
  target += beta_lpdf(p | 2, 2);           // weakly regularizing prior centered on equal weighting

  for (i in 1:N) {
    real alpha_post = 0.5 + p * blue1[i] + (1.0 - p) * blue2[i];   // posterior "success" pseudo-count
    real beta_post  = 0.5 + p * (total1[i] - blue1[i])             // posterior "failure" pseudo-count from direct evidence
                           + (1.0 - p) * (total2[i] - blue2[i]);   // posterior "failure" pseudo-count from social evidence

    target += beta_binomial_lpmf(choice[i] | choice_total, alpha_post, beta_post); // likelihood of the observed updated rating
  }
}

generated quantities {
  real lprior = beta_lpdf(p | 2, 2);       // log prior density for prior-sensitivity workflows
  real p_prior = beta_rng(2, 2);           // one draw from the prior for prior predictive simulation
  vector[N] log_lik;                       // pointwise log-likelihood values for model comparison
  array[N] int prior_pred;                 // prior predictive simulated ratings
  array[N] int posterior_pred;             // posterior predictive simulated ratings

  for (i in 1:N) {
    real alpha_post = 0.5 + p * blue1[i] + (1.0 - p) * blue2[i];         // fitted posterior alpha for trial i
    real beta_post  = 0.5 + p * (total1[i] - blue1[i])                   // fitted posterior beta from direct evidence
                           + (1.0 - p) * (total2[i] - blue2[i]);         // fitted posterior beta from social evidence
    real alpha_prior_post = 0.5 + p_prior * blue1[i] + (1.0 - p_prior) * blue2[i]; // prior-drawn alpha for prior PPC
    real beta_prior_post  = 0.5 + p_prior * (total1[i] - blue1[i])                 // prior-drawn beta from direct evidence
                                 + (1.0 - p_prior) * (total2[i] - blue2[i]);       // prior-drawn beta from social evidence

    log_lik[i] = beta_binomial_lpmf(choice[i] | choice_total, alpha_post, beta_post);            // store pointwise fit to the observed data
    prior_pred[i] = beta_binomial_rng(choice_total, alpha_prior_post, beta_prior_post);           // simulate one rating from the prior predictive distribution
    posterior_pred[i] = beta_binomial_rng(choice_total, alpha_post, beta_post);                   // simulate one rating from the posterior predictive distribution
  }
}
