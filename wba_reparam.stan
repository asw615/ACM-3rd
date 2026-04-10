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
  real<lower=0, upper=1> rho;              // relative allocation of total weight to direct evidence
  real<lower=0> kappa;                     // total amount of evidence weight used overall
}

model {
  real weight_direct = rho * kappa;                   // direct weight implied by the rho-kappa reparameterization
  real weight_social = (1.0 - rho) * kappa;           // social weight implied by the rho-kappa reparameterization

  target += beta_lpdf(rho | 2, 2);                    // prior centered on balanced direct vs social weighting
  target += lognormal_lpdf(kappa | log(2), 0.5);      // prior centered near the SBA-equivalent total weight of 2

  for (i in 1:N) {
    real alpha_post = 0.5 + weight_direct * blue1[i] + weight_social * blue2[i];   // posterior "success" pseudo-count
    real beta_post  = 0.5 + weight_direct * (total1[i] - blue1[i])                  // posterior "failure" pseudo-count from direct evidence
                           + weight_social * (total2[i] - blue2[i]);                // posterior "failure" pseudo-count from social evidence

    target += beta_binomial_lpmf(choice[i] | choice_total, alpha_post, beta_post);  // likelihood of the observed updated rating
  }
}

generated quantities {
  real lprior = beta_lpdf(rho | 2, 2) + lognormal_lpdf(kappa | log(2), 0.5); // joint log prior density
  real rho_prior = beta_rng(2, 2);                                             // prior draw for rho
  real kappa_prior = lognormal_rng(log(2), 0.5);                               // prior draw for kappa
  real wd_prior = rho_prior * kappa_prior;                                     // prior-implied direct weight
  real ws_prior = (1.0 - rho_prior) * kappa_prior;                             // prior-implied social weight
  vector[N] log_lik;                                                           // pointwise log-likelihood values for model comparison
  array[N] int prior_pred;                                                     // prior predictive simulated ratings
  array[N] int posterior_pred;                                                 // posterior predictive simulated ratings

  {
    real weight_direct = rho * kappa;                 // recompute fitted direct weight inside generated quantities
    real weight_social = (1.0 - rho) * kappa;         // recompute fitted social weight inside generated quantities

    for (i in 1:N) {
      real alpha_post = 0.5 + weight_direct * blue1[i] + weight_social * blue2[i];   // fitted posterior alpha for trial i
      real beta_post  = 0.5 + weight_direct * (total1[i] - blue1[i])                  // fitted posterior beta from direct evidence
                             + weight_social * (total2[i] - blue2[i]);                // fitted posterior beta from social evidence
      real alpha_prior_post = 0.5 + wd_prior * blue1[i] + ws_prior * blue2[i];        // prior-drawn alpha for prior PPC
      real beta_prior_post  = 0.5 + wd_prior * (total1[i] - blue1[i])                 // prior-drawn beta from direct evidence
                                   + ws_prior * (total2[i] - blue2[i]);               // prior-drawn beta from social evidence

      log_lik[i] = beta_binomial_lpmf(choice[i] | choice_total, alpha_post, beta_post);           // store pointwise fit to the observed data
      prior_pred[i] = beta_binomial_rng(choice_total, alpha_prior_post, beta_prior_post);          // simulate one rating from the prior predictive distribution
      posterior_pred[i] = beta_binomial_rng(choice_total, alpha_post, beta_post);                  // simulate one rating from the posterior predictive distribution
    }
  }
}
