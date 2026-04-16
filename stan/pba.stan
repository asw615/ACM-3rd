
// PROPORTIONAL BAYESIAN AGENT (PBA) 
 //  - Uses proportional updating based on discrepancy
 //  - Includes variants (k and k-1 style updating)

data {
  int<lower=1> N;                       
  array[N] int<lower=0, upper=7> choice;   // Observed second rating, recoded from 1:8 to 0:7.
  array[N] int<lower=0, upper=7> direct_count;   // First rating on the 0:7 pseudo-count scale.
  array[N] int<lower=0, upper=7> social_count;   // Group rating on the 0:7 pseudo-count scale.
  int<lower=1> choice_total;     // Total number of steps on the scale, fixed at 7.
}

parameters {
  real<lower=0, upper=1> p;  // PBA parameter: weight placed on the participant's own first rating.
}


model {
  // PRIOR
  target += beta_lpdf(p | 2, 2);  // Weakly regularizing prior centered on balanced self vs group weighting.
  // POST. LOOP
  for (i in 1:N) {
    real alpha_post = 0.5
      + p * direct_count[i]
      + (1 - p) * social_count[i];   // Posterior alpha: weighted "success" evidence from self and group.

    real beta_post = 0.5
      + p * (choice_total - direct_count[i])
      + (1 - p) * (choice_total - social_count[i]); // Posterior beta: weighted "failure" evidence from self and group.
  // LIKELIHOOOD OF SECOND RATING
    target += beta_binomial_lpmf(choice[i] | choice_total, alpha_post, beta_post); 
    
    
  }
}

generated quantities {
  
  real lprior = beta_lpdf(p | 2, 2);   // Log prior density of the fitted parameter draw, useful for prior-posterior checks.
  real p_prior = beta_rng(2, 2);    // draws one fake value of p from the prior, for prior predictive simulation.
  vector[N] log_lik;            // Pointwise log-likelihood values for LOO.
  array[N] int prior_pred;       // Prior predictive second ratings on the 0:7 scale.
  array[N] int posterior_pred;    // Posterior predictive second ratings on the 0:7 scale.



  for (i in 1:N) {
    // This builds the fitted beta-binomial 
    real alpha_post = 0.5 
      + p * direct_count[i] // p controls the share of direct evidence.
                // 1 - p automatically becomes the share of social evidence
      + (1 - p) * social_count[i];   // This builds the fitted alpha parameter for trial i.
    real beta_post = 0.5
      + p * (choice_total - direct_count[i])
      + (1 - p) * (choice_total - social_count[i]);   // Recompute fitted posterior beta for this trial.
          // PBA is only learning relative weighting.
    
    // Prior-predictive version of the above = one random draw from the prior.
    real alpha_prior_post = 0.5
      + p_prior * direct_count[i]
      + (1 - p_prior) * social_count[i];   // Alpha implied by the prior draw of p.
    real beta_prior_post = 0.5
      + p_prior * (choice_total - direct_count[i])
      + (1 - p_prior) * (choice_total - social_count[i]); // Beta implied by the prior draw of p.

    log_lik[i] = beta_binomial_lpmf(choice[i] | choice_total, alpha_post, beta_post);   // Store trial-level fit for model comparison.
    prior_pred[i] = beta_binomial_rng(choice_total, alpha_prior_post, beta_prior_post);   // Simulate one rating from the prior predictive distribution.
    posterior_pred[i] = beta_binomial_rng(choice_total, alpha_post, beta_post);   // Simulate one rating from the posterior predictive distribution.
  
  }
}

