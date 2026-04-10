
// WEIGHTED BAYES AGENT, REPARAMETRIZED

data {
  int<lower=1> N;                          
  array[N] int<lower=0, upper=7> choice;  // Observed second rating, recoded from 1:8 to 0:7.
  array[N] int<lower=0, upper=7> direct_count; // First rating on the 0:7 pseudo-count scale.
  array[N] int<lower=0, upper=7> social_count; // Group rating on the 0:7 pseudo-count scale.
  int<lower=1> choice_total;     // Total number of steps on the scale, fixed at 7.
}

parameters {
  real<lower=0, upper=1> rho;  // Rho splits total evidence weight between self and group.
  real<lower=0> kappa;    // Kappa controls how much total evidence weight is used overall.
}

model {
  real weight_direct = rho * kappa;   // Direct-evidence weight implied by the rho-kappa reparameterization.
  real weight_social = (1 - rho) * kappa;  // Social-evidence weight implied by the rho-kappa reparameterization.
  
  // PRIORS
  target += beta_lpdf(rho | 2, 2);    // Prior centered on balanced allocation between self and group.
  target += lognormal_lpdf(kappa | log(2), 0.5); // Prior keeping total evidence weight positive and reasonably regularized.
  
  // POSTERIOR LOOP
  for (i in 1:N) {
    real alpha_post = 0.5
      + weight_direct * direct_count[i]
      + weight_social * social_count[i];   // Posterior alpha: weighted "success" evidence from self and group.

    real beta_post = 0.5
      + weight_direct * (choice_total - direct_count[i])
      + weight_social * (choice_total - social_count[i]);   // Posterior beta: weighted "failure" evidence from self and group.
      
    // LIKELIHOOD OF THE SECOND RATING
    target += beta_binomial_lpmf(choice[i] | choice_total, alpha_post, beta_post);  
  }
}

generated quantities {
  real lprior = beta_lpdf(rho | 2, 2) + lognormal_lpdf(kappa | log(2), 0.5);   // Joint log prior density.
  real rho_prior = beta_rng(2, 2);      // One prior draw of rho.
  real kappa_prior = lognormal_rng(log(2), 0.5); // One prior draw of kappa.
  real wd_prior = rho_prior * kappa_prior; // Direct weight implied by the prior draws.
  real ws_prior = (1 - rho_prior) * kappa_prior; // Social weight implied by the prior draws.
  vector[N] log_lik;              // Pointwise log-likelihood values for LOO.
  array[N] int prior_pred;          // Prior predictive second ratings on the 0:7 scale.
  array[N] int posterior_pred;      // Posterior predictive second ratings on the 0:7 scale.

  {
    // Convert the reparameterized form back into the actual direct and social weights used by the model.
    real weight_direct = rho * kappa;   
    real weight_social = (1 - rho) * kappa; // This is the fitted version, meaning it uses the current posterior draw of rho and kappa

    for (i in 1:N) {
      // This builds the fitted beta-binomial 
        // 1. Fitted posterior alpha for this trial.
      real alpha_post = 0.5 // Bigger alpha_post means more mass toward larger second ratings.
        + weight_direct * direct_count[i] // Participant's first rating on the 0:7 scale
        + weight_social * social_count[i]; // Group rating on the 0:7 scale
        // 2. Fitted posterior beta for this trial.
      real beta_post = 0.5 // beta tracks "low-rating evidence"
        + weight_direct * (choice_total - direct_count[i])
        + weight_social * (choice_total - social_count[i]);  // Fitted posterior beta for this trial.
        
      // Prior-predictive version of the above = one random draw from the prior.
      real alpha_prior_post = 0.5
        + wd_prior * direct_count[i]
        + ws_prior * social_count[i];  // Alpha implied by the prior draws of rho and kappa.
      real beta_prior_post = 0.5
        + wd_prior * (choice_total - direct_count[i])
        + ws_prior * (choice_total - social_count[i]); // Beta implied by the prior draws of rho and kappa.
        // "If parameters were drawn just from the prior, before seeing the data, what second-rating distribution would the model imply for this trial?"

      log_lik[i] = beta_binomial_lpmf(choice[i] | choice_total, alpha_post, beta_post);  // This stores the pointwise log-likelihood for observation i for model comparison.
        // It uses the fitted trial distribution
      prior_pred[i] = beta_binomial_rng(choice_total, alpha_prior_post, beta_prior_post);  // Simulate one rating from the prior predictive distribution.
        // useful for checking whether the priors generate sensible behavior before seeing data.
      posterior_pred[i] = beta_binomial_rng(choice_total, alpha_post, beta_post);   // Simulate one rating from the posterior predictive distribution.
        // useful for posterior predictive checks: does the fitted model reproduce the kinds of ratings we actually observe?

    }
  }
}
