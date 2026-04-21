## Social Conformity Task: Bayesian Models of Trustworthiness Updating

Study group 2  
(Barbora, Daniel, Mattis, Niels & Søren)

---

## Research Question

This assignment analyzes social conformity in trustworthiness judgments using Bayesian models of cognition. Participants first rate a face on a 1-8 trustworthiness scale, then see a group rating, and finally rate the face again. The core question is whether participants integrate their own prior judgment and the social signal in a way that is better captured by a proportional Bayesian agent or a weighted Bayesian agent.

The empirical analysis in this scaffold uses `cogsci_clean.csv`, the pre-pandemic student dataset.

---

## Minimal Data Description

Each row in the dataset corresponds to one participant-face observation.

- `ID`: participant identifier
- `FaceID`: face identifier
- `FirstRating`: initial trustworthiness judgment
- `GroupRating`: social information from others
- `SecondRating`: updated trustworthiness judgment after social exposure
- `Feedback`: signed difference between group and initial ratings
- `Change`: signed difference between second and initial ratings

To stay aligned with the Chapter 10 lecture notes, the 1-8 ratings are recoded as pseudo-counts on a 0-7 scale. This lets the models treat the participant's initial judgment and the group judgment as two evidence sources that jointly determine a posterior belief over the final rating.

Rows with `GroupRating = 0` are removed before fitting, because they do not correspond to a valid social rating on the observed 1-8 scale.

---

## Models

### Proportional Bayesian Agent (PBA)

#### Verbal description

The proportional Bayesian agent assumes that direct and social evidence are both used, but only their relative allocation matters. A parameter `p` determines how much of the total evidence budget is assigned to the participant's own initial rating, while `1 - p` determines how much is assigned to the group rating.

#### Formal description

Let `k_direct` denote the recoded first rating, `k_social` the recoded group rating, and `y` the recoded second rating, all on a 0-7 scale.

The posterior pseudo-counts are:

`alpha_post = 0.5 + p * k_direct + (1 - p) * k_social`

`beta_post = 0.5 + p * (7 - k_direct) + (1 - p) * (7 - k_social)`

The final rating is modeled as:

`y ~ BetaBinomial(7, alpha_post, beta_post)`

### Weighted Bayesian Agent (Reparameterized)

#### Verbal description

The reparameterized weighted Bayesian agent extends the proportional model by separating relative allocation from total evidence strength. The parameter `rho` controls the balance between direct and social information, while `kappa` controls the total amount of evidence used overall. This follows the reparameterization recommended in the course notes to avoid the identifiability problems of the raw `w_d` and `w_s` parameterization.

#### Formal description

`w_direct = rho * kappa`

`w_social = (1 - rho) * kappa`

`alpha_post = 0.5 + w_direct * k_direct + w_social * k_social`

`beta_post = 0.5 + w_direct * (7 - k_direct) + w_social * (7 - k_social)`

`y ~ BetaBinomial(7, alpha_post, beta_post)`

---

## Simulation Plan

The simulation stage mirrors the lecture pipeline.

1. Simulate synthetic observations from both PBA and reparameterized WBA over the same 0-7 pseudo-count task representation.
2. Fit both Stan models back to the synthetic datasets.
3. Inspect divergences, effective sample sizes, and other sampling diagnostics.
4. Compare prior and posterior predictive behavior.
5. Compare prior and posterior parameter distributions.
6. Evaluate model recovery through cross-fitting and predictive comparison.

Optional extensions from the notes include prior sensitivity, SBC, parameter recovery, and multilevel population models.

---

## Empirical Analysis Plan

For the `cogsci_clean.csv` dataset, the scaffold is set up to support:

1. Basic data exploration and recoding checks
2. Fitting `PBA` and reparameterized `WBA`
3. Sampling diagnostics, including divergences and E-BFMI
4. Predictive checks
5. Prior-posterior updates
6. Prior/likelihood sensitivity checks
7. Model comparison using PSIS-LOO
8. Interpretation of parameter estimates

---

## Results Template

### Data exploration

Summarize the distribution of first, group, and second ratings, and how strongly second ratings move toward the group.

### Model quality checks

Report divergences, E-BFMI, R-hat, ESS, predictive checks, and prior-posterior updates for both models.

### Model comparison

Report PSIS-LOO comparison between `PBA` and reparameterized `WBA`, and discuss whether the more flexible weighted model earns its extra complexity.

### Parameter interpretation

For `PBA`, discuss whether `p` suggests stronger reliance on self or group information. For `WBA`, interpret both `rho` and `kappa`: whether participants allocate more weight to direct or social evidence, and whether they use evidence conservatively or strongly overall.

---

## Current Status

This report file is intentionally structured as a writing scaffold. The empirical narrative, figures, and numerical results can be filled in directly once `assignment3_main.R`, `pba.stan`, `wba_reparam.stan`, and the outputs in `plots/` and `results/` have been generated.
