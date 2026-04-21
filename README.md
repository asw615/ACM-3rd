# Advanced Cognitive Modelling — Assignment 3

This folder is structured around the third assignment in the course: fitting Bayesian models of cognition to real-world social conformity data.

The empirical dataset used for the main analysis is `data/cogsci_clean.csv`, which contains trustworthiness ratings from cognitive science students before and after exposure to social information. The two focal models are:

- `PBA` — Proportional Bayesian Agent
- `WBA (reparameterized)` — Weighted Bayesian Agent with `rho` and `kappa`

## Research Question

How do participants combine their own trustworthiness judgments with group information, and does a proportional or weighted Bayesian integration model better capture these rating updates?

## Modeling Adaptation

The lecture notes frame the agents as Beta-Binomial evidence integration models. To adapt them to the 1-8 trustworthiness ratings, this scaffold treats each rating as a pseudo-count on a 0-7 scale:

- `FirstRating - 1` becomes direct evidence
- `GroupRating - 1` becomes social evidence
- `SecondRating - 1` becomes the observed post-social rating

This keeps the inference pipeline close to the Chapter 10 notes while matching the actual task format.

Rows with `GroupRating == 0` are excluded during preprocessing, because they do not provide a valid 1-8 social rating and appear together with missing feedback in the student dataset.

## Folder Structure

- `assignment3_main.R` — the single main analysis script
- `pba.stan` — Stan model for the proportional Bayesian agent
- `wba_reparam.stan` — Stan model for the reparameterized weighted Bayesian agent
- `ACM_Portfolio_3.md` — report scaffold in the style of assignment 2
- `plots/` — generated figures
- `results/` — saved summaries, fitted objects, and model-comparison outputs

## Immediate Workflow

1. Run `Rscript assignment3_main.R` to prepare the empirical dataset and write exploratory summaries.
2. Turn on the sampling flags inside `assignment3_main.R` when you are ready to compile the Stan models and fit them.

## Main Files

- `assignment3_main.R`
- `pba.stan`
- `wba_reparam.stan`
