# Decision Log

## 2026-04-16

### Recovery setup
Chosen:
- Option 1 first-pass recovery
- Moderate parameter grid
- One simulated dataset per true parameter setting for the first pass

Not chosen:
- Option 2 task-faithful recovery
- Option 3 broad stress-aware recovery

Reason not chosen now:
- First pass is meant to check basic recovery and sampler stability before spending time on a larger grid.

### Formal recovery comparison
Chosen:
- Use `loo` for formal model comparison in recovery

Not chosen:
- Raw total log-likelihood only
- Diagnostics-only recovery without formal model comparison

### Escalation rule
Chosen:
- Run the broader stress-test only if the first-pass recovery shows warning signs

Warning signs defined as:
- Divergences
- High R-hat
- Low ESS
- Low E-BFMI
- Recovery failures under `loo`

### Outcome of first-pass recovery
Observed:
- Sampler diagnostics were clean enough to finish all fits
- One first-pass dataset was misclassified by `loo`
- WBA produced transient rejected-proposal warnings during sampling

Decision taken:
- Trigger the broader stress-test grid

Alternative not taken:
- Stop after the first pass and ignore the recovery warning
