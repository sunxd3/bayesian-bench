# 06 — SARS-CoV-2 Household Transmission (Final-Size Model in Stan)

## Citation

van Boven M, van Dorp CH,"; et al. (2024).
**Estimation of introduction and transmission rates of SARS-CoV-2 in a prospective household study.**
*PLOS Computational Biology* 20(1): e1011832.

## Links

| Resource | URL |
|----------|-----|
| Paper | https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011832 |
| GitHub | https://github.com/mvboven/sars2-households |
| Zenodo | https://doi.org/10.5281/zenodo.10534386 |

## Folder structure

```
06_sars2_household_transmission/
├── README.md                              (this file)
├── paper/
│   └── journal.pcbi.1011832.pdf           (1.7 MB — full paper)
├── data/
│   ├── data_finalsize_stan_05072022.csv   (2.3 KB, 59 rows — household outbreak final sizes)
│   └── data_escape_stan_05072022.csv      (18 KB, 1,209 rows — person-level at-risk periods)
└── code/
    ├── REPO_README.md                     (original repo README)
    ├── households_github.stan             (30 KB, 594 lines — monolithic Stan model)
    └── households_github.R                (9 KB, 330 lines — CmdStanR driver + plotting)
```

## Study

The CoKids study: a prospective household cohort in the Netherlands (Aug 2020 -- Jul 2021). 307 households with children were followed. When any member tested positive, the entire household was monitored to track secondary infections. 59 households experienced at least one introduction.

The core question: what are the age-specific transmission rates within households, and how does the external introduction hazard vary over time with the Dutch pandemic waves?

## Data format

### `data_finalsize_stan_05072022.csv` — Household outbreaks (59 rows)

One row per infected household. Each household has a multi-type epidemic final size.

| Column | Description |
|--------|-------------|
| `household` | Household ID (1--307) |
| `j1`, `j2`, `j3` | New infections by type (1=children, 2=adolescents, 3=adults) |
| `a1`, `a2`, `a3` | Primary/co-primary cases by type |
| `n1`, `n2`, `n3` | Initially uninfected members by type |
| `c1`, `c2`, `c3` | Co-primary cases by type |
| `outbreak_start` | Day index of first case in household |
| `outbreak_end` | Day index of last case or follow-up end |
| `conditioning` | 0 or 1 — whether to condition on at least one infection |

### `data_escape_stan_05072022.csv` — Person-level at-risk periods (1,209 rows)

One row per person. Captures the time window during which each individual was at risk of external introduction.

| Column | Description |
|--------|-------------|
| `household` | Household ID |
| `time_start` | Day index when at-risk period begins |
| `time_end` | Day index when at-risk period ends (infection day or censoring) |
| `Type` | 1=child, 2=adolescent, 3=adult |
| `ext_infected` | 0=not infected externally, 1=primary/co-primary case |

## The modelling challenge: a recursive combinatorial likelihood

### Ball (1986) final-size distribution

Standard epidemic models describe *dynamics* (how infection propagates over time). This model instead describes the *final outcome*: given a household of known composition, what is the probability of observing exactly the final infection pattern?

The probability of `j` new infections given `a` primary cases and `n` susceptibles:

```
q(j | a, n, β) = C(n, j) · p(j | a, n, β)
```

where `p` is computed recursively:

```
p(ω | a, n) = b^{n-ω} · φ(β'(n-ω))^{a+ω}
              - Σ_{τ < ω} p(τ | a, n) · C(ω, τ) · φ(β'(n-ω))^{ω-τ}
```

Here `φ` is the Laplace transform of the infectious period distribution, `β` is the transmission rate matrix, and `b` is the vector of external escape probabilities. The recursion iterates over all multi-indices `ω ≤ j`, which for a household with 3 types means up to `(j₁+1)(j₂+1)(j₃+1)` terms.

### The Laplace transform trick

Rather than simulating the stochastic epidemic, the model integrates out the infectious period analytically via the Laplace transform. For a Gamma-distributed infectious period with shape=scale=50 (mean 1, concentrated around 0.75--1.25):

```
φ(x) = (1 + x/50)^{-50}
```

This is differentiable in the transmission rate parameters, enabling gradient-based sampling.

### Time-varying introduction hazard

The external introduction hazard is modelled as a penalised B-spline (P-spline) with 52 basis functions (50 knots, degree 3). The weights follow a random walk prior:

```
weights[i] = weights[i-1] + increment[i] · √(RWvar)
hazard(t) = exp(weights · B(t))
```

Children and adolescents have multiplicative factors relative to adults: `hazard_children(t) = ext_hazard_children · hazard_adults(t)`.

### Two-component likelihood

**Component 1 — Household outbreaks** (59 households):
```
log_lik[hh] = log(prob_infect_pattern(J, A, N, β, ext_esc))
```
where `ext_esc[k] = exp(-Σ_{t=start}^{end} hazard(t, k))` is the probability of *not* being introduced from outside during the outbreak window.

**Component 2 — Survival process** (1,209 persons):
For each person at risk from `t_start` to `t_end` of type `k`:
- If not infected: `log_lik += -Σ_{t} hazard(t, k)` (survived)
- If primary case: `log_lik += log(hazard(t_end, k)) - Σ_{t<t_end} hazard(t, k)` (infected on last day)

### Parameters

| Parameter | Count | Prior | Description |
|-----------|-------|-------|-------------|
| `rel_susceptibility` | 2 | Flat (>0) | Susceptibility of children, adolescents relative to adults |
| `infectivity` | 3 | Flat (>0) | Absolute infectivity by age type |
| `extra_trans` | 1 | Flat (>0) | Separate child-to-child transmission rate |
| `ext_hazard_weights` | 52 | N(-7.5, 2.5) first; N(0, 10) rest | B-spline weights (log scale) |
| `RWvar` | 1 | InvGamma(1, 0.0005) | Random walk smoothing variance |
| `ext_hazard_children` | 1 | Flat (>0) | Children's introduction hazard multiplier |
| `ext_hazard_adolescents` | 1 | Flat (>0) | Adolescents' introduction hazard multiplier |

**Total: ~61 parameters** (52 spline + 9 epidemiological).

Note: most epidemiological parameters have *only* positivity constraints and no informative priors — the posterior is entirely likelihood-driven for these.

### Model selection

8 model variants compared via:
- **LOO-IC** (Pareto-smoothed importance sampling LOO)
- **WBIC** (Watanabe 2013) — computed by setting `mode=1`, which tempers the likelihood by `1/log(num_households)`

Winning model: variable infectivity with separate child-to-child transmission (`extra_trans`).

## Implementation details

- **PPL**: Stan via CmdStanR (requires CmdStan 2.29+)
- **Sampler**: NUTS, 10 chains, 1000 iterations (1000 warmup), `thin=10`, `adapt_delta=0.99`, `max_treedepth=15`
- **Effective samples**: 1,000 total (100 per chain after thinning)
- **R dependencies**: tidyverse, cmdstanr, loo, rstan (for extraction), jsonlite, zoo, gridExtra, ggplot2
- **External dependency**: Dutch hospitalisation data from `data.rivm.nl` (for plotting only, not for model fitting)

The aggressive MCMC settings (`adapt_delta=0.99`, `max_treedepth=15`, `thin=10`) indicate a challenging posterior geometry — likely due to the spline funnel (`RWvar` → small means spline increments must also be small) and strong correlations between infectivity and susceptibility parameters.

## Benchmark value

### The structural challenge: recursive likelihood in Stan

This model's likelihood is not a standard distribution — it's a custom recursive computation (`prob_infect_pattern`) that iterates over multi-type infection outcomes. The function uses:
- Multi-index enumeration via `ravel_idx`/`unravel_idx`
- Integer binomial coefficients (`choose()`)
- Dynamic loop bounds dependent on data (household composition)
- The Laplace transform as a differentiable surrogate for stochastic simulation

An agent cannot simply call `neg_binomial_2_log()` or `bernoulli_logit()`. It must understand the mathematical structure of the final-size distribution and implement the recursion.

### Agent evaluation dimensions

1. **Mathematical modelling** — can the agent formulate the household final-size distribution from the epidemiological description? This requires understanding SIR/SIS dynamics, the Laplace transform approach, and the multi-type extension.

2. **Custom Stan functions** — the model requires writing `functions { }` block code with recursive combinatorial logic, multi-index array manipulation, and careful use of Stan's integer/real type system.

3. **Spline infrastructure** — implementing B-splines in Stan (either from scratch or adapting Kharratzadeh's case study), with penalised random walk priors.

4. **Two-component likelihood** — correctly combining the discrete final-size probability with the continuous survival process on the log scale.

5. **Posterior geometry awareness** — recognising that the RW1 + `RWvar` structure creates a Neal's funnel, and that the infectivity/susceptibility parameterisation requires care (the authors note proportional mixing has 2n-1 not 2n identifiable parameters).

6. **Model comparison** — implementing WBIC (tempered likelihood) alongside standard LOO-IC for model selection across variants.

### Comparison to other benchmark papers

| Dimension | This paper | Paper 02 (DDM) | Paper 05 (Rabies) |
|-----------|-----------|-----------------|-------------------|
| PPL | Stan (CmdStanR) | Stan (CmdStanR) | brms + RStan |
| Likelihood | Custom recursive (Ball 1986) | Wiener first-passage | neg_binomial_2_log |
| Key challenge | Combinatorial recursion | Censoring/truncation | Nonlinear power mean predictor |
| Splines | Yes (P-splines, 52 basis) | No | No |
| Funnel risk | Yes (RWvar + spline) | No | Yes (village RE) |
| Data size | Tiny (59 households) | Moderate (5,600 trials) | Large (87 villages × 252 months) |
| Parameters | ~61 | ~20 per condition | ~100 (village model) |
