# 09 — Epidemic Modelling with Survey Covariates (Turing.jl)

## Citation

Koher A, Brinch CN, Hartvig H, et al. (2023).
**Epidemic modelling of monitoring public behavior using surveys during pandemic-induced lockdowns.**
*Communications Medicine* 3:116.

## Links

| Resource | URL |
|----------|-----|
| Paper | https://www.nature.com/articles/s43856-023-00310-z |
| GitHub | https://github.com/andreaskoher/Covid19Survey |
| Zenodo | https://doi.org/10.5281/zenodo.7818793 |

## Folder structure

```
09_epidemic_turing/
├── README.md                              (this file)
├── paper/
│   └── s43856-023-00310-z.pdf             (1.9 MB — full paper)
├── data/
│   ├── data_with_predictors.csv           (1,225 rows × 172 cols — epi data + survey + mobility)
│   └── data_without_predictors.csv        (1,225 rows × 20 cols — epi data only)
└── code/
    ├── REPO_README.md                     (original repo README)
    ├── Project.toml                       (Julia package dependencies)
    ├── src/
    │   ├── Covid19Survey.jl               (38 lines — module definition)
    │   ├── models.jl                      (306 lines — two @model definitions)
    │   ├── data.jl                        (230 lines — data loading, serial interval, delay)
    │   ├── utils.jl                       (254 lines — KLogistic, RandomWalk, NegBin2, hpdi)
    │   └── postprocessing.jl              (630 lines — posterior processing)
    ├── inference_with_predictors.jl       (106 lines — MCMC driver with covariates)
    └── inference_without_predictors.jl    (104 lines — MCMC driver, baseline)
```

## Study

COVID-19 in Denmark (Aug 2020 -- Feb 2021), covering the second wave and lockdown period. 307 households tracked via the CoKids study. The core question: do self-reported contact surveys predict hospitalizations better than passively collected mobility data (Google, Apple, telco)?

Data covers 5 Danish regions over ~245 days, with daily hospitalization counts as the outcome and multiple candidate predictors from three data streams:
- **Survey contacts**: fraction reporting >N contacts with family, colleagues, friends, strangers
- **Google mobility**: retail, grocery, workplaces, residential, transit
- **Apple mobility**: driving, walking
- **Telco mobility**: aggregated cellphone movement

## Data format

### `data_with_predictors.csv` (1,225 rows × 172 columns)

One row per region × date. Key columns:

| Column | Description |
|--------|-------------|
| `date` | YYYY-MM-DD |
| `region` | 5 Danish regions |
| `hospit` | Daily new hospitalisations (outcome) |
| `cases`, `deaths` | Alternative outcomes |
| `population` | Regional population |
| `rw` | Random walk step index |
| `inference` | Boolean — in inference window |
| `family`, `colleagues`, `friends`, `strangers` | Survey mean contacts (relative to Dec 2020) |
| `total-above0` ... `total-above25` | Fraction reporting >N total contacts |
| `google_retail`, `google_workplaces`, ... | Google mobility change |
| `apple_driving`, `apple_walking`, `apple` | Apple mobility change |
| `telco` | Telco mobility change |
| `Q1_threat_to_society`, `Q3_avoid_contacts`, ... | Survey attitude questions |

### `data_without_predictors.csv` (1,225 rows × 20 columns)

Same epidemiological data without covariates. Used for the nonparametric baseline model.

## Two Turing.jl models

### Shared epidemic mechanics

Both models use the same generative process for infections and hospitalisations:

**1. Renewal equation:**
```
newly_infected[t] = R_t × Σ_{τ} newly_infected[t-τ] × serial_interval[τ]
```
Serial interval: `Gamma(mean=5.06, std=2.11)`, discretised into 15 bins.

**2. Infection → hospitalisation convolution:**
```
expected_hospit[t] = IHR × Σ_{τ} newly_infected[t-τ] × inf2hosp[τ]
```
Delay distribution: convolution of incubation `Gamma(mean=5.1, cv=0.86)` and hospitalisation delay `Weibull(0.845, 5.506)`, discretised into 40 bins.

**3. Observation model:**
```
hospit[t] ~ NegativeBinomial2(expected_hospit[t], φ)
```
Mean-variance parameterisation: `Var = μ + μ²/φ`.

### `parametricmodel` — R_t driven by covariates

```julia
R_t[i] = KLogistic(4)( latent_rt[step[i]] + Σ_j covariate[i,j] × effect[j] )
```

The reproduction number is a bounded (0, 4) transformation of a latent random walk *plus* linear covariate effects. The `KLogistic(K)` link function is `K / (1 + exp(-x))`.

**Hierarchical covariate effects:**
```julia
grouped_effect ~ filldist(Laplace(0, 0.2), num_covariates)   # national
effect_std     ~ GammaMeanStd(0.03, 0.02)                    # regional variation scale
effects_z      ~ filldist(MvNormal(num_covariates, 1), num_regions)
effects[m]     = effects_z[:,m] × effect_std + grouped_effect  # regional
```

### `nonparametricmodel` — R_t as pure random walk

```julia
R_t[i] = KLogistic(4)( latent_rt[step[i]] )
```

Same as above but without covariate effects. Used as baseline for model comparison.

### Parameters (parametric model)

| Parameter | Prior | Description |
|-----------|-------|-------------|
| `R0s[m]` | Truncated Normal(1, 0.1; 0, 4) per region | Initial R_t |
| `σ_rt` | Truncated Normal(0.3s, 0.02s; 0, 0.5s) | Random walk step size |
| `latent_rts_z[m]` | RandomWalk(n, σ_rt, μ_rt) per region | Latent R_t trajectory |
| `grouped_effect` | Laplace(0, 0.2) per covariate | National-level effect sizes |
| `effect_std` | GammaMeanStd(0.03, 0.02) | Regional effect deviation scale |
| `effects_z[m]` | MvNormal(K, 1) per region | Standardised regional deviations |
| `ys[m]` | Exponential(3 × init_infected) per region | Initial infection seed |
| `ihr` | Truncated Normal(0.028, 0.002; 0, 0.05) | Infection-hospitalisation rate |
| `φ` | GammaMeanStd(50, 20) | NegBin overdispersion |

### Custom distributions

| Name | Definition | Purpose |
|------|-----------|---------|
| `KLogistic(K)` | `K / (1 + exp(-x))` | Bounded link for R_t ∈ (0, K) |
| `RandomWalk(n, σ, μ₀, f)` | AR(1) with initial `Normal(μ₀, σ)` + increments `MvNormal(n-1, σ)` | Latent R_t process |
| `NegativeBinomial2(μ, φ)` | `NegBin(r=φ, p=1/(1+μ/φ))` | Mean-variance parameterisation |
| `GammaMeanStd(μ, σ)` | `Gamma(α=(μ/σ)², θ=σ²/μ)` | Convenient Gamma parameterisation |

## Implementation details

- **PPL**: Turing.jl (Julia 1.6)
- **AD backend**: ReverseDiff with caching
- **Sampler**: NUTS, target_accept = 0.99, max_treedepth = 5
- **Default**: 1000 warmup + 1000 post-warmup, 4 chains (`MCMCThreads()`)
- **Model comparison**: pointwise log-likelihoods → PSIS-LOO via the `:loglikelihood` context mode
- **Dependencies**: Turing, Distributions, ReverseDiff, MCMCChains, Bijectors, DrWatson, DataFrames

## Benchmark value

### Cross-PPL ecosystem coverage

This is the only Turing.jl (Julia) benchmark in the collection. It tests whether an agent can work in the Julia ecosystem — `@model` macros, custom distribution types, Bijectors, ReverseDiff AD, and Julia's multiple dispatch patterns.

### Agent evaluation dimensions

1. **Julia/Turing.jl proficiency** — can the agent write `@model` functions, define custom `ContinuousMultivariateDistribution` subtypes with `logpdf` methods, and use `filldist`/`arraydist` for vectorised sampling?

2. **Semi-mechanistic epidemic model** — the renewal equation with convolution-based delay distributions is not a standard regression. The agent must understand how latent infections drive observed hospitalisations through a known delay kernel.

3. **Custom distributions** — `RandomWalk`, `NegativeBinomial2`, `KLogistic` are all custom types with hand-written `logpdf` and `bijector` methods. An agent must implement these or find equivalents.

4. **Hierarchical covariate structure** — the `grouped_effect + region_deviation × effect_std` pattern implements partial pooling across Danish regions. The Laplace prior on effects encourages sparsity.

5. **Three context modes** — the same model serves inference, pointwise log-likelihood extraction, and prediction through a `cntxt` flag. This pattern is unusual and tests whether the agent understands the multi-use design.

6. **Model comparison** — comparing survey vs. mobility predictors via LOO-CV is the paper's central scientific contribution.

### Comparison to other benchmark papers

| Dimension | This paper | Paper 06 (SARS-CoV-2) | Paper 05 (Rabies) |
|-----------|-----------|----------------------|-------------------|
| PPL | Turing.jl | Stan (CmdStanR) | brms + RStan |
| Language | Julia | R | R |
| Epidemic model | Renewal equation | Household final-size | Neg. binomial GLMM |
| Observation | NegBin hospitalisations | Final-size distribution | NegBin case counts |
| Covariates | Survey + mobility | None (spline hazard) | Vaccination + spatial |
| Hierarchy | 5 regions, partial pooling | 307 households | 87 villages, random intercepts |
| Key challenge | Custom distributions + Julia | Recursive combinatorial likelihood | Nonlinear power mean |
