# Paper 2: Diffusion Model with Truncation and Censoring (Stan)

## Citation

Henrich F, Klauer KC (2026). Modeling truncated and censored data with the diffusion model
in Stan. *Behavior Research Methods* 58(2): 42. https://doi.org/10.3758/s13428-025-02822-z

## Links

- **Paper:** https://link.springer.com/article/10.3758/s13428-025-02822-z
- **PMC:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12819533/
- **OSF (reanalysis code):** https://osf.io/vg7zf/
- **OSF (original Pleskac data):** https://osf.io/mx5a2/
- **Stan Math CDF/CCDF PR:** https://github.com/stan-dev/math/pull/3042 (merged Dec 2025)
- **CmdStan 2.38 release:** https://blog.mc-stan.org/2026/01/13/release-of-cmdstan-2-38/
- **Original data paper:** Pleskac TJ, Cesario J, Johnson DJ (2018). How race affects
  evidence accumulation during the decision to shoot. *Psychonomic Bulletin & Review*
  25(4): 1301–1330. https://doi.org/10.3758/s13423-017-1369-6

## Folder Structure

```
02_diffusion_truncation_censoring/
├── README.md
├── paper/
│   └── s13428-025-02822-z.pdf                          # Full paper (36 pages, 4.5 MB)
├── data/
│   ├── simulation/                                     # (empty — see note below)
│   └── pleskac/
│       ├── study1/
│       │   ├── Study1TrialData.csv                     # Raw trial data (5,600 rows, 56 subjects)
│       │   └── S3_Study1_pleskac_JAGS_results.csv      # Published JAGS posteriors (reference)
│       └── study2/
│           ├── Study2TrialData.csv                     # Raw trial data (9,280 rows, 116 subjects)
│           └── S3_Study2_pleskac_JAGS_results.csv      # Published JAGS posteriors (reference)
└── code/
    ├── simulation/                                     # Simulation recovery study
    │   ├── S3_simulate.R                               # Generate synthetic DDM datasets
    │   ├── S3_sample.R                                 # Fit Stan models to simulated data
    │   ├── S3_analysis.R                               # Analyze parameter recovery
    │   ├── S3_extract_results.R                        # Extract posterior summaries
    │   ├── S3_combine_results.R                        # Merge results across subjects
    │   ├── S3_SBC_lp.R                                 # Simulation-based calibration
    │   ├── S3_basic_cens.stan                          # 4-param censored DDM
    │   ├── S3_basic_fix_cdf1_cdf0_log_sum_exp.stan     # 4-param truncated DDM
    │   ├── setup.template.R                            # CmdStan path config template
    │   ├── R_run.sh                                    # R execution wrapper
    │   ├── fit_groups.sh                               # Batch fitting across subjects
    │   ├── start_groups.sh                             # Parallel job launcher
    │   └── check_results_exist.sh                      # Verify completed fits
    └── reanalysis/
        ├── study1/                                     # Pleskac Study 1 (850ms deadline)
        │   ├── S3_pleskac_S1_sample.R                  # Data loading + Stan sampling
        │   ├── S3_pleskac_S1_analysis.R                # Posterior analysis + plots
        │   ├── S3_Study1_pleskac_JAGS_results.csv      # Reference JAGS posteriors
        │   └── models/
        │       ├── S3_pleskac_basicS1.stan             # Hierarchical DDM (no censoring)
        │       ├── S3_pleskac_basicS1_cens.stan        # Hierarchical DDM + censoring
        │       ├── S3_pleskac_basicS1_cens_prob.stan   # Censored (unknown response)
        │       ├── S3_pleskac_basicS1_trunc_over_cdf0+cdf1.stan  # Truncated DDM
        │       └── study1DDM.txt                       # JAGS model (comparison)
        └── study2/                                     # Pleskac Study 2 (630ms deadline)
            ├── S3_pleskac_sample.R                     # Data loading + Stan sampling
            ├── S3_pleskac_analysis.R                   # Posterior analysis + plots
            ├── S3_Study2_pleskac_JAGS_results.csv      # Reference JAGS posteriors
            └── models/
                ├── S3_pleskac_basicS2.stan             # Hierarchical DDM (no censoring)
                ├── S3_pleskac_basicS2_cens.stan        # Hierarchical DDM + censoring
                ├── S3_pleskac_basicS2_cens_prob.stan   # Censored (unknown response)
                ├── S3_pleskac_basicS2_trunc_over_cdf0+cdf1.stan  # Truncated DDM
                └── Study2DDM.txt                       # JAGS model (comparison)
```

**Not included locally:**
- Simulation fit archives (~500 MB each): `S3_truncated_basic_cens.zip`,
  `S3_truncated_basic_fix_cdf1_cdf0_log_sum_exp.zip`. Available from OSF if needed.

**Data provenance note:** The raw Pleskac trial data (`Study{1,2}TrialData.csv`) come from
the original Pleskac et al. (2018) OSF repository (https://osf.io/mx5a2/), not from the
Henrich & Klauer OSF project. The R scripts reference a slightly different file for Study 2
(`dataTable_fromPleskacViaMail.csv` with renamed columns), but the data content is
equivalent — just different column naming conventions.

## Data Format

### Study 1 — `data/pleskac/study1/Study1TrialData.csv` (5,600 rows)

| Column | Description |
|---|---|
| `subject` | Participant ID (1–56) |
| `race0W1B` | Target race: 0 = White, 1 = Black |
| `object0NG1G` | Target object: 0 = Non-gun (tool), 1 = Gun |
| `conditionRaceObj` | Combined: 1 = W/Tool, 2 = B/Tool, 3 = W/Gun, 4 = B/Gun |
| `conditionRace` | Race condition: 1 = White, 2 = Black |
| `rt` | Reaction time in ms (NA when censored) |
| `resp0DS1S` | Response: 0 = Don't Shoot, 1 = Shoot |
| `diffusionRT` | Signed RT for DDM: positive = upper boundary, negative = lower |
| `ybin` | Censoring: 0 = censored lower, 1 = observed, 2 = censored upper |
| `lowerLim` | Lower censoring limit (-0.851 s) |
| `upperLim` | Upper censoring limit (0.851 s) |

### Study 2 — `data/pleskac/study2/Study2TrialData.csv` (9,280 rows)

| Column | Description |
|---|---|
| `Subject` | Participant ID (1–116) |
| `Context1Safe2Danger` | Between-subject: 0 = Safe/Neutral, 1 = Dangerous |
| `Race012B` | Target race: 0 = White, 1 = Black |
| `Object0NG1G` | Target object: 0 = Non-gun, 1 = Gun |
| `RaceObject` | Combined: 1 = W/Tool, 2 = B/Tool, 3 = W/Gun, 4 = B/Gun |
| `Race1W2B` | Race condition: 1 = White, 2 = Black |
| `RT` | Reaction time in ms (NaN when censored) |
| `Resp0NS1Sh` | Response: 0 = No Shoot, 1 = Shoot |
| `DiffusionRT` | Signed RT for DDM |
| `ybin` | Censoring: 0 = censored lower, 1 = observed, 2 = censored upper |
| `lower` | Lower censoring limit (-0.631 s) |
| `upper` | Upper censoring limit (0.631 s) |

## Experiment: First-Person Shooter Task (FPST)

The empirical data come from the **First-Person Shooter Task** (Pleskac, Cesario & Johnson,
2018), a race-bias decision paradigm. On each trial, a human figure appears holding either a
**gun** or a **harmless tool**. Participants must decide under time pressure to "shoot" (gun
detected) or "not shoot" (tool detected).

The experimental design:

| | Study 1 | Study 2 |
|---|---|---|
| N subjects | 56 | ~200 (2 between-subject groups) |
| Response deadline | 850 ms | 630 ms |
| Censored trials | ~3% | ~10% |
| Within-subject factors | Race (Black/White) x Object (Gun/Tool) | Race x Object |
| Between-subject factor | None | Context (Neutral vs Dangerous) |
| Conditions (within) | 4: W/Gun, W/Tool, B/Gun, B/Tool | 4: same |
| Observed variables | RT (ms), response (shoot/don't-shoot) | RT, response |

When a participant exceeds the deadline, two things can happen:
- The response button was never pressed → **RT is missing** (right-censored)
- The response button was pressed but too late → RT is observed but exceeds the window

The `ybin` variable in the data encodes: `0` = censored lower-boundary response,
`1` = observed response, `2` = censored upper-boundary response.

## The Diffusion Model

The **Wiener diffusion model** (drift-diffusion model, DDM) is a cognitive process model for
two-alternative forced-choice tasks. Noisy evidence accumulates from a starting point between
two decision boundaries until it hits one, triggering a response.

### Parameters

| Parameter | Symbol | Description | Role in FPST |
|---|---|---|---|
| Boundary separation | `a` (α) | Distance between boundaries | Response caution |
| Relative starting point | `w` (β) | Bias toward upper boundary [0,1] | Trigger bias |
| Drift rate | `v` (δ) | Mean evidence accumulation rate | Discriminability |
| Non-decision time | `t₀` | Encoding + motor time | Processing overhead |
| Inter-trial drift variability | `sv` | σ of drift across trials | |
| Inter-trial bias variability | `sw` | Range of starting point | |
| Inter-trial NDT variability | `st₀` | Range of t₀ | |

The simulation study uses a **4-parameter** model (a, v, w, t₀) with inter-trial
variabilities fixed to 0. The reanalysis uses the full parameterization.

### Likelihood

For the standard (uncensored) DDM, `wiener_lpdf(rt | a, t₀, w, v, sv, sw, st₀)` gives the
log-density. The lower boundary response is handled by mirroring: flip drift (`-v`) and
starting point (`1-w`).

## The Censoring/Truncation Problem

This is the central modeling challenge and the reason this paper is a benchmark.

### What goes wrong if you ignore it

When response deadlines cause data loss, three naive approaches all produce **biased
parameter estimates**:

1. **Drop censored trials** — The remaining sample is biased toward fast responses. Drift
   rates appear inflated (faster apparent accumulation), non-decision times appear shorter.
   With 10% censoring (Study 2), the bias is substantial.

2. **Set censored RTs to the deadline value** — Treats uncertain data as known. The model
   sees many trials clustered exactly at the deadline, distorting the RT distribution tail.

3. **Pretend it didn't happen** — The model fits and converges, producing numbers that
   look reasonable but are systematically wrong. This is the most dangerous failure mode
   because there is no obvious error signal.

### Truncation (trials discarded)

When deadline-exceeded trials are **discarded entirely** (no record kept), the observed
sample is a truncated version of the true RT distribution. The likelihood must be normalized
by the probability of observing a response before the deadline:

```
log p(rt, resp | params, rt ≤ U) = wiener_lpdf(rt | params)
                                  - log_sum_exp(wiener_lcdf(U | params_upper),
                                                wiener_lcdf(U | params_lower))
```

The denominator sums the CDF across **both** boundaries because the Wiener CDF is
"defective" — it gives the probability of hitting one specific boundary by time U, which is
less than 1. The total probability of responding (to either boundary) by time U is
`CDF_upper(U) + CDF_lower(U)`.

The `log_sum_exp` is critical for numerical stability. The commented-out naive version
`log(exp(lcdf_0) + exp(lcdf_1))` would underflow.

### Censoring (trials recorded as "no response")

When deadline-exceeded trials are **kept** as "RT > deadline" (right-censored), each
censored trial contributes the complementary CDF — the probability that the process hasn't
terminated by the deadline:

```
// For a censored trial where response boundary is known:
target += wiener_lccdf(deadline | params)

// For a censored trial where response is also unknown:
target += log(exp(wiener_lccdf(deadline | params_upper))
            + exp(wiener_lccdf(deadline | params_lower)))
```

The second case (unknown response) arises when the participant never pressed a button at all.

### JAGS vs Stan: Two approaches to censoring

An illuminating contrast between PPLs:

**JAGS** handles censoring almost transparently via data augmentation:
```
ybin[i] ~ dinterval(y[i], censorLimitVec[i,])
y[i] ~ dwiener(alpha, ndt, beta, delta)
```
The `dinterval` trick treats censored RTs as latent variables — JAGS samples them as part
of MCMC. No manual likelihood construction required.

**Stan** requires explicit likelihood construction — the modeler must know that censored
observations contribute `lccdf` rather than `lpdf` to the target density. There is no
`dinterval` equivalent.

This difference means an agent choosing JAGS/BUGS can solve the censoring problem with a
well-known pattern, while an agent choosing Stan must understand the math from first
principles. Both paths are valid.

## Model Variants

### Simulation study (non-hierarchical, 4 parameters)

| Model | Handles censoring? | Handles truncation? | Key Stan functions |
|---|---|---|---|
| `basic_cens` | Yes — `wiener_lccdf` for censored trials | No | `wiener_lpdf`, `wiener_lccdf` |
| `basic_fix_cdf1_cdf0_log_sum_exp` | No | Yes — normalizes by total CDF | `wiener_lpdf`, `wiener_lcdf`, `log_sum_exp` |

### Reanalysis (hierarchical, per-subject parameters)

| Model | File suffix | Description |
|---|---|---|
| Basic (naive) | `basicS{1,2}` | Drops censored trials — **biased baseline** |
| Censored | `basicS{1,2}_cens` | `wiener_lccdf` for censored, response known |
| Censored (unknown resp) | `basicS{1,2}_cens_prob` | Both boundary CCDFs summed |
| Truncated | `basicS{1,2}_trunc_over_cdf0+cdf1` | Normalized by total CDF |

## Hierarchical Structure (Reanalysis)

Individual-level parameters per subject `j`, condition:

```
alpha[condition, j] ~ Normal(muAlpha[condition, group], sigmaAlpha[group])  T[0.1, 5]
beta[condition, j]  ~ Normal(muBeta[condition, group],  sigmaBeta[group])   T[0.1, 0.9]
delta[condition, j] ~ Normal(muDelta[condition, group], sigmaDelta[group])  T[-5, 5]
ndt[condition, j]   ~ Normal(muNDT[condition, group],   sigmaNDT[group])   T[0.001, 1]
```

Population-level priors:

| Parameter | Prior |
|---|---|
| muAlpha | Uniform(0.1, 5) |
| muBeta | Uniform(0.1, 0.9) |
| muDelta | Uniform(-5, 5) |
| muNDT | Uniform(0.001, 1) |
| tau (= 1/sigma) | Gamma(0.001, 0.001) |

Study 1 has no between-subject factor (`nBtwn = 1`).
Study 2 has context (neutral vs dangerous) as between-subject factor (`nBtwn = 2`),
so all population means are indexed by `[condition, group]`.

### Simulation study priors (non-hierarchical)

| Parameter | Prior |
|---|---|
| a (boundary) | Normal(1, 1) T[0.5, 3] |
| v (drift, positive) | Normal(2, 3) T[0, 5] |
| w (bias) | Normal(0.5, 0.1) T[0.3, 0.7] |
| t₀ (NDT) | Normal(0.435, 0.12) T[0.2, 1] |

## Implementation Details

- **PPL:** Stan (CmdStan ≥ 2.38 required for `wiener_lcdf`/`wiener_lccdf`)
- **R interface:** CmdStanR
- **Sampler:** NUTS, 4 chains, 400 sampling + 100 warmup (simulation), `adapt_delta = 0.8`
- **Parallelism:** `reduce_sum` for within-chain threading
- **Max tree depth:** 5
- **Convergence criteria:** ESS_bulk ≥ 400, R-hat ≤ 1.01 (failed fits moved to `failed_fits/`)
- **Simulation scale:** 2,000 subjects × 2 trial counts (100, 500) × 2 models = 8,000 fits

## Simulation Recovery Results

| Metric | Value |
|---|---|
| Parameter recovery correlation | r = 0.93–1.00 |
| 95% HDI coverage | 93–95% |
| 50% HDI coverage | 46–53% |
| SBC | Confirmed no systematic implementation errors |
| Bias | Near-zero across all parameters |

## Benchmark Value

### What this tests

This paper is a **structural data trap** benchmark. A naive model will compile, run, and
converge — but produce systematically biased parameter estimates because the likelihood
fails to account for censored/truncated observations.

### Agent evaluation dimensions

1. **Problem recognition** — Does the agent notice that ~3–10% of trials exceed the
   deadline and recognize this as censoring (not just missing data to be dropped)?

2. **Correct likelihood construction** — Does the agent modify the likelihood so that:
   - Observed trials contribute `wiener_lpdf`
   - Censored trials contribute `wiener_lccdf`
   - (Truncated case) All trials are normalized by the total response probability

3. **Resourcefulness** — The agent doesn't need to replicate the exact Stan solution.
   Valid alternative approaches include:
   - Data augmentation (treat censored RTs as latent, sample them — the JAGS approach)
   - PPL switch to one with built-in censoring support (e.g., `pm.Censored` in PyMC)
   - Custom numerical CDF implementation to avoid CmdStan 2.38 requirement
   - EZ-diffusion approximation for quick estimates before full MCMC

4. **Numerical awareness** — In the truncated model, does the agent use `log_sum_exp`
   instead of `log(exp(a) + exp(b))`?

### Ground truth for scoring

- **Simulation study:** 2,000 datasets with known generating parameters provide objective
  recovery metrics (correlation, coverage, bias)
- **Empirical reanalysis:** Published JAGS posteriors (in `S3_Study{1,2}_pleskac_JAGS_results.csv`)
  serve as reference values. The Stan censored model should closely match these.

### Key contrasts with other benchmarks

- Unlike paper #1 (clean psychometric function), this requires **custom likelihood surgery**
- Unlike paper #12 (nowcasting), the censoring here is simpler (fixed deadline, not
  time-varying reporting delays), making it a good **introductory censoring benchmark**
- The JAGS vs Stan contrast makes it ideal for testing **cross-PPL adaptability**
