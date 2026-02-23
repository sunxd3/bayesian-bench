# 07 — Pain Learning and Statistical Prediction (Hierarchical RL/Kalman Filter in RStan)

## Citation

Onysk J, Varga S, Garro A, Allen M, Mancini F (2024).
**Statistical learning shapes pain perception and prediction independently of external cues.**
*eLife* 12:RP90634.

## Links

| Resource | URL |
|----------|-----|
| Paper | https://doi.org/10.7554/eLife.90634 |
| PMC | https://pmc.ncbi.nlm.nih.gov/articles/PMC11236420/ |
| GitHub | https://github.com/nox-lab/tsl_elife_24 |
| Zenodo | https://doi.org/10.5281/zenodo.11394627 |

## Folder structure

```
07_pain_learning/
├── README.md                          (this file)
├── paper/
│   └── elife-90634.pdf                (9.2 MB — full paper, v1)
├── data/
│   ├── a_coeffs.csv                   (675 B — per-subject linear transform slope)
│   ├── b_coeffs.csv                   (702 B — per-subject linear transform intercept)
│   └── pptIds.csv                     (112 B — participant IDs, 30 subjects)
└── code/
    ├── REPO_README.md                 (original repo README)
    ├── models/
    │   ├── m1_RL_sliced2_halfStudentT_reparam_confEs_norm.stan      (202 lines — eRL model)
    │   ├── m1_RL_sliced2_halfStudentT_reparam_confEs_np_norm.stan   (194 lines — RL model)
    │   ├── m2_KF_sliced2_halfStudentT_reparam_confEs_norm.stan      (270 lines — eKF model)
    │   ├── m2_KF_sliced2_halfStudentT_reparam_confEs_np_norm.stan   (267 lines — KF model)
    │   └── m4_C_sliced2_halfStudentT_reparam_confEs_norm.stan       (149 lines — Constant baseline)
    ├── fit_models_cs.R                (137 lines — RStan model fitting driver)
    ├── primary_analysis_cs.R          (1,327 lines — LOO comparison, diagnostics, plotting)
    ├── extra_analysis_cs.R            (27 KB — ESS, BFMI, additional diagnostics)
    ├── stan_utility.R                 (122 lines — MCMC diagnostic utilities)
    ├── hpc_spec.txt                   (HPC spec: 4 chains, 12000 iter, 6000 warmup, adapt_delta=0.95)
    └── preprocessing/
        └── preprocess_data_linear.py  (1,038 lines — MATLAB .mat → Python pickle → Stan data)
```

## Study

Participants (N=27, after exclusions from 33 recruited) received sequences of painful thermal stimuli on their forearm. On each trial they either rated perceived pain intensity (0--100 VAS) or predicted the next stimulus intensity, plus a confidence rating (0--1). The experiment had 4 conditions crossing two factors:

| Condition | Volatility | Signal strength |
|-----------|-----------|-----------------|
| LVLS | Low volatility | Low signal |
| LVHS | Low volatility | High signal |
| HVLS | High volatility | Low signal |
| HVHS | High volatility | High signal |

80 trials per condition (40 perception + 40 prediction, interleaved), 320 trials total.

The core question: does the brain use statistical learning (building expectations from past stimuli) to modulate pain perception, independent of external cues?

## Data format

### Stimulus transformation

Each participant has personalised linear transform coefficients (`a_coeffs.csv`, `b_coeffs.csv`) mapping raw stimulus intensity to a perceptual scale, calibrated via a pre-experiment psychophysics session. The preprocessing script (`preprocess_data_linear.py`) converts raw MATLAB `.mat` files into Stan-ready data.

### Stan data structure (per condition)

The Stan models expect these data fields:

| Field | Shape | Description |
|-------|-------|-------------|
| `N` | scalar | Number of subjects (27) |
| `Tn` | scalar | Number of trials per condition (80) |
| `N_obs_VAS` | scalar | Number of perception ratings (40) |
| `N_obs_CP` | scalar | Number of predictions (varies, ~36--40) |
| `noxInputTrans` | N × Tn | Linearly transformed stimulus intensities |
| `painRating` | N × N_obs_VAS | VAS pain ratings (0--100) |
| `painPred` | N × N_obs_CP | Pain predictions (0--100) |
| `RatingConf` | N × N_obs_VAS | Confidence ratings for perceptions (0--1) |
| `PredConf` | N × N_obs_CP | Confidence ratings for predictions (0--1) |
| `TrialType` | N × Tn | 1=perception, 2=prediction, -1=missing |
| `TrialTypeAllWMissing` | N × Tn | Same but missing coded as 2 (for indexing) |
| `PainValsAll` | N × Tn | All responses combined (-1 for missing) |
| `PercIndexArray` | N × N_obs_VAS | Trial indices of perception ratings |
| `PredIndexArray` | N × N_obs_CP | Trial indices of predictions |
| `IndexMissingAll` | N × N_missing | Trial indices of missing last-in-sequence predictions |

**Note:** Stan data is serialised as Python pickles (.pickle) in the full repo and read as R .rds in the fitting code. The raw .mat participant files and pickle intermediates are available on Zenodo (10.1 GB archive) but not included here due to size.

## Five competing models

All models share the same observation likelihood but differ in their learning dynamics:

### Observation model (shared)

```
response_t ~ Normal(μ_t, σ_t)
σ_t = ξ · exp((1 - conf_t) / cs)
```

where `μ_t` is either perceived pain P_t (perception trials) or expected pain E_{t+1} (prediction trials), selected by trial type. The noise σ_t is **heteroscedastic**: low confidence inflates observation noise exponentially.

### Model 1: eRL (expectation-weighted Reinforcement Learning) — 5 parameters

```
P_t = (1 - γ) · N_t + γ · E_t       ← perception is weighted mix of input and expectation
PE_t = P_t - E_t                      ← prediction error
E_{t+1} = E_t + α · PE_t              ← update rule
```

Parameters: α (learning rate, [0,1]), γ (perceptual weight, [0,1]), ξ (response noise), E₀ (initial expectation), cs (confidence scale).

### Model 2: RL (non-perceptive) — 4 parameters

Same as eRL but γ = 0: perception equals raw input. Tests whether expectation modulates perception.

### Model 3: eKF (expectation-weighted Kalman Filter) — 7 parameters

```
γ_t = ε² / (ε² + ψ² + w_t)           ← perceptual weight (adaptive)
P_t = (1 - γ_t) · N_t + γ_t · E_t
α_t = w_t / (ψ² + w_t)               ← learning rate (adaptive)
E_{t+1} = E_t + α_t · (P_t - E_t)
w_{t+1} = w_t · (ε² + ψ²)/(ε² + ψ² + w_t) + η²  ← uncertainty update
```

Parameters: ε (input noise), ψ (observation noise), η (volatility), ξ (response noise), E₀, w₀ (initial uncertainty), cs.

### Model 4: KF (non-perceptive) — 6 parameters

Same Kalman filter but without the perceptual weighting layer (ε removed, γ = 0).

### Model 5: Constant — 3 parameters

```
μ_t = C    (no learning)
```

Parameters: C (constant response), ξ, cs.

## Hierarchical structure

### Two-level hierarchy

**Group level:** `μ ~ Normal(0, 1)` for each parameter; scale priors use the half-Student-t(3, 0, 1) reparameterisation: `σ = σ_a / √σ_b` where `σ_a ~ Normal(0, 1)` and `σ_b ~ Gamma(1.5, 1.5)`.

**Subject level:** Non-centred parameterisation ("Matt trick"):
- Unbounded parameters (ξ, E₀, w₀, cs): `raw ~ Normal(0, 1); param = μ + σ · raw`
- Bounded parameters (α, γ ∈ [0,1]): `raw ~ Normal(0, 1); param = Φ_approx(μ + σ · raw)` where `Φ_approx` is Stan's fast probit approximation

### Missing data as latent parameters

The last trial in each sequence has a prediction with no subsequent perception to validate it. These missing responses and their confidence ratings are declared as Stan `parameters` (`MissingData`, `MissingConf`) and estimated jointly with model parameters — a sliced missing data pattern.

## Implementation details

- **PPL**: RStan (requires Stan 2.29+)
- **Sampler**: NUTS, 4 chains, 12,000 iterations (6,000 warmup), `adapt_delta = 0.95`, `max_treedepth = 11`
- **Model comparison**: PSIS-LOO via `loo` package (per-subject log-likelihood)
- **Fitting strategy**: 5 models × 4 conditions = 20 separate Stan fits
- **Computational cost**: HPC required; pre-fitted .rds objects total ~10 GB on Zenodo
- **Data pipeline**: MATLAB .mat → Python preprocessing → pickle → R .rds → RStan
- **Known data issue**: Label swap in raw JSON files (HV ↔ LV mislabeled), corrected during preprocessing

## Benchmark value

### A clean model comparison benchmark

Unlike the structural-trap papers, this benchmark tests an agent's ability to implement and compare a family of related cognitive models with progressively richer dynamics. The key question is not "will the wrong model mislead?" but "can the agent correctly implement the full hierarchy of models and draw the right scientific conclusion?"

### Agent evaluation dimensions

1. **Sequential dynamics in Stan** — each model involves a trial-by-trial recurrence relation computed inside a for-loop in `transformed parameters` / `model`. The KF models couple 4 evolving state variables (E, w, α, γ). These loops are not vectorisable.

2. **Heteroscedastic likelihood** — the confidence-dependent noise `ξ · exp((1 - conf) / cs)` is non-standard. An agent defaulting to constant noise would miss the key finding that confidence modulates the observation model.

3. **Non-centred parameterisation with Φ_approx** — bounded parameters use the probit transform `Phi_approx()` in the Matt trick. Using `inv_logit()` instead changes the effective prior. Getting the reparameterisation wrong leads to divergences.

4. **Missing data as parameters** — the sliced missing data pattern requires declaring observations as parameters with appropriate bounds, then inserting them into the data matrix in `transformed parameters`.

5. **Half-Student-t(3,0,1) scale prior** — the `σ_a / √σ_b` decomposition is a specific reparameterisation that avoids sampling a constrained half-Student-t directly. An agent must recognise this pattern or independently derive an equivalent.

6. **Cross-language pipeline** — raw data is MATLAB, preprocessing is Python, fitting is R. An agent must navigate or replicate this pipeline.

### Comparison to other benchmark papers

| Dimension | This paper | Paper 02 (DDM) | Paper 06 (SARS-CoV-2) |
|-----------|-----------|-----------------|----------------------|
| PPL | RStan | Stan (CmdStanR) | Stan (CmdStanR) |
| Model count | 5 competing models | 4 variants | 8 variants |
| Hierarchy | Subject-level (27 subjects) | Subject-level | Household-level (307 households) |
| Dynamics | Trial-by-trial RL/KF | Drift-diffusion (continuous) | Epidemic final-size (discrete) |
| Key challenge | Heteroscedastic noise + adaptive learning | Censoring/truncation | Recursive combinatorial likelihood |
| Missing data | Latent parameters | NA (censored = different) | No |
| Comparison method | PSIS-LOO | LOO-IC | LOO-IC + WBIC |
