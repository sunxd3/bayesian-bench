# Paper 4: Systems Biology Multimodel Inference (PyMC-SMC)

## Citation

Linden-Santangeli N, Zhang J, Kramer B, Rangamani P (2025). Increasing certainty in
systems biology models using Bayesian multimodel inference. *Nature Communications* 16:
7307. https://doi.org/10.1038/s41467-025-62415-4

## Links

- **Paper:** https://www.nature.com/articles/s41467-025-62415-4
- **PMC:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12339951/
- **GitHub:** https://github.com/RangamaniLabUCSD/multimodel-inference
- **Zenodo:** https://doi.org/10.5281/zenodo.15129141 (v1.1.0, 427 MB)
- **arXiv:** https://arxiv.org/abs/2406.11178
- **Experimental data origin:** Keyes et al. (eLife 2020) — subcellular ERK activity

## Folder Structure

```
04_systems_biology_multimodel/
├── README.md
├── paper/
│   └── s41467-025-62415-4.pdf                              # Full paper (16 pages, 3.8 MB)
├── data/
│   ├── HF_96_synthetic_data.csv                            # Synthetic dose-response (10 points)
│   ├── HF_96_traj_data.json                                # Synthetic trajectories (3 conditions)
│   ├── Keyes_et_al_2020-fig1-data1-v2-CYTO.json            # Experimental: cytoplasm ERK (76 cells)
│   ├── Keyes_et_al_2020-fig1-data1-v2-PM.json              # Experimental: plasma membrane ERK
│   └── model_info_supplemental_material.xlsx               # Model details supplemental table
└── code/
    ├── REPO_README.md
    ├── utils.py                                            # Core: build_pymc_model, smc_pymc, priors
    ├── diffrax_ODE_PyTensor.py                             # Custom PyTensor Op for JAX ODE solver
    ├── plotting_helper_funcs.py                            # Plotting utilities
    ├── process_data.ipynb                                  # Data generation and preprocessing
    ├── model_validation.ipynb                              # Model validation
    ├── models/                                             # 13 ODE model files (equinox.Module)
    │   ├── huang_ferrell_1996.py                           # 26 free params, mass action
    │   ├── kholodenko_2000.py                              # 6 free params, Michaelis-Menten
    │   ├── levchenko_2000.py                               # 25 free params
    │   ├── hornberg_2005.py                                # 22 free params
    │   ├── birtwistle_2007.py                              # 20 free params
    │   ├── orton_2009.py                                   # 5 free params, PI3K/Akt/Rap1
    │   ├── vonKriegsheim_2009.py                           # 18 free params, neg feedback
    │   ├── shin_2014.py                                    # 7 free params
    │   ├── ryu_2015.py                                     # 5 free params, Hill equations
    │   ├── kochanczyk_2017.py                              # 11 free params
    │   ├── ryu_2015_Rap1.py                                # Rap1 variant
    │   ├── shin_2014_Rap1.py                               # Rap1 variant
    │   └── vonKriegsheim_2009_Rap1.py                      # Rap1 variant
    ├── param_est/                                          # Parameter estimation via SMC
    │   ├── model_info.json                                 # Per-model config: params, priors, solver
    │   ├── inference_process_traj.py                       # CLI: trajectory inference
    │   ├── inference_process_dose_response.py              # CLI: dose-response inference
    │   ├── inference_process_location_diff_traj.py         # CLI: location-specific inference
    │   ├── HF96_DR_analyze.py                              # Analysis: synthetic dose-response
    │   ├── HF96_traj_analyze.py                            # Analysis: synthetic trajectories
    │   ├── Keyes_analyze.py                                # Analysis: experimental data
    │   ├── Keyes_data_len_analyze.py                       # Analysis: data length sensitivity
    │   ├── Keyes_data_quality_inference.py                 # Analysis: data quality sensitivity
    │   └── Keyes_location_diff_analyze.py                  # Analysis: location-specific
    └── multimodel_inference/                               # MMI computation
        ├── multimodel_inf_Keyes_data.py                    # Core MMI: BMA, pseudo-BMA, stacking
        ├── Keyes_MMI.ipynb                                 # Main experimental MMI analysis
        ├── HF96_DR_MMI.ipynb                               # Synthetic dose-response MMI
        └── HF96_traj_MMI.ipynb                             # Synthetic trajectory MMI
```

**Not included locally:** Pre-computed SMC results (~400 MB of JSON/NPY files in
`results/MAPK/param_est/`). Available from Zenodo. Julia identifiability analysis
files and GSA scripts also omitted.

## Study Overview

The ERK/MAPK signaling pathway (EGF → EGFR → Ras → Raf → MEK → ERK) is a core
intracellular cascade regulating cell proliferation and survival. The BioModels database
contains 125+ ODE models of this pathway, differing in complexity, kinetic formulation,
and feedback structure. No single model is definitively "correct."

This paper uses **Bayesian multimodel inference (MMI)** — fitting 10 competing ODE models
independently via SMC, then combining their predictions weighted by data support — to
produce ensemble predictions that are more robust and certain than any individual model.

## The 10 Competing Models

| Abbreviation | Reference | Free params | Key features |
|---|---|---|---|
| H' 1996 | Huang & Ferrell | 26 | Mass action, no feedback, 16 state vars |
| K' 2000 | Kholodenko | 6 | Michaelis-Menten, 6 state vars |
| L' 2000 | Levchenko | 25 | Long integration times (100,000s) |
| H' 2005 | Hornberg | 22 | Very small EGF conversion (1e-9) |
| B' 2007 | Birtwistle | 20 | Detailed receptor dynamics |
| O' 2009 | Orton | 5 | 26 state vars, PI3K/Akt/Rap1 crosstalk |
| vK' 2009 | von Kriegsheim | 18 | Negative feedback loop |
| S' 2014 | Shin/Sturm | 7 | Newton solver for steady-state |
| R' 2015 | Ryu | 5 | Neg + pos feedback, Hill equations |
| K' 2017 | Kochanczyk | 11 | Newton solver, feedback |

Each model is implemented as a JAX-compatible `equinox.Module` in `code/models/`.
Parameters range from 5 to 26, creating dramatic complexity asymmetry.

## Data

### Synthetic data (ground truth known)

- **Dose-response:** 10 EGF concentrations (0.001–0.106 nM), steady-state ERK
  normalized to % max. Generated from Huang-Ferrell 1996 model + Gaussian noise (σ=0.1).
- **Trajectories:** 3 time-dependent ERK activity curves at different EGF concentrations,
  30-minute duration.

### Experimental data (EKAR4 biosensor)

- **Source:** Keyes et al. (eLife 2020) — live-cell imaging of subcellular ERK
- **Measurement:** YFP/CFP emission ratios from EKAR4 kinase activity reporter
- **Locations:** Cytoplasm (CYTO) and plasma membrane (PM) — measured separately
- **Conditions:** With and without Rap1 inhibition
- **Sample size:** 76 individual cell recordings, 40-minute time series
- **Format:** JSON files with time points and normalized ERK activity

## Bayesian Inference: PyMC SMC

### Why SMC instead of NUTS

The key reason is that SMC provides the **log marginal likelihood** as a free byproduct
of the tempering sequence — the quantity needed for Bayesian model averaging. NUTS does
not provide marginal likelihoods.

Additionally, the ODE solver wrapped in a custom PyTensor Op does not provide gradients
(`grad()` raises `NotImplementedError`), making gradient-based samplers like NUTS
unusable. SMC uses Independent Metropolis-Hastings (IMH), which only requires likelihood
evaluations.

### How it works

1. **Stage 0 (β=0):** Sample particles from the prior
2. **Intermediate stages:** Gradually increase inverse temperature β from 0 to 1,
   targeting `p(θ)^{1-β} · p(data|θ)^β`
3. **Stage N (β=1):** Particles approximate the full posterior
4. **Marginal likelihood:** Accumulated normalizing constants across stages →
   `idata.sample_stats["log_marginal_likelihood"]`

### Configuration

```python
pm.smc.sample_smc(draws=500, chains=4, threshold=0.85,
                  correlation_threshold=0.01)
```

### Prior specification

Log-normal priors on all parameters, centered on literature values:
```
prior mean = log(θ_nominal)
prior σ = 2.350  (so 95% mass in [10^{-2} × θ_nominal, 10^2 × θ_nominal])
```

Computed using PreLiz's `maxent()` function. Model priors: uniform `p(M_k) = 1/K`.

### ODE integration stack

```
equinox.Module (model) → diffrax (JAX ODE solver) → custom PyTensor Op
    → PyMC model → SMC sampler
```

The `StimRespOp` class in `diffrax_ODE_PyTensor.py` + `jax_funcify` registration
in `utils.py` are the glue code. The forward pass calls `jax.jit(simulator)`;
the gradient is never computed.

## Model Comparison: Three Weighting Schemes

### 1. Bayesian Model Averaging (BMA) — from SMC marginal likelihoods

```python
w_BMA[k] ∝ p(M_k) · p(data | M_k)
```

Uses `log_sum_exp` for numerical stability. Favors models with high marginal likelihood
(evidence). Sensitive to prior specification (Lindley's paradox).

### 2. Pseudo-BMA — from PSIS-LOO-CV

```python
az.compare(idata_dict, ic='loo', method='BB-pseudo-BMA')
```

Uses ArviZ PSIS-LOO to estimate ELPD per model, then normalizes via softmax.
Favors out-of-sample predictive performance.

### 3. Stacking of Predictive Densities

```python
az.compare(idata_dict, ic='loo', method='stacking')
```

Jointly optimizes weights to maximize the combined LOO predictive density. Can
concentrate all weight on a single model (e.g., Orton 2009) while BMA distributes
weight more evenly.

### Combined predictive distribution

```
p_combined(y*) = Σ_k w_k · p(y* | data, M_k)
```

## Implementation Details

- **PPL:** PyMC 5.10.4 with PyTensor 2.18.6 backend
- **ODE solver:** diffrax 0.5.0 (JAX-based), equinox for model definitions
- **Sampler:** SMC with IMH kernel, 4 chains × 500 draws, threshold=0.85
- **Diagnostics:** ArviZ 0.16.1 for PSIS-LOO-CV, pseudo-BMA, stacking
- **Prior elicitation:** PreLiz 0.3.6 `maxent()`
- **Timeout:** `func_timeout` to handle ODE solver failures
- **Language:** Python 3.11.6 (Julia 1.10.0 for identifiability analysis only)

## Benchmark Value

### What this tests

This paper tests capabilities **completely absent from papers #1–3**:

1. **Multi-model orchestration.** An agent must run 10 independent Bayesian models
   with different structures, parameter counts (5–26), solver configurations, and
   unit conventions, then combine their results.

2. **SMC sampling (not NUTS).** The agent must know that `pm.smc.sample_smc()` exists,
   how to configure it, and that it provides marginal likelihoods. Most LLM training
   data emphasizes NUTS.

3. **Marginal likelihood extraction.** The key quantity (`log_marginal_likelihood`)
   is buried in `idata.sample_stats` — the agent must know where to find it.

4. **Three-way model comparison.** BMA (from marginal likelihoods) vs. pseudo-BMA
   (from PSIS-LOO) vs. stacking — conceptually different approaches with different
   strengths. The agent must implement and compare all three.

5. **ODE solver integration.** JAX/diffrax/equinox → custom PyTensor Op → PyMC.
   The agent must either reproduce this glue code or find an alternative (e.g., `sunode`,
   `scipy.integrate` with a PyMC potential).

6. **Prior sensitivity awareness.** Marginal likelihoods are notoriously prior-sensitive
   (Lindley's paradox). Diffuse/uniform priors will produce drastically different BMA
   weights than the paper's calibrated log-normal priors.

### Traps

1. **Prior specification dominates BMA weights.** An agent using uniform or overly wide
   priors will get completely different model rankings — and may not realize why.

2. **SMC marginal likelihood variance.** Each chain gives a different estimate. The agent
   must average across chains, not take a single value.

3. **Model complexity asymmetry.** Under BMA, the Occam factor in the marginal likelihood
   integral strongly favors simpler models (5 params vs 26 params), independent of actual
   fit quality. Stacking compensates for this, concentrating weight differently.

4. **ODE solver failures.** Some models with certain parameter draws will fail to integrate
   (stiff systems, Newton non-convergence). The agent must handle NaN/timeout gracefully.

5. **Unit conversion heterogeneity.** EGF conversion factors range from 1e-9 to 6048
   across models. Time conversion factors differ. Getting these wrong produces nonsensical
   likelihoods that won't cause an error — just wrong results.

### Comparison to papers #1–3

| Dimension | Papers #1–3 | Paper #4 |
|---|---|---|
| Sampler | NUTS/HMC | SMC (no gradients) |
| Models compared | 1–7 variants of one model | 10 structurally different ODE models |
| Model comparison | PSIS-LOO | BMA + pseudo-BMA + stacking |
| Key quantity | ELPD | Marginal likelihood |
| Likelihood | Analytical | Requires numerical ODE solve |
| Failure modes | Divergences, R-hat | ODE solver crash, timeout, NaN |
| Prior sensitivity | Moderate | Critical (Lindley's paradox) |
| Dependencies | PyMC or Stan | PyMC + JAX + diffrax + equinox + PreLiz |
