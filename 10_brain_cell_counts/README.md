# 10 — Hierarchical Bayesian Brain Cell Counts (Stan)

## Citation

Dimmock S, Exley BMS, Moore G, Menage L, Delogu A, Schultz SR, Warburton EC, Houghton CJ, O'Donnell C (2025).
**Hierarchical Bayesian modeling of multiregion brain cell count data.**
*eLife* 14:RP102391.

## Links

| Resource | URL |
|----------|-----|
| Paper | https://elifesciences.org/articles/102391 |
| bioRxiv | https://doi.org/10.1101/2024.07.20.603979 |
| GitHub (paper) | https://github.com/BayesianCellCounts/DimmockEtAl2025 |
| GitHub (tutorial) | https://github.com/BayesianCellCounts/R_example |
| Project site | https://bayesiancellcounts.github.io/ |
| Zenodo (code) | https://doi.org/10.5281/zenodo.16340994 |
| Zenodo (CS1 data) | https://doi.org/10.5281/zenodo.12787211 |
| Zenodo (CS2 data) | https://doi.org/10.5281/zenodo.12787287 |

## Folder structure

```
10_brain_cell_counts/
├── README.md                                      (this file)
├── paper/
│   └── elife-102391.pdf                           (11 MB — full eLife paper)
├── data/
│   ├── case_study_1/
│   │   └── data.csv                               (4,533 rows × 14 cols — cFos rat data)
│   └── case_study_2/
│       └── data.csv                               (501 rows × 10 cols — Sox14 mouse data)
└── code/
    ├── case_study_1/
    │   ├── models/
    │   │   └── model_poiss.stan                   (39 lines — baseline Poisson)
    │   ├── sample_posterior.r
    │   ├── plot_diagnostics.r
    │   ├── plot_ppc.r
    │   └── plot_results.r
    ├── case_study_2/
    │   ├── models/
    │   │   ├── model_hs.stan                      (50 lines — horseshoe shrinkage)
    │   │   ├── model_inflate.stan                 (97 lines — zero-inflated Poisson, optimised)
    │   │   └── model_inflate_base.stan            (similar — unoptimised reference)
    │   ├── sample_posterior.r
    │   ├── plot_diagnostics_hs.r
    │   ├── plot_diagnostics_inflate.r
    │   ├── plot_ppc.r
    │   └── plot_results.r
    └── R_example/
        ├── models/
        │   └── model_hs.stan                      (tutorial horseshoe, sparse-obs format)
        ├── sample_posterior.r
        ├── plot_diagnostics.r
        ├── plot_ppc.r
        ├── plot_results.r
        └── data.csv                               (tutorial data)
```

## Study

Post-mortem cFos immunostaining or genetic labelling lets neuroscientists count active cells across entire brains. The data is inherently **under-sampled**: each animal provides one snapshot, but there are 20–50 brain regions to compare across only ~10 animals per group. Classical region-by-region t-tests are both underpowered and fail to share statistical strength across regions.

### Case Study 1: cFos in rats (2×2 design)

Nucleus reuniens lesion vs. sham surgery × novel vs. familiar stimulus. ~10 rats per group, up to 23 brain regions. Cell counts from cFos immunostaining (a marker of recent neuronal activity). The scientific question: which brain regions respond differently to novelty, and does the lesion disrupt this?

### Case Study 2: Sox14-expressing neurons in mice

Sox14/GFP heterozygotes vs. Sox14 null mice, ~50 brain regions. Counts are sparse — many exact zeros where a cell type is simply absent in a region. This drives the need for zero-inflated models.

## Data format

### `case_study_1/data.csv` (4,533 rows × 14 columns)

| Column | Description |
|--------|-------------|
| `Brain` | Brain section ID |
| `BrainRegion` | Abbreviation (ACC, PRL, DPC, IFC, M2C, VDG, ...) |
| `BregmanRegion` | Anterior-posterior coordinate (mm from Bregma) |
| `Cfos_Counts` | Raw integer cFos cell counts |
| `counts_per_area` | Normalised counts per area |
| `rat_ID` | Animal identifier |
| `group` | Factorial: `Lesion_Novel`, `Sham_Novel`, etc. |
| `surgery` | `Lesion` or `Sham` |
| `stimulus` | `Novel` or `Familiar` |
| `VOID_PERCENTAGE` | Fraction of section that was void/damaged |

### `case_study_2/data.csv` (501 rows × 10 columns)

| Column | Description |
|--------|-------------|
| `counts` | Integer cell counts (many zeros) |
| `regions` | ~50 brain region abbreviations |
| `group` | `het` (heterozygote) or `null` |
| `hemisphere` | `L` or `R` |
| `ID` | Animal ID |
| `group_idx`, `animal_idx`, `region_idx` | Stan-ready integer indices |

## Three Stan models

### Shared structure

All models use a hierarchical log-linear Poisson framework:

```
gamma[a] = theta[g] + tau[g] .* gamma_raw[a] (.* kappa[a])
y[a,r] ~ Poisson(exp(gamma[a,r]))
```

where `theta[g][r]` is the group×region population mean on the log-count scale, `tau[g][r]` is the between-animal scale, and `gamma_raw[a]` captures standardised animal-level deviations (non-centred parameterisation).

### `model_poiss.stan` — baseline hierarchical Poisson

The simplest model. Uses a sparse observation format with an exposure offset `E`:

```stan
y ~ poisson_log(E + gamma)
```

No zero-inflation, no shrinkage. Serves as a comparison baseline.

### `model_hs.stan` — horseshoe shrinkage

Adds a **local shrinkage multiplier** `kappa[a][r]` per animal per region:

```stan
gamma[a] = theta[g] + tau[g] .* gamma_raw[a] .* kappa[a]
kappa[a] ~ HalfNormal(0, 1)   // via <lower=0> constraint
```

When `kappa` is near zero, the animal's deviation is shrunk toward the group mean. When `kappa` is large, the animal retains its individual effect. This is a multiplicative horseshoe variant — `tau[g]` plays the global scale role and `kappa[a]` provides local adaptation.

### `model_inflate.stan` — zero-inflated Poisson

Handles the excess zeros in Case Study 2 via a mixture model:

```stan
pi ~ Beta(1, 5)

if y[n] == 0:
    target += log_sum_exp(log(pi), log1m(pi) + poisson_log_lpmf(0 | gamma[n]))
else:
    target += log1m(pi) + poisson_log_lpmf(y[n] | gamma[n])
```

A zero can arise from the "structural zero" component (probability `pi` — cell type absent from this region) or from the Poisson component (probability `(1-pi) * Poisson(0|lambda)` — cell type present but count happens to be zero). The optimised version pre-indexes zero and non-zero observations in `transformed data` for efficiency.

### Parameters

| Parameter | Prior | Description |
|-----------|-------|-------------|
| `theta[g][r]` | Normal(5, 2) | Group×region population mean (log-rate) |
| `tau[g][r]` | HalfNormal(0, log(1.05)) | Between-animal scale |
| `gamma_raw[a]` | Normal(0, 1) | Standardised animal random effects |
| `kappa[a][r]` | HalfNormal(0, 1) | Local shrinkage (horseshoe model only) |
| `pi` | Beta(1, 5) | Structural zero probability (ZI model only) |

## Implementation details

- **PPL**: Stan via CmdStanR
- **Language**: R (71.6%), Stan (16.6%)
- **Sampler**: NUTS (default CmdStan settings)
- **Model comparison**: Posterior predictive checks (y_rep vs. observed)
- **License**: BSD-3-Clause (code), CC-BY 4.0 (data)

## Benchmark value

### Teaching paper with progressive model building

This is designed as a tutorial — the paper walks through three models of increasing sophistication, each motivated by a specific data pathology. The benchmark tests whether an agent can recognise each pathology and select the appropriate model.

### Agent evaluation dimensions

1. **Recognising under-sampling** — A × R matrix with A < R (fewer animals than regions). Standard frequentist tests are underpowered; the agent must recognise the need for partial pooling.

2. **Sparse-to-dense data format** — Case Study 1 uses sparse format (one row per observation), Case Study 2 uses dense format (A × R matrix). The agent must adapt the Stan model's data block accordingly.

3. **Zero-inflation diagnosis** — Case Study 2 has massive exact zeros. The agent must recognise that a standard Poisson model cannot distinguish "cell type absent" from "count happens to be zero" and implement the ZI mixture.

4. **Horseshoe shrinkage** — The multiplicative `kappa` mechanism is unusual. An agent must understand that `kappa ≈ 0` shrinks an animal toward the group mean while `kappa ≫ 0` lets it deviate freely. This is not the standard Piironen-Vehtari horseshoe with a global-local decomposition.

5. **Posterior predictive checks** — All models generate `y_rep` for PPC. The agent should use these to diagnose model misfit (e.g., the baseline Poisson under-predicting zeros).

### Comparison to other benchmark papers

| Dimension | This paper | Paper 05 (Rabies) | Paper 06 (SARS-CoV-2) |
|-----------|-----------|-------------------|----------------------|
| PPL | Stan (CmdStanR) | brms + Stan | Stan (CmdStanR) |
| Hierarchy | Group × region × animal | District × village | Household × individual |
| Key challenge | Zero-inflation + shrinkage | Power mean custom Stan | Recursive combinatorial |
| Zeros | Structural zeros (cell absent) | Not an issue | Not applicable |
| Model comparison | PPC (y_rep) | brms loo() | LOO-IC + WBIC |
