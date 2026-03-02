# Bayesian Paper Collection

Bayesian modeling papers with open data and reference code. Each entry has the paper PDF, datasets, and original code. Organised by publication date relative to LLM knowledge cutoffs (~Aug 2025).

## Selection criteria

Each paper was included because it meets all of the following:

1. **Open data and code** — paper PDF, datasets, and reference implementation are all publicly available.
2. **Non-trivial Bayesian workflow** — the modeling goes beyond textbook examples: custom likelihoods, reparameterizations, numerical stability tricks, or multi-model comparison.
3. **Human insight required** — the code contains modeling decisions that are hard to derive from first principles (e.g., data augmentation for censoring, recursive combinatorial likelihoods, power-mean aggregation with estimated exponent).
4. **Complete reference material** — the paper, data, and code together provide enough information to replicate the modeling approach, even if the original scripts aren't runnable out of the box.
5. **Diverse PPL coverage** — the collection spans Stan, PyMC, NumPyro, Turing.jl, and brms to cover the major probabilistic programming ecosystems.

Papers were dropped if they lacked reference code entirely or offered only standard textbook-level modeling.

### `post_cutoff/`

| Directory | Paper | PPL |
|-----------|-------|-----|
| `rat_multisensory_psychophysics/` | Zanzi et al. (Oct 2025) — Rat visual perceptual space | Stan |
| `diffusion_truncation_censoring/` | Henrich & Klauer (2026) — Diffusion model with censoring | Stan |
| `brain_cell_counts/` | Dimmock et al. (Nov 2025) — Hierarchical brain cell counts | Stan |

### `around_cutoff/`

| Directory | Paper | PPL |
|-----------|-------|-----|
| `systems_biology_multimodel/` | Linden-Santangeli et al. (Aug 2025) — Systems biology multimodel inference | PyMC |
| `robust_astronomy/` | Martin & Mortlock (Aug 2025) — Robust Bayesian regression in astronomy | NumPyro/Stan |

### `pre_cutoff/`

| Directory | Paper | PPL |
|-----------|-------|-----|
| `rabies_vaccination/` | Ferguson et al. (May 2025) — Rabies vaccination coverage heterogeneity | brms/Stan |
| `nowcasting/` | Lison et al. (Apr 2024) — Generative Bayesian nowcasting | Stan |
| `sars2_household_transmission/` | van Boven et al. (2024) — SARS-CoV-2 household transmission | Stan |
| `pain_learning/` | Onysk et al. (2024) — Pain learning and statistical prediction | RStan |
| `supernova_hierarchical/` | Grayling et al. (2024) — Hierarchical supernova dust inference | NumPyro/JAX |
| `epidemic_turing/` | Koher et al. (2023) — Epidemic modelling with survey covariates | Turing.jl |
