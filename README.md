# Bayesian Paper Collection

A curated collection of Bayesian modeling papers with open data and reproducible reference models. Each entry packages the full paper PDF, original datasets, and original code across diverse probabilistic programming ecosystems (Stan, PyMC, NumPyro, Turing.jl, brms).

## Papers

Papers are organised by publication date relative to the latest LLM knowledge cutoff (Aug 31, 2025 — GPT-5.2 / Sonnet 4.6).

### `post_cutoff/` — published after Aug 2025

| Directory | Paper | Published | PPL |
|-----------|-------|-----------|-----|
| `rat_multisensory_psychophysics/` | Zanzi et al. — Rat visual perceptual space | Oct 2025 | Stan |
| `diffusion_truncation_censoring/` | Henrich & Klauer — Diffusion model with censoring | 2026 | Stan |
| `h5n1_viral_kinetics/` | Eales et al. — H5N1 viral kinetics in dairy cattle | Jan 2026 | Stan |
| `brain_cell_counts/` | Dimmock et al. — Hierarchical brain cell counts | Nov 2025 | Stan |

### `around_cutoff/` — published within weeks of Aug 2025

| Directory | Paper | Published | PPL |
|-----------|-------|-----------|-----|
| `systems_biology_multimodel/` | Linden-Santangeli et al. — Systems biology multimodel inference | Aug 2025 | PyMC |
| `robust_astronomy/` | Martin & Mortlock — Robust Bayesian regression in astronomy | Aug 2025 | NumPyro/Stan |

### `pre_cutoff/` — published before Aug 2025

| Directory | Paper | Published | PPL |
|-----------|-------|-----------|-----|
| `rabies_vaccination/` | Ferguson et al. — Rabies vaccination coverage heterogeneity | May 2025 | brms/Stan |
| `zero_inflated_beta/` | Pietila — Zero-inflated plant cover, left-censored Beta regression | Jun 2025 | Stan/brms |
| `nowcasting/` | Lison et al. — Generative Bayesian nowcasting | Apr 2024 | Stan |
| `sars2_household_transmission/` | van Boven et al. — SARS-CoV-2 household transmission | 2024 | Stan |
| `pain_learning/` | Onysk et al. — Pain learning and statistical prediction | 2024 | RStan |
| `supernova_hierarchical/` | Grayling et al. — Hierarchical supernova dust inference | 2024 | NumPyro/JAX |
| `epidemic_turing/` | Koher et al. — Epidemic modelling with survey covariates | 2023 | Turing.jl |

## LLM knowledge cutoffs (as of Feb 2026)

| Model | Knowledge Cutoff | Training Data Cutoff | Released |
|-------|-----------------|---------------------|----------|
| Claude Opus 4.6 | May 2025 | Aug 2025 | Feb 5, 2026 |
| Claude Sonnet 4.6 | Aug 2025 | Jan 2026 | Feb 17, 2026 |
| GPT-5.2 (all variants) | Aug 31, 2025 | — | Dec 11, 2025 |
| GPT-5.3-Codex | ~Aug 2025 (inherits 5.2) | — | Feb 5, 2026 |
| Gemini 3.1 Pro | Jan 2025 | — | Feb 19, 2026 |

Pre-cutoff papers remain valuable for testing modelling ability rather than memorisation — but post-cutoff papers have the strongest guarantee against data contamination.
