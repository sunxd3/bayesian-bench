# Bayesian Paper Collection

A curated collection of Bayesian modeling papers with open data and reproducible reference models. Each entry packages the full paper PDF, original datasets, and original code across diverse probabilistic programming ecosystems (Stan, PyMC, NumPyro, Turing.jl, brms).

## Papers

| Directory | Paper | Year | PPL |
|-----------|-------|------|-----|
| `rat_multisensory_psychophysics/` | Zanzi et al. — Rat visual perceptual space | 2025 | Stan |
| `diffusion_truncation_censoring/` | Henrich & Klauer — Diffusion model with censoring | 2026 | Stan |
| `h5n1_viral_kinetics/` | Eales et al. — H5N1 viral kinetics in dairy cattle | 2026 | Stan |
| `systems_biology_multimodel/` | Linden-Santangeli et al. — Systems biology multimodel inference | 2025 | PyMC |
| `rabies_vaccination/` | Ferguson et al. — Rabies vaccination coverage heterogeneity | 2025 | brms/Stan |
| `sars2_household_transmission/` | van Boven et al. — SARS-CoV-2 household transmission | 2024 | Stan |
| `pain_learning/` | Onysk et al. — Pain learning and statistical prediction | 2024 | RStan |
| `supernova_hierarchical/` | Grayling et al. — Hierarchical supernova dust inference | 2024 | NumPyro/JAX |
| `epidemic_turing/` | Koher et al. — Epidemic modelling with survey covariates | 2023 | Turing.jl |
| `brain_cell_counts/` | Dimmock et al. — Hierarchical brain cell counts | 2025 | Stan |
| `nowcasting/` | Lison et al. — Generative Bayesian nowcasting | 2024 | Stan |
| `robust_astronomy/` | Martin & Mortlock — Robust Bayesian regression in astronomy | 2025 | NumPyro/Stan |
| `zero_inflated_beta/` | Pietila — Zero-inflated plant cover, left-censored Beta regression | 2025 | Stan/brms |

## Publication dates vs. LLM knowledge cutoffs (as of Feb 2026)

| Model | Knowledge Cutoff | Training Data Cutoff | Released |
|-------|-----------------|---------------------|----------|
| Claude Opus 4.6 | May 2025 | Aug 2025 | Feb 5, 2026 |
| Claude Sonnet 4.6 | Aug 2025 | Jan 2026 | Feb 17, 2026 |
| GPT-5.2 (all variants) | Aug 31, 2025 | — | Dec 11, 2025 |
| GPT-5.3-Codex | ~Aug 2025 (inherits 5.2) | — | Feb 5, 2026 |
| Gemini 3.1 Pro | Jan 2025 | — | Feb 19, 2026 |

### Data contamination risk by paper

Assessed against the latest knowledge cutoff (Aug 31, 2025 — GPT-5.2 / Sonnet 4.6).

**Post-cutoff** (published after Aug 31, 2025 — lowest contamination risk):

| Paper | Published |
|-------|-----------|
| Zanzi et al. — Rat psychophysics | Oct 29, 2025 |
| Henrich & Klauer — Diffusion censoring | 2026 |
| Eales et al. — H5N1 viral kinetics | Jan 2026 |
| Dimmock et al. — Brain cell counts | Nov 21, 2025 |

**Near cutoff** (published within weeks of Aug 31, 2025):

| Paper | Published |
|-------|-----------|
| Linden-Santangeli et al. — Systems biology | Aug 11, 2025 |
| Martin & Mortlock — Robust astronomy | Aug 13, 2025 |

**Pre-cutoff** (published before Aug 2025 — may appear in training data):

| Paper | Published |
|-------|-----------|
| Ferguson et al. — Rabies vaccination | May 5, 2025 |
| van Boven et al. — SARS-CoV-2 household | 2024 |
| Onysk et al. — Pain learning | 2024 |
| Grayling et al. — Supernova | 2024 |
| Koher et al. — Epidemic Turing | 2023 |
| Lison et al. — Nowcasting | Apr 2024 |
| Pietila — Zero-inflated beta | Jun 2025 |

Pre-cutoff papers remain valuable for testing modelling ability rather than memorisation — but post-cutoff papers have the strongest guarantee against data contamination.
