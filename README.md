# Bayesian Bench

A benchmark for evaluating Bayesian modelling agents against human experts. Each paper contains raw materials (PDF, data, code) and a comprehensive README documenting the model structure, data format, structural data traps, and evaluation dimensions.

## Papers

| # | Directory | Paper | Year | PPL |
|---|-----------|-------|------|-----|
| 01 | `01_rat_multisensory_psychophysics/` | Zanzi et al. — Rat visual perceptual space | 2025 | Stan |
| 02 | `02_diffusion_truncation_censoring/` | Henrich & Klauer — Diffusion model with censoring | 2026 | Stan |
| 03 | `03_h5n1_viral_kinetics/` | Eales et al. — H5N1 viral kinetics in dairy cattle | 2026 | Stan |
| 04 | `04_systems_biology_multimodel/` | Linden-Santangeli et al. — Systems biology multimodel inference | 2025 | PyMC |
| 05 | `05_rabies_vaccination/` | Ferguson et al. — Rabies vaccination coverage heterogeneity | 2025 | brms/Stan |
| 06 | `06_sars2_household_transmission/` | van Boven et al. — SARS-CoV-2 household transmission | 2024 | Stan |
| 07 | `07_pain_learning/` | Onysk et al. — Pain learning and statistical prediction | 2024 | RStan |
| 08 | `08_supernova_hierarchical/` | Grayling et al. — Hierarchical supernova dust inference | 2024 | NumPyro/JAX |
| 09 | `09_epidemic_turing/` | Koher et al. — Epidemic modelling with survey covariates | 2023 | Turing.jl |
| 10 | `10_brain_cell_counts/` | Dimmock et al. — Hierarchical brain cell counts | 2025 | Stan |
| 11 | `11_nowcasting/` | Lison et al. — Generative Bayesian nowcasting | 2024 | Stan |
| 12 | `12_robust_astronomy/` | Martin & Mortlock — Robust Bayesian regression in astronomy | 2025 | NumPyro/Stan |
| 13 | `13_zero_inflated_beta/` | Pietila — Zero-inflated plant cover, left-censored Beta regression | 2025 | Stan/brms |

## Frontier model knowledge cutoffs (as of Feb 2026)

| Model | Knowledge Cutoff | Training Data Cutoff | Released |
|-------|-----------------|---------------------|----------|
| Claude Opus 4.6 | May 2025 | Aug 2025 | Feb 5, 2026 |
| Claude Sonnet 4.6 | Aug 2025 | Jan 2026 | Feb 17, 2026 |
| GPT-5.2 (all variants) | Aug 31, 2025 | — | Dec 11, 2025 |
| GPT-5.3-Codex | ~Aug 2025 (inherits 5.2) | — | Feb 5, 2026 |
| Gemini 3.1 Pro | Jan 2025 | — | Feb 19, 2026 |

Papers from 2026 (#02 Henrich, #03 Eales) and late 2025 fall near or after these cutoffs, reducing the chance that models have memorised the solutions.
