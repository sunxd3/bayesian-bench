# Candidate Papers for Bayesian Modeling Agent Benchmark

## 1. Rat Multisensory Psychophysics (PyMC)
- **Title:** Seeing what you hear: Compression of rat visual perceptual space by task-irrelevant sounds
- **Authors:** Zanzi, Rinaldi, Fornasaro, Piasini, Zoccolan
- **Journal:** PLOS Computational Biology (Oct 2025)
- **PPL:** PyMC v5.2.0, NUTS (4 chains, 2000 draws, 1000 tuning, target_accept=0.99)
- **Data/Code:** Zenodo DOI 10.5281/zenodo.17280352
- **Paper:** https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013608
- **PMC:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12571299/
- **Model:** Bayesian hierarchical ideal observer. Trial-level choice likelihood per rat, hierarchical priors (partial pooling across rats), NUTS sampling. Prior predictive checks, posterior predictive checks, LOO/ELPD model comparison (PSIS-LOO).
- **Benchmark value:** Very clean reproduce target. Full pipeline on Zenodo. Good for reproduction (run existing code) and re-implementation (rebuild from paper text).

## 2. Diffusion Model with Truncation/Censoring (Stan)
- **Title:** Modeling truncated and censored data with the diffusion model in Stan
- **Authors:** Henrich & Klauer
- **Journal:** Behavior Research Methods (accepted Aug 2025, issue 2026)
- **PPL:** Stan
- **Data/Code:** OSF (project listed in paper) + GitHub
- **PMC:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12819533/
- **Model:** Hierarchical Bayesian cognitive (diffusion) model handling censoring and truncation. Careful likelihood construction in Stan.
- **Benchmark value:** Stan-centric alternative. Good for testing agent's ability to handle censoring/truncation in likelihood.

## 3. H5N1 Viral Kinetics in Dairy Cattle (Stan)
- **Title:** Modeling of H5N1 influenza virus kinetics during dairy cattle infection suggests the timing of infectiousness
- **Journal:** 2026
- **PPL:** Stan, NUTS
- **Data/Code:** GitHub, archived on Zenodo
- **PMC:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12782433/
- **Model:** Bayesian hierarchical model for Ct trajectories.
- **Benchmark value:** Hierarchical time-series / censoring target.

## 4. Systems Biology Multimodel Inference (PyMC-SMC)
- **Title:** Increasing certainty in systems biology models using Bayesian multimodel inference
- **Journal:** Nature Communications (2025)
- **PPL:** PyMC (SMC sampler for posterior samples + marginal likelihood estimates)
- **Data/Code:** Zenodo + source data with paper
- **Paper:** https://www.nature.com/articles/s41467-025-62415-4
- **Model:** ODE-based systems biology models. PyMC SMC for marginal likelihood and model comparison.
- **Benchmark value:** Model comparison / marginal likelihood target. Tests agent on SMC and Bayes factors.

## 5. Synaptic Dynamics Across Sleep-Wake (PyMC)
- **Title:** A unified framework to model synaptic dynamics during the sleep-wake cycle
- **Journal:** PLOS Biology (Jun 2025)
- **PPL:** PyMC + ArviZ
- **Data/Code:** Zenodo (code) + CRCNS (spike-pattern records) + supporting spreadsheets
- **Paper:** https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3003198
- **Model:** Bayesian MCMC for synaptic dynamics parameters.
- **Benchmark value:** Neuroscience domain with publicly available electrophysiology data.

## 6. Rabies Vaccination Campaign Effectiveness (Stan/brms)
- **Title:** Improved effectiveness of vaccination campaigns against rabies by reducing spatial heterogeneity in coverage
- **Journal:** PLOS Biology (May 2025)
- **PPL:** Stan via RStan and brms
- **Data/Code:** Zenodo (deidentified data + code)
- **Paper:** https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002872
- **Model:** Bayesian GLMMs for campaign effectiveness.
- **Benchmark value:** Applied GLMM target. Tests agent on brms/Stan GLMM workflow.

## 7. SARS-CoV-2 Household Transmission (Stan)
- **Title:** Estimation of introduction and transmission rates of SARS-CoV-2 in a prospective household study
- **Journal:** PLOS Computational Biology (Jan 2024)
- **PPL:** Stan via CmdStanR
- **Data/Code:** GitHub, archived on Zenodo
- **Paper:** https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011832
- **Model:** Bayesian transmission model for household epidemiology.
- **Benchmark value:** Epidemiological modeling with clear data+code pipeline.

## 8. Pain Learning and Statistical Prediction (RStan)
- **Title:** Statistical learning shapes pain perception and prediction independently of external cues
- **Journal:** eLife (2024)
- **PPL:** RStan (NUTS/HMC)
- **Data/Code:** Zenodo (all code and data openly available)
- **PMC:** https://pmc.ncbi.nlm.nih.gov/articles/PMC11236420/
- **Model:** Hierarchical Bayesian parameter estimation for pain learning/prediction.
- **Benchmark value:** Quick reproduce target. Clean Zenodo deposit with RStan.

## 9. Hierarchical Supernova Population Modeling (NumPyro)
- **Title:** Scalable hierarchical BayeSN inference: Investigating dependence of SN Ia host galaxy dust properties on stellar mass and redshift
- **Journal:** arXiv 2024 (accepted MNRAS)
- **PPL:** NumPyro (JAX), NUTS
- **Data/Code:** BayeSN GitHub repo; inputs from Pantheon+ public survey
- **Paper:** https://arxiv.org/html/2401.08755v1
- **Model:** Large-scale hierarchical Bayesian model for supernova standardization. Moved from Stan to NumPyro for scalability.
- **Benchmark value:** "Hard mode" hierarchical model. Good for testing agent on large-scale inference and NumPyro/JAX.

## 10. Epidemic Modelling with Survey Covariates (Turing.jl)
- **Title:** Epidemic modelling of monitoring public behavior using surveys during pandemic-induced lockdowns
- **Journal:** Communications Medicine (2023)
- **PPL:** Julia + Turing.jl
- **Data/Code:** GitHub (source code) + Zenodo (hospitalization + predictor data)
- **Paper:** https://www.nature.com/articles/s43856-023-00310-z
- **Model:** Bayesian epidemic model with survey-derived behavioral covariates.
- **Benchmark value:** Julia/Turing.jl target. Tests agent across PPL ecosystems.

---

## Tutorial Papers (Methodological Benchmarks)

### 11. Hierarchical Bayesian Modeling of Brain Cell Counts (Stan) — "Noisy / Undersampled Data" Benchmark
- **Title:** Hierarchical Bayesian modeling of multiregion brain cell count data
- **Authors:** Dimmock et al.
- **Journal:** eLife (2025)
- **PPL:** Python or R / Stan
- **Data/Code:** Openly attached to the eLife publication
- **Model:** Multiregion brain cell counts with undersampled data (high dimensionality, few animals), massive zero-counts, overdispersion. Requires partially pooled hierarchical structure with zero-inflated Poisson likelihood and Horseshoe prior.
- **Benchmark value:** Tests agent's ability to handle messy biological data — recognize that standard models will fail, handle zeros and overdispersion, choose appropriate likelihood (zero-inflated Poisson), and implement regularizing priors (Horseshoe). Designed as a teaching paper so the workflow is meticulously documented.

---

## Structural Data Trap Papers

Papers where naive/standard models will compile and run but produce **scientifically wrong** conclusions. The agent must deduce real-world physical constraints and invent a bespoke generative architecture.

### 12. Generative Bayesian Nowcasting from Line Lists (Stan) — "Right-Truncation & Missingness" Trap
- **Title:** Generative Bayesian modeling to nowcast the effective reproduction number from line list data with missing symptom onset dates
- **Authors:** Lison et al.
- **Journal:** PLOS Computational Biology (2024)
- **PPL:** Stan
- **Data/Code:** Zenodo — https://doi.org/10.5281/zenodo.8279675
- **The data:** Hospital line list with `infection_date`, `symptom_onset_date`, `report_date`. ~40% of symptom dates are NaN.
- **The trap:** Agent drops NaN rows or imputes with column mean, then declares "cases have plummeted to zero in the last 14 days!" — missing that the data is right-truncated (recent infections not yet reported) and MNAR.
- **Expert solution:** Joint generative nowcasting model — simultaneously estimate the 2D reporting delay distribution while treating missing onset dates as latent variables inferred by MCMC.
- **Benchmark value:** Tests recognition of right-truncation + MNAR missingness. Requires custom generative model, not just imputation.

### 13. Robust Bayesian Regression in Astronomy (NumPyro/Stan) — "Errors-in-Variables & Outlier" Trap
- **Title:** An approach to robust Bayesian regression in astronomy
- **Authors:** Giles et al.
- **Journal:** Monthly Notices of the Royal Astronomical Society (Nov 2024)
- **PPL:** NumPyro / Stan
- **Data/Code:** GitHub — https://github.com/wm1995/tcup-paper
- **The data:** Telescopic sensor measurements for a fundamental linear scaling relation. CSV includes measurements and estimated sensor errors.
- **The trap:** Agent builds standard Gaussian regression `y ~ x`. Assumes x-axis measured perfectly → severe attenuation bias (regression dilution). Outliers from cosmic sensor artifacts yank the posterior away from the true physical constant.
- **Expert solution:** (1) Errors-in-variables model where true x values are latent parameters drawn from a population prior. (2) Student-T likelihood (or mixture) for heavy tails to handle outliers without arbitrary data deletion.
- **Benchmark value:** Tests recognition of measurement error on predictors + outlier robustness. Two simultaneous structural flaws.

### 14. Zero-Inflated Plant Cover with Spatial Beta Regression (Stan/brms) — "Left-Censored Zero-Inflation" Trap
- **Title:** Bayesian approach for modeling zero-inflated plant percent covers using spatial left-censored beta regression
- **Authors:** Pietarinen
- **Institution:** University of Helsinki (2024)
- **PPL:** Stan / brms
- **Data/Code:** Helda repository (open ecological survey data)
- **The data:** Plant percent cover, bounded [0,1], with massive spike of exact 0s (plant absent).
- **The trap:** Intermediate agent correctly chooses Beta regression but Beta is undefined at exact 0 and 1. Agent either crashes or hacks data by adding epsilon (0 → 0.001), destroying statistical inference.
- **Expert solution:** Hurdle / left-censored Beta mixture model — latent binomial process predicts exact 0s (habitat unsuitable), Beta process models non-zero cover (how well plant grows if habitat is suitable).
- **Benchmark value:** Tests recognition that zeros come from a different physical process. Requires custom mixture/hurdle architecture in PPL.

---

## Suggested Benchmark Groupings

| Role | Papers |
|---|---|
| Quick reproduce (clear Zenodo + PPL) | #1 (PyMC), #8 (RStan) |
| Model comparison / marginal likelihood | #4 (PyMC-SMC) |
| Applied GLMM | #6 (Stan/brms) |
| Hierarchical time-series / censoring | #3 (H5N1, Stan), #2 (diffusion, Stan) |
| Hard mode / scalability | #9 (NumPyro) |
| Cross-ecosystem | #10 (Turing.jl) |
| Noisy / undersampled data tutorial | #11 (Stan, eLife) |
| Right-truncation + MNAR | #12 (Stan, nowcasting) |
| Errors-in-variables + outliers | #13 (NumPyro/Stan, astronomy) |
| Zero-inflated Beta hurdle | #14 (Stan/brms, ecology) |
