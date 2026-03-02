# 05 — Rabies Vaccination Coverage Heterogeneity (Power Mean Models)

## Citation

Ferguson EA, Lugelo A, Czupryna A, Anderson D, Lankester F, Sikana L, Dushoff J, Hampson K (2025).
**Improved effectiveness of vaccination campaigns against rabies by reducing spatial heterogeneity in coverage.**
*PLOS Biology* 23(1): e3002872.

## Links

| Resource | URL |
|----------|-----|
| Paper | https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002872 |
| Zenodo archive | https://doi.org/10.5281/zenodo.15249730 |
| GitHub | https://github.com/ElaineAFerguson/Serengeti_vaccination_heterogeneity |

## Folder structure

```
05_rabies_vaccination/
├── README.md                          (this file)
├── paper/
│   └── journal.pbio.3002872.pdf       (2.3 MB — full paper)
├── data/
│   ├── serengeti_rabid_dogs.csv       (831 KB, 3,362 rows — rabid dog cases with locations/dates)
│   ├── tree_ct_data.csv               (206 KB, 3,489 rows — transmission tree case data)
│   ├── DogdemographydataPlos2016.csv  (1.4 MB, 7,617 rows — dog demography panel data)
│   ├── vaccinationCoverageByVillageMonth_Jan2002_Dec2022.csv
│   │                                  (423 KB, 87 villages × 252 months — monthly village coverage)
│   ├── vcVillByYear.csv               (25 KB, 87 villages × 21 years — annual village coverage)
│   └── vaccination_table.csv          (450 B, 3 rows — district-level annual summary)
└── code/
    ├── REPO_README.md                 (original repo README)
    ├── stan/
    │   ├── power_mean_model_village.stan              (4.7 KB — village-level, unstandardised)
    │   ├── power_mean_model_village_standardised.stan (4.9 KB — village-level, standardised)
    │   ├── power_mean_model_district.stan             (2.4 KB — district-level, unstandardised)
    │   └── power_mean_model_district_standardise.stan (2.5 KB — district-level, standardised)
    └── R/
        ├── workflow.R                 (2.5 KB — master script, sources 25 subscripts)
        ├── incidence_coverage_models.R (77 KB — brms models + diagnostics + plotting)
        ├── incidence_coverage_models_annual.R (46 KB — annual-resolution brms models)
        ├── model_data_prep.R          (19 KB — monthly data preparation, spatial weights)
        ├── model_data_prep_annual.R   (17 KB — annual data preparation)
        ├── Power_mean_model_village_stan.R  (12 KB — fits village Stan models, WAIC)
        ├── Power_mean_model_district_stan.R (5.9 KB — fits district Stan models)
        ├── Plot_stan_district_monthly.R (8.6 KB — district model plots)
        ├── Plot_stan_village_monthly.R  (33 KB — village model plots)
        ├── Simulate_stan_village_monthly.R (37 KB — counterfactual simulations)
        ├── explore_powers.R           (2.2 KB — explores power parameter behaviour)
        ├── process_CT_data.R          (13 KB — contact tracing data processing)
        └── Functions/
            ├── Lmeans.R               (power mean implementation: powTrans, invTrans, powMean)
            └── computePeriod.R        (utility for time period calculation)
```

## Study

Rabies vaccination campaigns in Serengeti District, Tanzania (2002--2022). The district contains ~87 villages with annual mass dog vaccination campaigns. The core question: does *spatial heterogeneity* in vaccination coverage matter beyond the district-level average? If 70% of dogs are vaccinated overall but all in the same 5 villages, the district is far less protected than if coverage is evenly spread.

The data span 21 years of campaigns, 3,362 confirmed rabid dog cases, and monthly records of dog populations and vaccination coverage across all villages.

## Data format

### `serengeti_rabid_dogs.csv` — Rabid dog case records (3,362 rows)

| Column | Description |
|--------|-------------|
| `District`, `Ward`, `Village` | Administrative location |
| `UTM.Easting`, `UTM.Northing` | Spatial coordinates (jittered ±1km for de-identification) |
| `ID`, `Biter.ID` | Case ID and ID of biting animal (0 = unknown/incursion) |
| `Species` | Almost all "Domestic dog" |
| `Rabid` | Confirmed rabid (TRUE/FALSE) |
| `Vaccination` | Vaccination status of the case animal |
| `Date.bitten`, `Symptoms.started` | Epidemiological dates |
| `Incubation.period`, `Infectious.period` | Clinical periods |
| `month` | Month index (1 = Jan 2002, ..., 252 = Dec 2022) |

### `vaccinationCoverageByVillageMonth_Jan2002_Dec2022.csv` — Monthly village coverage (87 × 252)

Matrix of vaccination coverage proportions (0--1). Rows = villages, columns = months (Jan 2002 through Dec 2022). This is the key predictor variable — it captures the *spatial heterogeneity* that motivates the power mean approach.

### `vcVillByYear.csv` — Annual village coverage (87 × 21)

Same structure but aggregated annually (2002--2022). Row per village, column per year.

### `vaccination_table.csv` — District summary (3 rows)

| Row | Description |
|-----|-------------|
| `n_dogs_vax_dist` | Total dogs vaccinated in district per year |
| `percent_dogs_vax` | District-level vaccination coverage (%) |
| `percent_villages_vax` | Percentage of villages receiving any vaccination |

### `DogdemographydataPlos2016.csv` — Dog demography panel (7,617 rows)

Longitudinal dog census from household surveys. Each row = one dog × one survey round. Key columns: `Village`, `dogid`, `sex`, `birthdate`, `deathdate`, `deathcause`, `alive`, `dead`, `age`, `lagvacc` (vaccinated at last round), survival hazard components (`sch1`--`sch17`, `sca1`--`sca17`).

### `tree_ct_data.csv` — Transmission tree data (3,489 rows)

Compact version of case data for transmission tree reconstruction: `id_case`, `id_biter`, `x_coord`, `y_coord`, `owned`, `date_symptoms`, `days_uncertain`.

## The modelling challenge: quantifying heterogeneity with power means

### The problem with simple averages

Standard vaccination models use district-level mean coverage as a predictor: if 30% of dogs are vaccinated, susceptibility = 0.70. But this ignores *where* those dogs are. Consider two scenarios both averaging 30%:

- **Uniform**: every village at 30% coverage → susceptibility ≈ 0.70 everywhere
- **Heterogeneous**: half the villages at 60%, half at 0% → the unprotected villages act as transmission reservoirs

The arithmetic mean treats these identically. The paper's innovation is to replace the simple mean with a **power mean** (also called generalised mean) that can upweight extreme values.

### Power mean definition

For susceptibilities s_1, ..., s_n with weights w_1, ..., w_n:

```
M_p(s) = ( Σ w_i s_i^p / Σ w_i )^{1/p}
```

Special cases:
- p = 1: arithmetic mean (heterogeneity doesn't matter)
- p = 2: root-mean-square (penalises variance)
- p → ∞: maximum (only the worst village matters)
- p → 0: geometric mean (requires special handling)

The key insight: **p is estimated from data**, not chosen a priori. If the posterior for p concentrates well above 1, that's direct evidence that spatial heterogeneity in coverage impacts disease transmission.

### Numerical stability near p = 0

The raw formula `(x^p - 1)/p` is undefined at p = 0. The code uses a smooth first-order Taylor expansion when |p| < ε:

```
powTrans(x, p) ≈ log(x) + p·log(x)²/2     when |p| < ε
invTrans(y, p) ≈ exp(y - p·y²/2)           when |p| < ε
```

This allows the sampler to smoothly traverse p = 0 without discontinuities. In Stan, this is implemented as a conditional branch in `transformed parameters`.

Additionally, the implementation normalises susceptibilities by their geometric mean before computing the power mean (divides by `exp(mean(log(x)))`), which keeps intermediate values near 1 and prevents numerical overflow for large p.

## Bayesian models

### Two-tier approach: brms baseline → custom Stan

The paper uses a **dual modelling strategy**:

1. **brms models** — standard negative binomial GLMMs using `brm()` with simple vaccination coverage as a linear predictor. These are the baseline models that an analyst comfortable with R formula syntax would naturally write.

2. **Custom Stan models** — the power mean susceptibility models that cannot be expressed as brms formulas, because the power parameter `p` enters nonlinearly into the computation of the predictor variable itself.

### brms district model (baseline)

```r
brm(formula = cases ~ vax_last2monthMean +
                      log_case_rate_last2monthMean +
                      log_dog_density +
                      offset(log(dogs)),
    data    = data_dist,
    family  = negbinomial(),
    prior   = c(set_prior("normal(0,100000)", class = "Intercept"),
                set_prior("normal(0,100000)", class = "b")),
    warmup = 1500, iter = 3000, chains = 4,
    control = list(adapt_delta = 0.99, max_treedepth = 15))
```

This treats mean vaccination coverage as a linear effect. The question is: can we do better by accounting for spatial heterogeneity?

### brms village model

Adds village-level random intercepts `(1|village)` and separates spatial effects into neighbour/non-neighbour case rates. Still limited to *linear* effects of mean coverage.

### Stan power mean model (the innovation)

The custom Stan model replaces the simple `vax_last2monthMean` with power mean susceptibility computed inside `transformed parameters`:

**Likelihood:**
```
Y_vt ~ NegBinomial2(μ_vt, φ)
log(μ_vt) = X_vt · β + γ_v + log(dogs_vt)
```

where `X_vt` includes the power mean susceptibility of neighbouring and non-neighbouring villages, computed jointly with parameter `p`.

**Power mean in transformed parameters (village-level model):**

The Stan code computes separate power means for:
- **Neighbouring villages** (shared borders) — weighted by border length × dog population
- **Non-neighbouring villages** — weighted by dog population only

For each village v at month t:
```
PM_n[v,t] = ( Σ_{j∈neighbours(v)} w_j · s_j^p / Σ w_j )^{1/p} · norm_n[v,t]
```

where `norm_n` is the geometric-mean normalisation constant, and `s_j = S_n[j,t]` is the normalised susceptibility.

The two-month rolling average is then computed:
```
PPM_n[v,m] = (PM_n[v,m] + PM_n[v,m+1]) / 2
```

**Parameters:**

| Parameter | Prior | Description |
|-----------|-------|-------------|
| `beta[1..K+2]` | Normal(0, 100000) | Regression coefficients (intercept + covariates + neighbour/non-neighbour power mean susceptibility) |
| `p` | Normal(1, 2) | Power parameter for generalised mean |
| `phi` | Gamma(0.01, 0.01) | Negative binomial dispersion |
| `sigma_village` | Exponential(0.001) | Village random effect SD |
| `gamma[1..Nv]` | Normal(0, 1) | Standardised village random effects (non-centred) |

The non-centred parameterisation `gamma_t = sigma_village * gamma` avoids the funnel geometry that plagues centred random effects in hierarchical models.

**Design matrix structure (village model, K=7 covariates):**
- Intercept
- Mean susceptibility over last 2 months (village-level)
- Log case rate over last 2 months (own village)
- Log case rate in neighbouring villages
- Log case rate in non-neighbouring villages
- Log dog density
- Human-to-dog ratio
- Power mean susceptibility (neighbours) — computed in `transformed parameters`
- Power mean susceptibility (non-neighbours) — computed in `transformed parameters`

### Model variants

Six Stan model fits at each spatial resolution (district/village):
1. **Full model** — all covariates + power mean susceptibility
2. **Without distant cases** — drops neighbour/non-neighbour case rates
3. **Without any case history** — drops all prior incidence

Each fit in both unstandardised and standardised covariate versions, for a total of 12 Stan models. Model comparison via WAIC.

## Implementation details

- **PPL**: brms (R) for baseline; RStan for power mean models
- **Sampler**: NUTS via RStan, 4 chains, 3000 iterations (1500 warmup), `adapt_delta = 0.95`
- **Runtime warning**: "village fits take a day each!!" (per `workflow.R` line 66)
- **Spatial data**: requires shapefiles for village boundaries (not included — GIS data at Zenodo)
- **Preprocessing pipeline**: 25 R scripts orchestrated by `workflow.R`, including dog population estimation, vaccination coverage interpolation, transmission tree reconstruction, spatial weight matrices
- **Data de-identification**: coordinates jittered ±1km, so distance kernels and village-level assignments will differ slightly from published results

## What we provide vs. the full repo

The full Zenodo archive includes shapefiles, pre-computed outputs, and all 25+ scripts. We include only:
- The **data CSVs** needed for the incidence models
- The **4 Stan models** (the core Bayesian contribution)
- The **key R scripts** for model fitting, data prep, and plotting
- Helper functions (`Lmeans.R`, `computePeriod.R`)

Missing from our subset (available at Zenodo):
- GIS shapefiles (`data/GIS/`, `output/SD_vill/`)
- Pre-computed intermediate outputs (`output/*.csv`, `output/*.rds`, `output/*.RData`)
- Scripts for transmission trees, distance kernels, dog population estimation
- Simulation scripts and figures

## Benchmark value

### The structural trap: brms comfort zone vs. custom Stan necessity

An agent asked to "model the effect of vaccination on rabies incidence" would naturally reach for brms:

```r
brm(cases ~ vaccination_coverage + ..., family = negbinomial())
```

This is *correct* but *incomplete*. It assumes vaccination coverage enters linearly — that a 10% increase from 20→30% has the same effect as 30→40%, regardless of how coverage is distributed across villages. The paper shows this misses the key finding: spatial heterogeneity matters (p ≈ 3--5 in the posterior), and the power mean model substantially outperforms the simple linear model.

### Agent evaluation dimensions

1. **Recognising the limitation of standard tools** — does the agent identify that brms cannot express a nonlinear predictor where the nonlinearity parameter is itself estimated?

2. **Custom Stan implementation** — can the agent write a Stan model where a derived quantity (power mean susceptibility) feeds back into the likelihood, with the derivation parameter (`p`) jointly estimated?

3. **Numerical stability** — does the agent handle the p → 0 singularity in the power mean, either via the Taylor expansion trick or an alternative (e.g., log-space computation)?

4. **Hierarchical spatial structure** — can the agent correctly set up the neighbour/non-neighbour weight matrices and compute spatially-weighted power means?

5. **Non-centred parameterisation** — does the agent use `gamma ~ N(0,1); gamma_t = sigma * gamma` for the village random effects, or fall into the centred parameterisation trap?

### Comparison to other benchmark papers

| Dimension | This paper | Paper 02 (DDM) | Paper 03 (H5N1) | Paper 04 (ODE) |
|-----------|-----------|-----------------|------------------|----------------|
| PPL | brms + raw Stan | Stan (CmdStanR) | Stan | PyMC |
| Trap type | Tool limitation → custom Stan | Censoring in likelihood | Data augmentation censoring | ODE-in-likelihood |
| Key challenge | Nonlinear predictor with estimated exponent | `wiener_lccdf` for censored RTs | Latent censored values as parameters | JAX/diffrax ODE inside PyTensor |
| Spatial structure | Yes (87 villages, neighbour weights) | No | No | No |
| Hierarchical | Village random effects | Subject-level | Cow-level | No (per-model fits) |
