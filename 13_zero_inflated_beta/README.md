# 13 — Zero-Inflated Plant Cover with Left-Censored Beta Regression (Stan/brms)

## Citation

Pietila J (2025).
**Bayesian approach for modeling zero-inflated plant percent covers using spatial left-censored beta regression.**
MSc thesis, University of Helsinki. Faculty of Science, Master's Programme in Mathematics and Statistics.

### Reference paper

Tang B, Frye HA, Gelfand AE, Silander JA Jr (2023).
**Zero-Inflated Beta Distribution Regression Modeling.**
*Journal of Agricultural, Biological and Environmental Statistics* 28: 117–137.

## Links

| Resource | URL |
|----------|-----|
| Thesis (Helda) | http://hdl.handle.net/10138/598194 |
| Reference paper (Springer) | https://doi.org/10.1007/s13253-022-00516-z |
| Reference paper (arXiv) | https://arxiv.org/abs/2112.07249 |
| Reference data | https://github.com/beckytang/ZIB_data |
| VELMU data paper | https://doi.org/10.1038/s41597-024-04092-4 |
| VELMU data (Zenodo) | https://doi.org/10.5281/zenodo.13822190 |

## Folder structure

```
13_zero_inflated_beta/
├── README.md                                      (this file)
├── paper/
│   ├── Pietila_Juho_Mastersthesis_2025.pdf       (1.7 MB — MSc thesis)
│   └── Tang_et_al_2023_ZIB.pdf                   (2.7 MB — reference paper)
└── data/
    ├── cfr_zib_data.csv                           (179 rows × 10 cols — Cape Floristic Region)
    └── ZIB_README.md                              (data description)
```

## Study

Species distribution models combine species observations with environmental data to predict where species occur. When recording plants, the natural measure is **percent cover** — a value in [0, 1] representing how much of a quadrat is covered by a species. This data has a fundamental statistical challenge: **massive excess zeros** where the species is absent.

### The Cape Floristic Region data (Tang et al. reference)

179 sampling plots across the Cape Floristic Region of South Africa, recording percent cover for two plant families (Crassulaceae and Restionaceae) alongside environmental covariates.

### The Baltic Sea data (Pietila thesis)

Underwater macroalgae, vascular plants, and invertebrates from the VELMU Finnish inventory programme. 194,000+ spatially explicit observations across 280+ genera. Percent cover assessed by scuba diving quadrats.

## Data format

### `cfr_zib_data.csv` (179 rows × 10 columns)

| Column | Description |
|--------|-------------|
| `plot` | Plot identifier (e.g., BK_1, BK_10) |
| `latitude` | Decimal degrees (South Africa) |
| `longitude` | Decimal degrees |
| `region` | Geographic region (baviaanskloof, ...) |
| `CRASSULACEAE` | Percent cover [0, 1] — many exact 0s |
| `RESTIONACEAE` | Percent cover [0, 1] — many exact 0s |
| `apan_mean_an` | Mean annual pan evaporation |
| `tminave07c` | Average minimum temperature, July |
| `RFL_CONC` | Rainfall concentration index |
| `BioClimMAP30s` | Mean annual precipitation (BioClim) |

## The structural data trap

This is a particularly clean example of a trap where the **naive model compiles and runs but produces wrong results**.

### Step 1: Agent chooses Beta regression (correct intuition)

Percent cover is bounded [0, 1] — the Beta distribution is the natural choice. So far, so good.

### Step 2: Beta crashes on exact zeros

The Beta distribution has support on (0, 1) — it is **undefined at 0 and 1**. With many species having 0% cover at most sites, the likelihood evaluation fails.

### Step 3: Agent hacks the data (the trap)

The common workaround is to add epsilon: `0 → 0.001`, `1 → 0.999`. This lets the model compile and run, but:
- The arbitrary epsilon choice affects inference
- Zeros from "habitat unsuitable" and zeros from "species present but rare" are conflated
- Parameter estimates are biased, especially the precision parameter phi
- Uncertainty is underestimated

### The correct approach: two sources of zeros

Zeros arise from **two distinct ecological processes**:

1. **Structural zeros** (habitat unsuitable): The species cannot grow here regardless of sampling effort. Modelled by a Bernoulli suitability process.
2. **Sampling zeros** (present but not detected / very low cover): The species could grow here but cover is below the detection threshold. Modelled by left-censoring of the latent Beta regression.

## The model

### Left-censored Beta mixture

For observation y_i at site i:

```
z_i ~ Bernoulli(psi_i)          # Is the habitat suitable?
c_i | z_i=1 ~ Beta(mu_i, phi)   # Latent cover if suitable

y_i = 0  if  z_i = 0            # Structural zero (unsuitable)
y_i = 0  if  z_i = 1, c_i < L   # Sampling zero (left-censored)
y_i = c_i if z_i = 1, c_i >= L  # Observed cover
```

The likelihood contribution:
```
P(y_i = 0) = (1 - psi_i) + psi_i * P(c_i < L | mu_i, phi)
P(y_i > 0) = psi_i * f_Beta(y_i | mu_i, phi)
```

### Covariate structure

Both the suitability probability and the cover mean depend on environmental covariates:
```
logit(psi_i) = x_i' * beta_psi
logit(mu_i)  = x_i' * beta_mu
```

### Spatial random effects

The thesis extends the model with spatial random effects to account for spatial autocorrelation:
```
logit(mu_i) = x_i' * beta_mu + w_i
w ~ GP(0, K(s_i, s_j))    # or CAR/ICAR prior
```

### Variable precision

A key finding of the thesis: modelling the Beta precision parameter phi as a function of covariates (rather than constant) is "highly beneficial for predictive ability and essential to make the model catch the patterns." This means:
```
log(phi_i) = x_i' * beta_phi
```

### Key parameters

| Parameter | Prior | Description |
|-----------|-------|-------------|
| `beta_mu` | Normal(0, ...) | Cover regression coefficients |
| `beta_psi` | Normal(0, ...) | Suitability regression coefficients |
| `phi` (or `beta_phi`) | ... | Beta precision (constant or covariate-dependent) |
| `L` | Fixed (near 0) | Left-censoring threshold |
| `w_i` | GP or CAR | Spatial random effects |

## Implementation details

- **PPL**: Stan (likely via brms or custom Stan code — Vanhatalo group)
- **Supervisor**: Jarno Vanhatalo (EnvStat group, University of Helsinki)
- **Thesis code**: Not deposited separately (code is in the PDF appendix or available on request)
- **Reference code**: Tang et al. data at https://github.com/beckytang/ZIB_data
- **VELMU data**: https://doi.org/10.5281/zenodo.13822190 (194,000+ observations)
- **License**: CC BY 4.0 (thesis)

## Benchmark value

### The "epsilon hack" detector

This paper tests whether an agent falls into the most common trap in Beta regression: adding epsilon to zeros. The correct solution requires understanding that zeros come from a **different physical process** than non-zero observations and building a mixture/hurdle architecture.

### Agent evaluation dimensions

1. **Beta distribution domain awareness** — Can the agent recognise that Beta(a,b) has support (0,1) exclusive, and that exact 0s will cause numerical failures?

2. **Two-process zero recognition** — The key ecological insight: "habitat unsuitable" zeros and "low cover" zeros are different. An agent must implement a mixture model, not just a workaround.

3. **Left-censoring vs. zero-inflation** — This is subtler than standard zero-inflation. The left-censoring component says "the species is present but cover is below detection threshold L." The hurdle/structural component says "the species cannot grow here at all."

4. **Variable precision** — The thesis shows that constant phi is inadequate. An agent should consider modelling precision as a function of covariates, which is non-trivial in brms (requires distributional regression via `bf(y ~ x, phi ~ x)`).

5. **Spatial structure** — The data has spatial coordinates. An agent should consider spatial random effects (GP, CAR, ICAR) but must balance model complexity with computational cost.

### Comparison to other benchmark papers

| Dimension | This paper | Paper 10 (Brain Cells) | Paper 05 (Rabies) |
|-----------|-----------|----------------------|-------------------|
| PPL | Stan/brms | Stan (CmdStanR) | brms + Stan |
| Response | Continuous [0,1] | Count (integer) | Count (integer) |
| Zero problem | Beta undefined at 0 | Structural zeros | Not applicable |
| Solution | Left-censored Beta mixture | Zero-inflated Poisson | Power mean |
| Spatial | GP/CAR random effects | None | Neighbour weights |
| Key insight | Zeros from two processes | Shrinkage via kappa | Custom likelihood |
