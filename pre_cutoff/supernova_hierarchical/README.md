# 08 — Hierarchical Supernova Dust Inference (NumPyro/JAX)

## Citation

Grayling M, Thorp S, Uzsoy AS, Narayan G, Mandel KS (2024).
**Scalable hierarchical BayeSN inference: Investigating dependence of SN Ia host galaxy dust properties on stellar mass and redshift.**
*Monthly Notices of the Royal Astronomical Society* (accepted). arXiv:2401.08755.

## Links

| Resource | URL |
|----------|-----|
| Paper (arXiv) | https://arxiv.org/abs/2401.08755 |
| GitHub (BayeSN) | https://github.com/bayesn/bayesn |
| PyPI | `pip install bayesn` |
| Documentation | https://bayesn.readthedocs.io/en/stable/ |
| Foundation DR1 data | https://github.com/djones1040/Foundation_DR1 |
| Pantheon+ data | https://github.com/PantheonPlusSH0ES/DataRelease |

## Folder structure

```
08_supernova_hierarchical/
├── README.md                              (this file)
├── paper/
│   └── 2401.08755v1.pdf                   (2.8 MB — arXiv preprint)
├── data/
│   └── example_lcs/
│       └── Foundation_DR1_2016W.txt       (SNANA-format example light curve)
└── code/
    ├── REPO_README.md                     (original repo README)
    ├── bayesn_model.py                    (3,524 lines — all NumPyro models)
    ├── bayesn_io.py                       (132 lines — SNANA data I/O)
    ├── spline_utils.py                    (174 lines — JAX-compatible spline interpolation)
    ├── zltn_utils.py                      (242 lines — zero-lower-truncated normal utilities)
    └── model_files/
        ├── T21_model/                     (Thorp+2021 pre-trained SED model)
        │   ├── l_knots.txt, tau_knots.txt (wavelength/phase knot locations)
        │   ├── W0.txt, W1.txt             (mean SED + 1st principal component)
        │   ├── L_Sigma_epsilon.txt        (Cholesky factor of residual covariance)
        │   └── M0_sigma0_RV_tauA.txt     (global calibration parameters)
        └── W22_model/                     (Ward+2023 pre-trained SED model)
            └── (same 6 files)
```

## Study

Type Ia supernovae (SNe Ia) are "standardisable candles" — their intrinsic brightness variations can be corrected using light curve shape and color, enabling distance measurements for cosmology. A key source of color variation is **dust** in the host galaxy, characterised by extinction AV and the reddening law parameter RV.

This paper uses 475 SNe Ia from three surveys (Foundation DR1: 157, DES 3-year: 119, Pan-STARRS: 199) to ask: do dust properties depend on host galaxy stellar mass and/or redshift? The answer matters because systematic biases in dust correction propagate directly into cosmological parameter estimates (the Hubble constant, dark energy equation of state).

## The modelling challenge: physics-forward hierarchical inference

### SED model (pre-trained, frozen)

Each supernova's spectral energy distribution (SED) at rest-frame wavelength λ and phase t is modelled as:

```
S(λ, t) = M₀ · Hsiao(λ, t) · 10^{-0.4 · [W₀(λ,t) + θ·W₁(λ,t) + ε(λ,t)]}
```

where Hsiao is a template spectrum, W₀ and W₁ are cubic spline surfaces (the "warping" matrices trained on data), θ is a stretch-like parameter, and ε captures residual intrinsic variation. The pre-trained model parameters (W₀, W₁, L_Σ) are frozen during dust inference.

### Dust extinction

The SED is reddened by the Fitzpatrick99 dust law:

```
S_obs(λ) = S(λ) · 10^{-0.4 · AV · [k(λ; RV) / RV + 1]}
```

where AV is the V-band extinction and RV controls the shape of the extinction curve (how "gray" vs. "chromatic" the dust is).

### Observer-frame flux

The model SED is integrated through photometric filter transmission curves to produce synthetic photometry, then scaled by distance modulus:

```
flux = ∫ S_obs(λ/(1+z)) · T_band(λ) · dλ / ∫ T_band(λ) · dλ · 10^{-0.4·Ds}
```

### Hierarchical dust population model

**Population-level hyperparameters:**

| Parameter | Prior | Description |
|-----------|-------|-------------|
| `mu_R` | Uniform(1.2, 6) | Mean of RV population distribution |
| `sigma_R` | HalfNormal(2) | Spread of RV across SNe |
| `tauA` | Half-Cauchy (via `tan(Uniform(0, π/2))`) | Scale of exponential AV distribution |
| `sigma0` | 0.1 · tan(Uniform(0, π/2)) | Intrinsic achromatic scatter |

**Individual SN latent parameters (per SN, within `numpyro.plate`):**

| Parameter | Prior | Dimension | Description |
|-----------|-------|-----------|-------------|
| `theta` | Normal(0, 1) | 1 | Shape/stretch |
| `AV` | Exponential(1/tauA) | 1 | Dust extinction |
| `RV` | TruncNormal(mu_R, sigma_R; ≥1.2) | 1 | Dust law parameter |
| `eps_tform` | MvNormal(0, I) | ~24 | Intrinsic SED residuals (non-centred) |
| `Ds` | Normal(muhat, Ds_err) | 1 | Distance modulus |

**Total latent dimensions:** ~29 per SN × 475 SNe = **~13,800 individual parameters** + 4 hyperparameters.

### Truncated normal reparameterisation for RV

RV must be ≥ 1.2 (physical constraint). Instead of sampling from a truncated distribution directly, the code uses the CDF-inverse-CDF trick:

```python
phi_alpha_R = norm.cdf((1.2 - mu_R) / sigma_R)
Rv_tform ~ Uniform(0, 1)
Rv = mu_R + sigma_R * ndtri(phi_alpha_R + Rv_tform * (1 - phi_alpha_R))
```

This maps a uniform sample through the inverse CDF of the truncated Gaussian, enabling smooth gradients for HMC.

### Model variants

| Model | Hyperparameters | Key difference |
|-------|----------------|----------------|
| `dust_model` | 4 | Single population, no mass/redshift dependence |
| `dust_model_split_mag` | 10 | Separate dust params for high/low host mass, + magnitude step |
| `dust_model_split_sed` | 8 + SED perturbations | Mass-split with intrinsic SED differences (ΔW₀) |
| `dust_redshift_model` | 6 | Linear redshift evolution of mu_R and tauA |

## Implementation details

- **PPL**: NumPyro (JAX-based)
- **Sampler**: NUTS, 4 chains in parallel
- **Hardware**: 4 NVIDIA A100 GPUs (paper); CPU fallback supported
- **Runtime**: ~15-30 minutes on A100 for 475 SNe
- **Package**: `pip install bayesn` (MIT license, Python ≥ 3.10)
- **Dependencies**: jax, numpyro, sncosmo, astropy, extinction, h5py, arviz

### Data assembly note

The specific 475-SN sample is not bundled as a ready-made file. Reproducing the paper requires:
1. Downloading SNANA light curves from Foundation DR1 and Pantheon+ DataRelease
2. Applying selection cuts (redshift, light curve quality, host mass availability)
3. Constructing a data table with SNID, peak MJD, redshift, host stellar mass, file paths

For benchmarking, simpler alternatives exist:
- Use Foundation DR1 alone (157 SNe, single survey, no cross-calibration issues)
- Use BayeSN's simulation tools to generate synthetic data with known parameters
- Fit a single SN as a minimal end-to-end test

## Benchmark value

### "Hard mode" hierarchical inference at scale

This is the most computationally demanding benchmark in the collection. With ~13,800 latent parameters and a physics-forward likelihood (SED computation → dust extinction → filter integration → flux), it pushes beyond what simple PPL usage can handle.

### Agent evaluation dimensions

1. **NumPyro/JAX proficiency** — can the agent work with `numpyro.plate`, `numpyro.handlers.mask`, JAX's `vmap`/`jit`, and the NumPyro NUTS API? This is a different ecosystem from Stan and PyMC.

2. **Staged inference** — the model conditions on pre-trained SED parameters (W₀, W₁, L_Σ) that are loaded from files and frozen. The agent must understand this "train once, then infer dust" pipeline.

3. **Reparameterisation techniques** — tan-transformed half-Cauchy priors, CDF-inverse-CDF for truncated normals, Cholesky decomposition for correlated residuals. All essential for efficient HMC.

4. **Physics-in-the-loop** — the likelihood involves real astrophysics: SED warping, Fitzpatrick99 dust extinction, photometric filter integration, cosmological distance. An agent must either understand or correctly delegate this physics.

5. **GPU-scale inference** — with 475 SNe × 29 latent parameters, this requires GPU acceleration. An agent must handle JAX device placement and memory management.

6. **Model extension** — the paper's key contribution is extending the base dust model to account for host mass and redshift dependence. Can the agent implement these extensions from a description?

### Comparison to other benchmark papers

| Dimension | This paper | Paper 04 (ODE) | Paper 07 (Pain) |
|-----------|-----------|-----------------|------------------|
| PPL | NumPyro (JAX) | PyMC | RStan |
| Scale | ~13,800 latent params | ~30 per model | ~135 per condition |
| GPU required | Yes (practical) | No | No |
| Physics | SED + dust + cosmology | ODE systems biology | Cognitive RL/KF |
| Hierarchy | 2-level (pop + SN) | None (per-model) | 2-level (group + subject) |
| Model comparison | 4 dust variants | 10 ODE models via BMA | 5 cognitive models via LOO |
| Key challenge | Scale + physics | ODE solver in likelihood | Sequential dynamics |
