# Paper 1: Rat Multisensory Psychophysics (PyMC)

## Citation

Zanzi M, Rinaldi FG, Fornasaro S, Piasini E, Zoccolan D (2025). Seeing what you hear:
Compression of rat visual perceptual space by task-irrelevant sounds. *PLOS Computational
Biology* 21(10): e1013608. https://doi.org/10.1371/journal.pcbi.1013608

## Links

- **Paper:** https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013608
- **PMC:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12571299/
- **Data/Code:** https://zenodo.org/records/17280352 (DOI: 10.5281/zenodo.17280352)

## Folder Structure

```
01_rat_multisensory_psychophysics/
├── README.md
├── paper/
│   └── pcbi.1013608.pdf              # Full paper (35 pages)
├── data/
│   ├── AV_data.csv                   # Visual + AM sound trials (64,227 rows)
│   ├── absent_noise_data.csv         # Visual-only trials (6,873 rows)
│   ├── steady_noise_data.csv         # Visual + fixed-amplitude sound trials (6,925 rows)
│   └── White_noise_sounds/           # WAV audio stimuli for sound calibration
│       ├── 1_gaussian_scaled_1X.wav  # Calibration white noise
│       ├── White_Noise.wav           # Fixed-amplitude stimulus
│       └── noise_*.wav               # AM-modulated stimuli (9 frequencies)
└── code/
    └── Ideal_Observer_analysis_and_figures.ipynb  # Full analysis notebook
```

**Not included locally** (1.3 GB): `Idata_saved.zip` — pre-computed MCMC traces as ArviZ
InferenceData NetCDF files. Available from Zenodo if needed to skip re-sampling.

## Experiment

10 Long Evans rats classify drifting sinusoidal gratings by visual temporal frequency (TF)
into "Low TF" (< 2.12 Hz) vs "High TF" (> 2.12 Hz). Gratings are paired with
task-irrelevant white noise bursts of varying type and intensity:

- **V only** — no sound (visual-only baseline)
- **V + fixed-amplitude sound** — constant white noise (~74.8 dB)
- **V + AM sounds** — amplitude-modulated noise at 9 TFs (0.25–4 Hz), yielding 3 intensity
  clusters (~48.5, ~64.2, ~70.7 dB)

The auditory stimuli are uninformative about the correct visual classification.

## Data Format

All three CSVs share the same columns:

| Column | Description |
|---|---|
| `rat_ID` | Rat identifier (1–10) |
| `Vis_TFs` | Visual temporal frequency (0.25, 0.82, 1.32, 1.75, 2.12, 2.48, 2.91, 3.41, 4 Hz) |
| `Aud_TFs` | Auditory modulation frequency (0 = no sound; 0.25–4 = AM; steady = fixed-amp) |
| `resp_low_TF` | 1 if rat chose "Low TF", else 0 |
| `resp_high_TF` | 1 if rat chose "High TF", else 0 (= `Response` in the model) |
| `RT` | Response time in seconds from trial onset |

Total: ~78,000 trial-level observations.

## Bayesian Model

### Core psychometric equation (Eq. 49 in paper)

```
p(report H | s, n) = ε_H + (1 - ε_H - ε_L) * Φ((γ_n * s - s₀) / σ)
```

- **Φ** — standard Normal CDF (probit link)
- **s** — visual TF of the presented grating
- **s₀** = 2.12 Hz — fixed category boundary
- **γ_n** — perceptual gain factor (depends on sound intensity condition *n*)
- **σ** — sensory noise parameter (controls psychometric slope)
- **ε_H, ε_L** — lapse rates (stimulus-independent "High"/"Low" response probabilities)

### Hierarchical structure

Each rat *r* has its own parameters drawn from population-level distributions:

```
Lapse rates:    (ε_L^r, ε_H^r, 1-ε_L^r-ε_H^r) ~ Dirichlet(1, 1, 1)
Sensory noise:  σ^r ~ Gamma(μ=σ_μ, σ=σ_σ)
Perceptual gain: γ_n^r ~ Normal(μ=γ_μ_n, σ=γ_σ_n)   for each sound condition n
```

Population hyperpriors:
```
σ_μ   ~ Exponential(1/3)
σ_σ   ~ Exponential(1/3)
γ_μ_n ~ Normal(0, 3)        for each sound condition n
γ_σ_n ~ Exponential(1/3)    for each sound condition n
```

Likelihood: `Choice_t ~ Bernoulli(p_t)` per trial.

### Reference model: 5γ 1σ

- 5 γ values (one per sound intensity group), 1 shared σ
- Best model by LOO-ELPD comparison
- Posterior R² = 0.9973 ± 0.0002

### Model variants compared via PSIS-LOO

| Model | Free γ | Free σ | Description |
|---|---|---|---|
| **5γ 1σ** | 5 | 1 | Reference (best) — gain varies by intensity group |
| 5γ 5σ | 5 | 5 | Both gain and noise vary |
| 1γ 5σ | 1 | 5 | Only noise varies |
| 11γ 1σ | 11 | 1 | Gain per AM frequency (overfits) |
| 1γ 1σ | 1 | 1 | Simplest — no sound effect |
| 5γ 1σ λ | 5 | 1 | Adds offset λ parameter (Eq. 47) |
| Linear | 5 | 1 | γ_μ = α + β × measured_intensity |

## Implementation

- **PPL:** PyMC v5.2.0 (pytensor backend)
- **Sampler:** NUTS — 4 chains × 2000 draws (1000 tuning + 1000 posterior), `target_accept=0.99`
- **Diagnostics:** ArviZ 0.15.1 — R-hat, ESS (bulk/tail), MCSE
- **Model comparison:** `arviz.compare()` with PSIS-LOO for ELPD
- **Checks:** Prior predictive (1000 samples), posterior predictive

## Benchmark Value

- **Complete reproducible pipeline** — data + code + pre-computed traces on Zenodo
- **Clean hierarchical Bayesian structure** — probit psychometric function with partial pooling
- **Multiple model variants** for LOO comparison
- **Fully derived math** in paper (Eqs. 10–69) with clear notation
- **Moderate scale** — ~78k Bernoulli obs, 10 rats, ~60 parameters per variant
- Good target for both **reproduction** (run existing code) and **re-implementation** (rebuild from paper)
