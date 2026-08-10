# Galaxy Clustering Correlation Lab

Two-point correlation functions, galaxy bias and BAO-scale structure, in a zero-build browser lab.

Created and maintained by Biswajit Jana.

## Scientific Background

### Large-scale structure and the two-point correlation function

Galaxies are not distributed randomly through space. Gravitational instability acting on
primordial density fluctuations grows a "cosmic web" of filaments, sheets and voids, and
galaxies trace the underlying matter density field. The standard statistical tool for
quantifying this clustering is the **two-point correlation function**, xi(r), introduced in
its modern form by Peebles (1980). It measures the excess probability, relative to a random
(Poisson) distribution, of finding a galaxy pair separated by comoving distance r:

```
dP = n [1 + xi(r)] dV
```

where n is the mean galaxy number density and dV is a volume element at separation r from a
given galaxy. xi(r) > 0 means galaxies cluster more than random; xi(r) = 0 means no excess
clustering; xi(r) < 0 means galaxies are anti-correlated (under-dense) at that scale.

On small to intermediate scales (roughly 1-20 Mpc/h), xi(r) is well described empirically by a
power law,

```
xi(r) = (r / r0)^(-gamma)
```

with a correlation length r0 (typically ~5 Mpc/h for optical galaxy samples) and slope gamma
(typically ~1.6-1.8). Because galaxies are a biased tracer of the underlying dark-matter field,
the observed amplitude is scaled by the **linear bias** b, so that the galaxy correlation
function relates to the matter correlation function via xi_gal(r) ~ b^2 * xi_matter(r).

### Baryon Acoustic Oscillations (BAO) as a standard ruler

Before recombination, the photon-baryon plasma in the early universe supported acoustic waves
sourced by the same primordial perturbations that seeded structure formation. When
recombination occurred (z ~ 1100), the photons decoupled from the baryons and these sound
waves "froze in," leaving a faint but statistically detectable overdensity shell in the later
matter distribution at the comoving sound horizon scale, r_d ~ 150 Mpc. This appears in xi(r) as
a small, broad bump near r ~ 100-110 Mpc/h -- the **BAO feature** -- superimposed on the
smooth power-law clustering signal.

Because the sound-horizon scale can be computed precisely from early-universe physics (the
baryon and photon densities, and the expansion history up to recombination), the BAO bump acts
as a **standard ruler**: measuring its apparent size at different redshifts calibrates the
cosmic distance scale independently of supernova-based standard candles. Eisenstein et al.
(2005) made the first clear detection of this peak in the correlation function of SDSS luminous
red galaxies, and it has since become one of the pillars of precision cosmology.

### Connection to LCDM and dark energy

BAO measurements constrain the angular diameter distance D_A(z) and the Hubble parameter H(z)
through the ratio of the sound horizon to the observed peak scale. Combined with the cosmic
microwave background and Type Ia supernovae, BAO distance measurements over a range of
redshifts are one of the main probes used to constrain the matter density Omega_m, the
curvature of space, and the equation of state of dark energy in flat LCDM and its extensions.
Large spectroscopic surveys -- SDSS, BOSS, eBOSS and now DESI -- exist largely to map galaxy
positions over ever larger cosmological volumes so that xi(r) (and its Fourier-space
counterpart, the power spectrum P(k)) can be measured with enough precision to resolve the BAO
feature and track its evolution with redshift.

## How It Works

This is a client-side, dependency-free simulation and reference-comparison tool, not a fit to
real survey catalogs. It has two layers that are shown together:

1. **Reference anchors** (`data/reference.json`): four representative xi(r) data points --
   one-halo scale, correlation length, two-halo scale and the BAO bump -- labelled and plotted
   as fixed markers so the model curve can be visually compared against literature-motivated
   benchmarks (citing Eisenstein et al. 2005).
2. **Adjustable model** (`physicsWorker.js`, `clustering()` function): given the current slider
   parameters, a Web Worker evaluates

   ```
   xi(r) = bias^2 * (r / r0)^(-gamma) * exp(-r / 65) + bao * exp(-0.5 * ((r - 105) / 10)^2)
   ```

   over r in [0.5, 150] Mpc/h. This is a power-law clustering term with a soft large-scale
   cutoff, plus a Gaussian bump centered at r = 105 Mpc/h representing the BAO feature. The
   worker also produces a normalized 2D heatmap ("projected pair density") as an illustrative
   visualization of clustering strength on the sky-projected plane, and summary metrics
   (xi at r = 5 Mpc/h, the BAO contrast, and the bias value).

The model runs entirely in a Web Worker (`physicsWorker.js`) so slider interaction never blocks
the UI thread. `app.js` owns UI state, builds the parameter sliders from the `LAB.controls`
definition, posts parameter updates to the worker, and draws both the xi(r) curve (with the
reference points overlaid) and the heatmap on `<canvas>` elements. `styles.css` provides the
dark "mission-control" dashboard styling.

### Adjustable parameters

| Parameter | Meaning | Range |
|---|---|---|
| `r0` | Correlation length [Mpc/h] | 2.5 - 9 |
| `gamma` | Small-scale power-law slope | 1.2 - 2.4 |
| `bias` | Linear galaxy bias | 0.8 - 3 |
| `bao` | BAO bump contrast (amplitude) | 0 - 0.08 |

## Usage

```bash
python -m http.server 8080
```

Open `http://localhost:8080` and move the sliders. The plot updates live: the model curve
(green) is drawn against the fixed literature reference points (yellow markers), and the
heatmap panel shows the corresponding projected pair-density field.

## Validate

```bash
npm run check
```

`scripts/validate.js` checks that required files exist, that `data/reference.json` is well
formed, that `physicsWorker.js` parses, that citations are present, and that no unfinished
scaffold tokens remain. `scripts/validate_repository.mjs` performs the extended
research-quality checks described in [RESEARCH_QUALITY.md](RESEARCH_QUALITY.md), validating
`data/research-reference.json` and the optional `research-overlay.js` mission-control panel.

## Math Appendix

**Two-point correlation function (definition):**

```
dP(pair) = n^2 [1 + xi(r)] dV1 dV2
```

**Power-law clustering (small/intermediate scales):**

```
xi(r) = (r / r0)^(-gamma)
```

**Galaxy bias:**

```
xi_galaxy(r) ~ b^2 * xi_matter(r)
```

**Model implemented in this repository:**

```
xi(r) = bias^2 * (r / r0)^(-gamma) * exp(-r / 65) + bao * exp(-0.5 * ((r - 105) / 10)^2)
```

**BAO standard ruler:** the comoving sound horizon at the drag epoch, r_d ~ 150 Mpc, sets the
characteristic scale of the correlation-function bump; comparing the apparent angular/redshift
size of this feature across redshift constrains D_A(z) and H(z) in the assumed cosmology.

## Architecture

- `index.html`: mission-control interface.
- `styles.css`: dense dark scientific dashboard.
- `app.js`: UI state, Canvas rendering and worker orchestration.
- `physicsWorker.js`: numerical model and heatmap generation.
- `data/reference.json`: small auditable reference-data bundle.
- `data/research-reference.json`: extended benchmark anchors for the validation layer.
- `research-overlay.js`: optional mission-control quality panel (validation status, telemetry).
- `scripts/validate.js`: no-dependency repository validation.
- `scripts/validate_repository.mjs`: extended research-quality validation.

## References

- Peebles, P.J.E., 1980. The Large-Scale Structure of the Universe. Princeton University Press.
- Eisenstein, D.J. et al., 2005. Detection of the baryon acoustic peak in the large-scale
  correlation function of SDSS luminous red galaxies. The Astrophysical Journal, 633(2),
  pp.560-574.
- Wilson, G. et al., 2017. Good enough practices in scientific computing. PLOS Computational
  Biology, 13(6), p.e1005510.

## Research Quality Upgrade

See [RESEARCH_QUALITY.md](RESEARCH_QUALITY.md) for the validation layer, reference anchors,
equations and research boundaries added to this repository.
