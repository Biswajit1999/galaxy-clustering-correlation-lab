# Synthetic BAO peak location, AR(1) GLS and ruler identifiability

This module reproduces the `clustering()` toy function in `physicsWorker.js`: `ξ(r)=b²(r/r₀)^(-γ)exp(-r/65)+A exp[-(r−105)²/(2×10²)]`. The anchors in `data/reference.json` are literature-motivated **illustrations**, not digitized survey measurements. No SDSS/DESI galaxy catalogue, window function, cosmological covariance or actual BAO likelihood is used.

Three complementary checks are implemented: (1) raw versus broadband-subtracted peak location, including sub-grid quadratic interpolation and grid convergence; (2) profile-likelihood fitting of a Gaussian BAO amplitude and centre using exact stationary AR(1) innovations whitening under a *known* broadband; and (3) the conditional dependence of the algebraic `H₀=100 r_peak/r_d` unit conversion on an externally assumed physical sound horizon. The third is **not an independent Hubble-constant determination**. The browser's fixed 105 Mpc/h bump cannot by itself measure H₀. A survey inference requires measured anisotropic BAO dilation, a fiducial cosmology, the sound horizon, selection effects and covariance.

## Reproduce

From the repository root with Python 3.11+ and NumPy:

```bash
python -m unittest discover -s tests -p 'test_bao_robustness.py' -v
python research/bao_robustness.py
```

## Deterministic results

For the illustrative model, the raw peak is at 104.84743 Mpc/h and the known-broadband-subtracted peak at 105.00000 Mpc/h. Under hypothetical pointwise σ=0.003 and AR(1) ρ=0/0.5/0.8, the noise-free centre remains 105 Mpc/h, while discrete Δχ²≤1 half-widths are 0.1/0.2/0.3 Mpc/h (grid-limited). Changing the **assumed** ruler from 147.09 to 150 Mpc changes the algebraic conditional H₀ from 71.3849 to 70.0 km/s/Mpc for the same fixed 105 Mpc/h bump. These numbers are synthetic stress tests, not observational estimates or confidence intervals with validated frequentist coverage.

The AR(1) model assumes equally spaced, stationary errors. The profile uses the exact known broadband and a fixed Gaussian width; misspecified broadband or real survey covariance could change both centre and uncertainty. The reported Δχ² widths are coarse discrete profile-grid summaries, not calibrated coverage claims.