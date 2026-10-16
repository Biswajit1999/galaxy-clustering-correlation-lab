# Galaxy-clustering scientific diagnostics (Python)

This module adds independent, testable numerical methods to the existing illustrative browser model. It does **not** contain SDSS, BOSS or DESI galaxy catalogs, and no cosmological parameter is fitted from observations.

## Reproduce

```bash
python -m pip install -r requirements-science.txt
python -m unittest discover -s tests -v
```

### Three scientific analyses

1. **Normalized Landy–Szalay estimator** — `landy_szalay(data, random_catalog, bins)` computes normalized `DD`, `DR`, `RR` and `xi=(DD-2DR+RR)/RR` from Cartesian comoving coordinates in Mpc/h. The random catalog must reproduce the survey's selection function and geometry. `RR=0` returns NaN, not an artificial zero. The direct pair-count method is O(N²) in memory/time, suitable for small catalogs only.
2. **Delete-one-region jackknife covariance** — `region_jackknife` removes each region from both data and random catalogs, recomputes the estimator, and forms `(K-1)/K * sum((xi_k-mean)(xi_k-mean)^T)`. Region definitions should be spatially coherent and larger than the scales of interest. Covariance can be rank-deficient for fewer jackknife regions than bins. This is **not** a calibrated survey covariance without mocks/systematic checks.
3. **Fixed-template BAO amplitude GLS** — `fit_bao_amplitude` fits the amplitude of a Gaussian bump with user-supplied center, width, smooth baseline, and full positive-definite covariance. The conditional `sigma_amplitude` and `delta_chi2` ignore look-elsewhere effects, nuisance-model uncertainty, and survey selection effects. A fixed 105 Mpc/h template does **not** determine `H0`.

### Scientific correction: BAO and H0

The earlier browser code displayed `bao_implied_H0=100*(105/147.09)`. This ratio is **not** an independent measurement of the Hubble constant: the conversion between Mpc and Mpc/h already assumes `h`, the bump position is hard-coded rather than measured, and a BAO-only distance inference requires a fiducial cosmology, redshift-space modeling, Alcock–Paczynski dilation parameters, and a calibrated sound horizon or complementary data. The revised browser telemetry omits the misleading number. Real BAO analyses constrain combinations such as `D_M(z)/r_d` and `D_H(z)/r_d`, not a direct division of two hard-coded distances.

Reference: Landy & Szalay (1993), *ApJ* 412, 64; Eisenstein et al. (2005), *ApJ* 633, 560. All tests are deterministic and synthetic.
