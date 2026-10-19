# Offline 3-D pair-count validation (synthetic only)

Run `python -m pytest tests/test_pair_diagnostics.py` and `PYTHONPATH=. python scripts/run_pair_diagnostics.py` from the repository root after installing `numpy` and `pytest`. The existing browser simulation is unchanged. These offline tools are **not** a re-analysis of SDSS, BOSS, or DESI.

Three independently inspectable studies:

1. **Landy–Szalay pair counting** with distinct DD, DR, RR opportunity normalizations, optional periodic minimum-image distances, explicit empty-RR bins, and permutation/full-bin invariants. `xi=(DD-2DR+RR)/RR` uses *normalized* pair counts. A finite cube with uniform randoms reproduces the simple window in the synthetic fixture, not a survey footprint.
2. **Spatial delete-one-cell jackknife** deletes the same cubic subvolume from data and random catalogues; the covariance uses `(K-1)/K * Σ(xi_k-mean)(xi_k-mean)^T`. Its small eight-cell estimate is illustrative and may be rank deficient; no inverse-covariance significance should be inferred.
3. **Random-catalogue sampling sensitivity** holds the data fixed and varies the size/seed of the uniform random catalogue, recording the empirical standard deviation. This does **not** quantify cosmic variance or observational systematics.

Units: positions and bin edges are comoving Mpc/h. No peculiar velocities, redshift-space distortions, selection functions, fiber collisions, survey geometry, shot-noise weighting, or pair-completeness corrections are modeled. Pair counting is O(N²) and intended for small offline checks; a production implementation should use a validated pair-count engine and actual random masks. The browser's illustrative Gaussian BAO bump must not be used to infer H0 from a fixed displayed peak location: ruler calibration and anisotropic distance constraints require real data and a consistent cosmological distance model.

The generated CSV/JSON are seeded synthetic validation artifacts. Do not cite them as observed BAO measurements.
