# Regional pair-count sensitivity (17 November 2026)

This module provides (1) the Landy–Szalay correlation estimator for **already normalized** DD, DR and RR counts, and (2) delete-one-region jackknife covariance, its eigenvalue spectrum and rank. Pair-count normalization must account for distinct data/random catalogue sizes, self-pairs, weights and the survey selection function **before** calling the estimator. The leave-one-region-out pair counts must be recomputed for each omitted spatial region; simply subtracting one region's internal pairs is incorrect because cross-region pairs also disappear.

The jackknife covariance uses (R−1)/R times the replicate scatter. Region independence, spatial geometry, edge effects and covariance invertibility require separate validation. Rank is at most R−1, so a small number of regions cannot yield a full-rank covariance across many separation bins. No DESI/SDSS data or BAO/H0 estimate is made.

Reproduce with `cd research/2026-11-17 && python -m pytest -q test_pair_count_audit.py`. All tests are synthetic. AI coding assistance was used; Biswajit Jana remains the requested authenticated Git commit author.
