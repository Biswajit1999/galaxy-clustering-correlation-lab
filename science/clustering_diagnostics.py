"""Pair-count and BAO diagnostics for *synthetic* or user-supplied catalogs.

Distances and separations are in comoving Mpc/h. The browser's Gaussian
BAO feature is illustrative, not a direct determination of H0.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


def _catalog(points):
    a = np.asarray(points, dtype=float)
    if a.ndim != 2 or a.shape[1] != 3 or len(a) < 2 or not np.all(np.isfinite(a)):
        raise ValueError("catalog must be a finite (N>=2, 3) array in Mpc/h")
    return a


def pair_histogram(left, right, bin_edges, autocorrelation=False):
    """Count unordered same-catalog pairs or all cross-catalog pairs.

    Uses Euclidean separations in the supplied Cartesian comoving frame;
    no periodic boundary or angular mask correction is implied.
    """
    a, b = _catalog(left), _catalog(right)
    bins = np.asarray(bin_edges, dtype=float)
    if bins.ndim != 1 or len(bins) < 2 or not np.all(np.isfinite(bins)) or np.any(np.diff(bins) <= 0) or bins[0] < 0:
        raise ValueError("bin edges must be increasing, finite, nonnegative")
    if autocorrelation:
        if a.shape != b.shape or not np.array_equal(a, b):
            raise ValueError("autocorrelation requires the same catalog")
        d = np.linalg.norm(a[:, None, :] - a[None, :, :], axis=-1)
        distances = d[np.triu_indices(len(a), k=1)]
    else:
        distances = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1).ravel()
    return np.histogram(distances, bins=bins)[0]


def landy_szalay(data, random_catalog, bin_edges):
    """Landy-Szalay xi=(DD-2DR+RR)/RR with normalized pair counts.

    RR=0 bins return NaN, never a misleading zero; finite-volume and
    survey-selection effects require an appropriate random catalog.
    """
    d, r = _catalog(data), _catalog(random_catalog)
    dd = pair_histogram(d, d, bin_edges, True) / (len(d) * (len(d) - 1) / 2)
    dr = pair_histogram(d, r, bin_edges) / (len(d) * len(r))
    rr = pair_histogram(r, r, bin_edges, True) / (len(r) * (len(r) - 1) / 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        xi = (dd - 2 * dr + rr) / rr
    xi[rr == 0] = np.nan
    return {"xi": xi, "DD": dd, "DR": dr, "RR": rr}


def region_jackknife(data, random_catalog, data_regions, random_regions, bin_edges):
    """Delete-one-region covariance: (K-1)/K * sum((xi_k-mean)^2).

    Each region must retain >=2 data and random objects in every replicate;
    all separation bins must have nonzero RR in every replicate.
    """
    d, r = _catalog(data), _catalog(random_catalog)
    dr = np.asarray(data_regions)
    rr = np.asarray(random_regions)
    if dr.shape != (len(d),) or rr.shape != (len(r),):
        raise ValueError("region labels must match catalog lengths")
    regions = np.union1d(dr, rr)
    if len(regions) < 3 or not np.array_equal(np.unique(dr), regions) or not np.array_equal(np.unique(rr), regions):
        raise ValueError("at least three regions must occur in both catalogs")
    replicates = []
    for k in regions:
        sub_d, sub_r = d[dr != k], r[rr != k]
        if min(len(sub_d), len(sub_r)) < 2:
            raise ValueError("insufficient objects after deleting a region")
        xi = landy_szalay(sub_d, sub_r, bin_edges)["xi"]
        if not np.all(np.isfinite(xi)):
            raise ValueError("zero RR bin in a jackknife replicate; rebin or add randoms")
        replicates.append(xi)
    values = np.asarray(replicates)
    mean = values.mean(axis=0)
    residual = values - mean
    covariance = (len(regions) - 1) / len(regions) * (residual.T @ residual)
    return {"replicates": values, "mean": mean, "covariance": covariance}


@dataclass(frozen=True)
class BAOAmplitude:
    amplitude: float
    sigma_amplitude: float
    chi2: float
    delta_chi2_against_zero: float


def fit_bao_amplitude(separation, xi, covariance, smooth_model,
                      center=105.0, width=10.0) -> BAOAmplitude:
    """GLS amplitude for a *fixed-center, fixed-width* Gaussian BAO template.

    No cosmological distance/H0 inference: both the coordinate conversion
    and the sound horizon are model-dependent. Errors are conditional on the
    supplied smooth baseline and covariance, not a detection significance.
    """
    x = np.asarray(separation, dtype=float)
    y = np.asarray(xi, dtype=float)
    baseline = np.asarray(smooth_model, dtype=float)
    cov = np.asarray(covariance, dtype=float)
    if x.ndim != 1 or x.size < 3 or y.shape != x.shape or baseline.shape != x.shape:
        raise ValueError("separation, xi and smooth model must have equal 1D lengths >=3")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)) or not np.all(np.isfinite(baseline)) or np.any(x <= 0):
        raise ValueError("finite positive separations and finite models required")
    if cov.shape != (len(x), len(x)) or not np.all(np.isfinite(cov)) or not np.allclose(cov, cov.T):
        raise ValueError("covariance must be finite symmetric NxN")
    if not np.isfinite(center) or not np.isfinite(width) or width <= 0:
        raise ValueError("invalid template center or width")
    template = np.exp(-.5 * ((x - center) / width)**2)
    residual = y - baseline
    try:
        weighted_template = np.linalg.solve(cov, template)
        weighted_residual = np.linalg.solve(cov, residual)
    except np.linalg.LinAlgError as exc:
        raise ValueError("covariance must be invertible and positive definite") from exc
    if np.linalg.eigvalsh(cov).min() <= 0:
        raise ValueError("covariance must be positive definite")
    information = float(template @ weighted_template)
    if information <= 0:
        raise ValueError("template has zero information")
    amplitude = float(template @ weighted_residual) / information
    remaining = residual - amplitude * template
    chi2 = float(remaining @ np.linalg.solve(cov, remaining))
    return BAOAmplitude(amplitude, float(1 / np.sqrt(information)), chi2,
                        float(amplitude**2 * information))
