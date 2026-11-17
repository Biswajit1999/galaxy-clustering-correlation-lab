"""Conditional two-point pair-count and regional covariance diagnostics.

Counts must be normalized externally to the same pair-count convention.
No catalogue, survey mask or selection function is inferred.
"""
import numpy as np


def landy_szalay(dd, dr, rr):
    """Landy-Szalay xi=(DD-2DR+RR)/RR for normalized pair counts."""
    dd, dr, rr = np.broadcast_arrays(np.asarray(dd, float),
                                      np.asarray(dr, float), np.asarray(rr, float))
    if np.any(~np.isfinite(dd)) or np.any(~np.isfinite(dr)) or np.any(~np.isfinite(rr)):
        raise ValueError("nonfinite pair counts")
    if np.any(dd < 0) or np.any(dr < 0) or np.any(rr <= 0):
        raise ValueError("invalid normalized pair counts")
    return (dd - 2*dr + rr) / rr


def leave_region_out(dd, dr, rr):
    """Evaluate LS on independently supplied leave-one-region-out normalized counts."""
    dd, dr, rr = (np.asarray(a, float) for a in (dd, dr, rr))
    if dd.ndim != 2 or dd.shape != dr.shape or dd.shape != rr.shape or dd.shape[0] < 3:
        raise ValueError("need >=3 regional replicates with matching shapes")
    return landy_szalay(dd, dr, rr)


def jackknife_covariance(replicates):
    """Delete-one-region jackknife covariance, (R-1)/R sum deviations outer products."""
    x = np.asarray(replicates, float)
    if x.ndim != 2 or x.shape[0] < 3 or x.shape[1] < 1 or not np.all(np.isfinite(x)):
        raise ValueError("invalid replicates")
    delta = x - x.mean(axis=0)
    return (len(x)-1)/len(x) * (delta.T @ delta)


def covariance_spectrum(replicates):
    """Return eigenvalues, numerical rank, and condition of regional covariance."""
    c = jackknife_covariance(replicates)
    eig = np.linalg.eigvalsh(c)
    return {"covariance": c, "eigenvalues": eig,
            "rank": int(np.linalg.matrix_rank(c)),
            "condition": float(np.inf if eig[0] <= 0 else eig[-1]/eig[0])}
