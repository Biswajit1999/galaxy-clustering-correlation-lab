"""Offline 3-D pair-count diagnostics for *synthetic* galaxy catalogues.

Coordinates and bin edges must share units (e.g. comoving Mpc/h). No redshift
selection, mask, fiber collisions, FKP weights, or survey window are modeled.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class PairResult:
    edges: np.ndarray
    dd: np.ndarray
    dr: np.ndarray
    rr: np.ndarray
    xi: np.ndarray


def _points(a: np.ndarray) -> np.ndarray:
    p = np.asarray(a, dtype=float)
    if p.ndim != 2 or p.shape[1] != 3 or len(p) < 2 or not np.isfinite(p).all():
        raise ValueError("catalogue must be a finite (N,3) array with N>=2")
    return p


def _edges(edges: np.ndarray) -> np.ndarray:
    e = np.asarray(edges, dtype=float)
    if e.ndim != 1 or len(e) < 2 or not np.isfinite(e).all() or e[0] < 0 or np.any(np.diff(e) <= 0):
        raise ValueError("edges must be finite, nonnegative, and strictly increasing")
    return e


def _histogram(a: np.ndarray, b: np.ndarray, edges: np.ndarray,
               *, same: bool, box_size: float | None, chunk: int = 128) -> np.ndarray:
    """Count unordered auto pairs and ordered cross pairs, without self-pairs."""
    hist = np.zeros(len(edges) - 1, dtype=np.int64)
    for start in range(0, len(a), chunk):
        delta = a[start:start + chunk, None, :] - b[None, :, :]
        if box_size is not None:
            delta -= box_size * np.rint(delta / box_size)
        distances = np.sqrt(np.sum(delta * delta, axis=-1))
        if same:
            ii, jj = np.indices(distances.shape)
            valid = jj > ii + start
            distances = distances[valid]
        hist += np.histogram(distances, bins=edges)[0]
    return hist


def landy_szalay(data: np.ndarray, random: np.ndarray, edges: np.ndarray,
                  *, box_size: float | None = None) -> PairResult:
    """Landy–Szalay xi=(DD-2DR+RR)/RR with finite-catalogue normalization.

    Bins with RR=0 return NaN, never a misleading zero. This estimator uses
    unweighted positions; correlated uncertainties require separate treatment.
    """
    d, r, e = _points(data), _points(random), _edges(edges)
    if box_size is not None:
        if not np.isfinite(box_size) or box_size <= 0:
            raise ValueError("box_size must be positive and finite")
        if np.any(d < 0) or np.any(r < 0) or np.any(d >= box_size) or np.any(r >= box_size):
            raise ValueError("periodic coordinates must be in [0, box_size)")
    dd = _histogram(d, d, e, same=True, box_size=box_size)
    dr = _histogram(d, r, e, same=False, box_size=box_size)
    rr = _histogram(r, r, e, same=True, box_size=box_size)
    dd_n = dd / (len(d) * (len(d) - 1) / 2)
    dr_n = dr / (len(d) * len(r))
    rr_n = rr / (len(r) * (len(r) - 1) / 2)
    xi = np.full(len(rr), np.nan)
    np.divide(dd_n - 2 * dr_n + rr_n, rr_n, out=xi, where=rr_n > 0)
    return PairResult(e, dd, dr, rr, xi)


def spatial_jackknife(data: np.ndarray, random: np.ndarray, edges: np.ndarray,
                      *, box_size: float, divisions: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """Delete-one-3D-cell jackknife for xi and its covariance.

    Removes the same spatial cell from both catalogues. Assumes a cubic footprint;
    for actual surveys use a physically motivated angular/redshift jackknife.
    """
    d, r = _points(data), _points(random)
    if not np.isfinite(box_size) or box_size <= 0 or divisions < 2:
        raise ValueError("box_size > 0 and divisions >= 2 required")
    if np.any(d < 0) or np.any(r < 0) or np.any(d >= box_size) or np.any(r >= box_size):
        raise ValueError("coordinates must be in the cube")

    def labels(points: np.ndarray) -> np.ndarray:
        cell = np.floor(points / box_size * divisions).astype(int)
        return cell[:, 0] + divisions * (cell[:, 1] + divisions * cell[:, 2])

    dl, rl = labels(d), labels(r)
    replicates = []
    for cell in range(divisions ** 3):
        ds, rs = d[dl != cell], r[rl != cell]
        if len(ds) < 2 or len(rs) < 2:
            raise ValueError(f"insufficient points after deleting cell {cell}")
        result = landy_szalay(ds, rs, edges)
        if not np.isfinite(result.xi).all():
            raise ValueError("empty RR bin in jackknife; widen bins or increase random catalogue")
        replicates.append(result.xi)
    samples = np.stack(replicates)
    mean = samples.mean(axis=0)
    diff = samples - mean
    covariance = (len(samples) - 1) / len(samples) * (diff.T @ diff)
    return samples, covariance


def random_catalogue_sensitivity(data: np.ndarray, edges: np.ndarray,
                                 *, box_size: float, random_sizes: tuple[int, ...],
                                 repeats: int = 8, seed: int = 20261019) -> dict[int, dict[str, np.ndarray]]:
    """Quantify variation from independent uniform random catalogues at fixed data.

    Does not include cosmic variance, data Poisson uncertainty, or survey masks.
    """
    _points(data)
    _edges(edges)
    if box_size <= 0 or repeats < 2 or not random_sizes or any(n < 2 for n in random_sizes):
        raise ValueError("invalid box_size, repeats, or random_sizes")
    rng = np.random.default_rng(seed)
    output = {}
    for n in random_sizes:
        xi = []
        for _ in range(repeats):
            random = rng.uniform(0, box_size, size=(n, 3))
            xi.append(landy_szalay(data, random, edges).xi)
        stack = np.asarray(xi)
        output[n] = {"mean": np.nanmean(stack, axis=0),
                     "std": np.nanstd(stack, axis=0, ddof=1),
                     "finite_fraction": np.isfinite(stack).mean(axis=0)}
    return output
