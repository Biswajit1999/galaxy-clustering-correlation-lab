"""Small-catalogue Euclidean Landy-Szalay and spatial delete-one jackknife.

Distances are in a single consistent Cartesian unit. No survey selection function
or cosmological redshift-to-distance conversion is supplied by this module.
"""
import math


def pair_histogram(points_a, edges, points_b=None):
    """Count unordered within-sample pairs, or all cross-sample pairs."""
    if len(edges) < 2 or any(not math.isfinite(v) for v in edges) or any(
        a >= b for a, b in zip(edges, edges[1:])
    ):
        raise ValueError("edges must be finite and strictly increasing")
    if not points_a or (points_b is not None and not points_b):
        raise ValueError("samples must not be empty")
    dimensions = len(points_a[0])
    samples = list(points_a) + (list(points_b) if points_b is not None else [])
    if dimensions == 0 or any(len(p) != dimensions or
                              any(not math.isfinite(x) for x in p) for p in samples):
        raise ValueError("points require matching finite coordinates")
    counts = [0] * (len(edges) - 1)
    other = points_a if points_b is None else points_b
    for i, p in enumerate(points_a):
        for j, q in enumerate(other):
            if points_b is None and j <= i:
                continue
            r = math.dist(p, q)
            for k, (lo, hi) in enumerate(zip(edges, edges[1:])):
                if lo <= r < hi:
                    counts[k] += 1
                    break
    return counts


def landy_szalay(data, randoms, edges):
    """Return xi=(DD-2DR+RR)/RR with pair-count normalization.

    A zero RR bin is undefined and is returned as NaN, not as zero.
    """
    nd, nr = len(data), len(randoms)
    if nd < 2 or nr < 2:
        raise ValueError("at least two data and random points required")
    dd = pair_histogram(data, edges)
    dr = pair_histogram(data, edges, randoms)
    rr = pair_histogram(randoms, edges)
    dd_norm = nd * (nd - 1) / 2
    rr_norm = nr * (nr - 1) / 2
    dr_norm = nd * nr
    return [
        ((d / dd_norm) - 2 * (x / dr_norm) + (r / rr_norm)) / (r / rr_norm)
        if r else math.nan
        for d, x, r in zip(dd, dr, rr)
    ]


def delete_one_jackknife(estimates):
    """Mean and jackknife covariance from equal-weight region-deletion vectors.

    Region membership and re-estimation must be performed upstream; each row is
    the statistic with one spatial region removed. All values must be finite.
    """
    n = len(estimates)
    if n < 3 or not estimates[0]:
        raise ValueError("at least three nonempty leave-region-out vectors")
    m = len(estimates[0])
    if any(len(row) != m or any(not math.isfinite(x) for x in row)
           for row in estimates):
        raise ValueError("ragged or nonfinite estimates")
    mean = [sum(row[j] for row in estimates) / n for j in range(m)]
    covariance = [
        [(n - 1) / n * sum((row[i] - mean[i]) * (row[j] - mean[j])
                           for row in estimates) for j in range(m)]
        for i in range(m)
    ]
    return mean, covariance
