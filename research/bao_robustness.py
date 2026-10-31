"""Synthetic BAO peak, correlated-noise and ruler-identifiability diagnostics.

No galaxy survey catalogue or calibrated distance measurement is used.
"""
from __future__ import annotations

import numpy as np


def broadband(r: np.ndarray, r0: float = 5, gamma: float = 1.8,
              bias: float = 1.5) -> np.ndarray:
    r = np.asarray(r, dtype=float)
    if np.any(~np.isfinite(r)) or np.any(r <= 0) or r0 <= 0 or gamma <= 0 or bias <= 0:
        raise ValueError("r, r0, gamma and bias must be positive finite")
    return bias**2 * (r / r0) ** (-gamma) * np.exp(-r / 65)


def xi_model(r: np.ndarray, amplitude: float = 0.04, center: float = 105,
             width: float = 10) -> np.ndarray:
    if not np.isfinite(amplitude) or amplitude < 0 or width <= 0 or center <= 0:
        raise ValueError("invalid BAO parameters")
    r = np.asarray(r, dtype=float)
    return broadband(r) + amplitude * np.exp(-0.5 * ((r - center) / width)**2)


def parabolic_peak(r: np.ndarray, y: np.ndarray, lo: float = 80,
                   hi: float = 130) -> float:
    """Sub-grid quadratic peak; only valid for a well-sampled interior maximum."""
    r, y = np.asarray(r, float), np.asarray(y, float)
    if r.ndim != 1 or r.size != y.size or r.size < 5 or not np.all(np.diff(r) > 0):
        raise ValueError("r must be a sorted grid matching y")
    mask = (r >= lo) & (r <= hi)
    inds = np.flatnonzero(mask)
    if inds.size < 3 or np.any(~np.isfinite(y)):
        raise ValueError("invalid peak window")
    k = inds[np.argmax(y[mask])]
    if k <= inds[0] or k >= inds[-1]:
        raise ValueError("peak is at the search boundary")
    coeff = np.polyfit(r[k-1:k+2] - r[k], y[k-1:k+2], 2)
    if coeff[0] >= 0:
        raise ValueError("no local concave maximum")
    return float(r[k] - coeff[1] / (2 * coeff[0]))


def ar1_whiten(v: np.ndarray, sigma: float, rho: float) -> np.ndarray:
    """Exact innovations whitening for stationary equally spaced AR(1) noise."""
    v = np.asarray(v, float)
    if v.ndim != 1 or v.size < 2 or sigma <= 0 or not np.isfinite(sigma) or abs(rho) >= 1:
        raise ValueError("invalid AR(1) inputs")
    w = np.empty_like(v)
    w[0] = v[0] / sigma
    w[1:] = (v[1:] - rho * v[:-1]) / (sigma * np.sqrt(1 - rho**2))
    return w


def profile_peak(r: np.ndarray, y: np.ndarray, centers: np.ndarray,
                 sigma: float, rho: float = 0, width: float = 10) -> dict[str, float]:
    """Profile Gaussian BAO amplitude on a fixed known broadband, with AR(1) GLS."""
    r, y, centers = np.asarray(r, float), np.asarray(y, float), np.asarray(centers, float)
    if r.ndim != 1 or y.shape != r.shape or centers.ndim != 1 or centers.size < 3:
        raise ValueError("invalid grids")
    if not np.all(np.diff(r) > 0) or not np.all(np.diff(centers) > 0) or width <= 0:
        raise ValueError("grids must be increasing and width positive")
    if not np.all(np.isfinite(y)):
        raise ValueError("data must be finite")
    z = ar1_whiten(y - broadband(r), sigma, rho)
    chi2 = []
    amplitudes = []
    for center in centers:
        g = ar1_whiten(np.exp(-0.5 * ((r - center) / width)**2), sigma, rho)
        a = float(np.dot(z, g) / np.dot(g, g))
        amplitudes.append(a)
        chi2.append(float(np.dot(z - a * g, z - a * g)))
    chi2 = np.asarray(chi2)
    best = int(np.argmin(chi2))
    valid = centers[chi2 <= chi2[best] + 1]
    return {"center_hmpc": float(centers[best]), "amplitude": amplitudes[best],
            "chi2": float(chi2[best]), "profile_68_width_hmpc": float((valid[-1] - valid[0]) / 2)}


def conditional_h0(peak_hmpc: float, sound_horizon_mpc: float,
                   sigma_peak: float = 0, sigma_ruler: float = 0,
                   correlation: float = 0) -> dict[str, float]:
    """Unit conversion H0=100*r_peak[h^-1 Mpc]/r_d[Mpc], conditional on r_d.

    Not a BAO distance likelihood or an independent determination of H0.
    """
    if any(not np.isfinite(x) for x in (peak_hmpc, sound_horizon_mpc,
                                         sigma_peak, sigma_ruler, correlation)):
        raise ValueError("inputs must be finite")
    if peak_hmpc <= 0 or sound_horizon_mpc <= 0 or min(sigma_peak, sigma_ruler) < 0 or abs(correlation) > 1:
        raise ValueError("invalid ruler or uncertainty")
    h0 = 100 * peak_hmpc / sound_horizon_mpc
    a = 100 / sound_horizon_mpc
    b = -100 * peak_hmpc / sound_horizon_mpc**2
    variance = a*a*sigma_peak**2 + b*b*sigma_ruler**2 + 2*a*b*correlation*sigma_peak*sigma_ruler
    return {"conditional_h0_km_s_mpc": float(h0),
            "linearized_sigma_h0": float(np.sqrt(max(0, variance)))}


def diagnostic_study() -> dict[str, object]:
    r = np.linspace(70, 140, 281)
    y = xi_model(r)
    peak_raw = parabolic_peak(r, y)
    peak_subtracted = parabolic_peak(r, y - broadband(r))
    centers = np.linspace(95, 115, 201)
    # Explicitly noise-free synthetic fit: sigma is a hypothetical weight scale.
    fits = {str(rho): profile_peak(r, y, centers, sigma=0.003, rho=rho)
            for rho in (0, 0.5, 0.8)}
    return {"raw_peak_hmpc": peak_raw, "broadband_subtracted_peak_hmpc": peak_subtracted,
            "conditional_gls_noise_free": fits,
            "ruler_assumptions": {
                "rd_147.09": conditional_h0(105, 147.09, sigma_peak=2, sigma_ruler=0.3),
                "rd_150": conditional_h0(105, 150, sigma_peak=2, sigma_ruler=0.3),
            }}


if __name__ == "__main__":
    import json
    print(json.dumps(diagnostic_study(), indent=2))