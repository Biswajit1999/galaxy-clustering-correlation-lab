import unittest

import numpy as np

from research.bao_robustness import (ar1_whiten, broadband, conditional_h0,
                                     diagnostic_study, parabolic_peak, profile_peak, xi_model)


class TestBAORobustness(unittest.TestCase):
    def test_exact_worker_model(self):
        r = np.array([5., 105.])
        expected = 1.5**2 * (r / 5)**(-1.8) * np.exp(-r/65) + 0.04*np.exp(-0.5*((r-105)/10)**2)
        np.testing.assert_allclose(xi_model(r), expected, rtol=1e-14)

    def test_broadband_subtraction_corrects_peak(self):
        r = np.linspace(70, 140, 281)
        raw = parabolic_peak(r, xi_model(r))
        clean = parabolic_peak(r, xi_model(r) - broadband(r))
        self.assertLess(raw, 105)
        self.assertAlmostEqual(clean, 105, places=8)

    def test_peak_grid_convergence(self):
        vals = []
        for n in (141, 281, 561):
            r = np.linspace(70, 140, n)
            vals.append(parabolic_peak(r, xi_model(r) - broadband(r)))
        self.assertLess(max(vals)-min(vals), 1e-7)

    def test_gls_exact_center_noiseless(self):
        r = np.linspace(70, 140, 141)
        fit = profile_peak(r, xi_model(r), np.linspace(95,115,201), sigma=.003, rho=.6)
        self.assertAlmostEqual(fit["center_hmpc"], 105)
        self.assertAlmostEqual(fit["amplitude"], .04, places=9)

    def test_ar1_whitening_covariance(self):
        rho, sigma, n = .6, .01, 6
        cov = sigma**2 * rho**np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
        whiten = np.column_stack([ar1_whiten(np.eye(n)[:,i],sigma,rho) for i in range(n)])
        np.testing.assert_allclose(whiten @ cov @ whiten.T, np.eye(n), atol=1e-12)

    def test_correlated_noise_inflates_center_uncertainty(self):
        r = np.linspace(70, 140, 141)
        centers = np.linspace(95, 115, 401)
        a = profile_peak(r, xi_model(r), centers, sigma=.003, rho=0)
        b = profile_peak(r, xi_model(r), centers, sigma=.003, rho=.7)
        self.assertGreater(b["profile_68_width_hmpc"], a["profile_68_width_hmpc"])

    def test_conditional_ruler_and_uncertainty(self):
        a = conditional_h0(105,147.09,sigma_peak=2,sigma_ruler=.3)
        self.assertAlmostEqual(a["conditional_h0_km_s_mpc"], 100*105/147.09)
        self.assertGreater(a["linearized_sigma_h0"], 0)
        self.assertNotAlmostEqual(a["conditional_h0_km_s_mpc"], conditional_h0(105,150)["conditional_h0_km_s_mpc"])

    def test_bad_inputs(self):
        with self.assertRaises(ValueError):
            conditional_h0(105,0)
        with self.assertRaises(ValueError):
            ar1_whiten(np.ones(5),.01,1)
        with self.assertRaises(ValueError):
            parabolic_peak(np.array([1,2,3]),np.array([0,1,0]))


if __name__ == "__main__":
    unittest.main()