"""Run with python -m unittest discover -s tests -v."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "science"))
from clustering_diagnostics import pair_histogram, landy_szalay, region_jackknife, fit_bao_amplitude


class ClusteringDiagnosticsTest(unittest.TestCase):
    def test_unordered_pairs(self):
        a = np.array([[0, 0, 0], [1, 0, 0], [3, 0, 0]], float)
        np.testing.assert_array_equal(pair_histogram(a, a, [0, 1.5, 2.5, 3.5], True), [1, 1, 1])

    def test_cross_pairs(self):
        a = np.array([[0, 0, 0], [1, 0, 0]], float)
        b = np.array([[2, 0, 0], [4, 0, 0]], float)
        self.assertEqual(pair_histogram(a, b, [0, 3, 5]).sum(), 4)

    def test_landy_szalay_independent_uniform_catalogs(self):
        rng = np.random.default_rng(23)
        data = rng.uniform(0, 10, (450, 3))
        randoms = rng.uniform(0, 10, (650, 3))
        result = landy_szalay(data, randoms, [0, 3, 5, 7, 10, 18])
        self.assertTrue(np.all(np.isfinite(result["xi"])))
        self.assertLess(np.max(np.abs(result["xi"])), .12)

    def test_landy_szalay_empty_rr_is_nan(self):
        a = np.array([[0, 0, 0], [1, 0, 0]], float)
        r = landy_szalay(a, a, [0, .1, 2])
        self.assertTrue(np.isnan(r["xi"][0]))

    def test_jackknife_covariance_psd(self):
        rng = np.random.default_rng(32)
        d = rng.uniform(0, 10, (180, 3))
        r = rng.uniform(0, 10, (400, 3))
        dr = np.repeat(np.arange(6), 30)
        rr = np.arange(400) % 6
        out = region_jackknife(d, r, dr, rr, [0, 4, 7, 11, 18])
        self.assertEqual(out["replicates"].shape, (6, 4))
        self.assertGreaterEqual(np.linalg.eigvalsh(out["covariance"]).min(), -1e-12)

    def test_jackknife_invalid_regions(self):
        a = np.random.default_rng(4).uniform(size=(10, 3))
        with self.assertRaises(ValueError):
            region_jackknife(a, a, np.zeros(10), np.zeros(10), [0, 2])

    def test_bao_template_recovery(self):
        x = np.linspace(70, 140, 45)
        baseline = .005 * np.exp(-(x - 70) / 30)
        template = np.exp(-.5 * ((x - 105) / 10)**2)
        cov = np.diag(np.full(len(x), .001**2))
        out = fit_bao_amplitude(x, baseline + .03 * template, cov, baseline)
        self.assertAlmostEqual(out.amplitude, .03, places=12)
        self.assertAlmostEqual(out.chi2, 0, places=9)
        self.assertGreater(out.delta_chi2_against_zero, 0)

    def test_bao_correlated_noise_uncertainty(self):
        x = np.linspace(80, 130, 15)
        cov = np.eye(len(x)) * .0001
        cov += .00003
        out = fit_bao_amplitude(x, np.zeros(len(x)), cov, np.zeros(len(x)))
        naive = fit_bao_amplitude(x, np.zeros(len(x)), np.diag(np.diag(cov)), np.zeros(len(x)))
        self.assertGreater(out.sigma_amplitude, naive.sigma_amplitude)

    def test_bao_singular_covariance(self):
        x = np.linspace(80, 130, 5)
        with self.assertRaises(ValueError):
            fit_bao_amplitude(x, np.ones(5), np.zeros((5, 5)), np.zeros(5))


if __name__ == "__main__":
    unittest.main()
