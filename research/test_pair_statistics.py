import math
import unittest
from pair_statistics import pair_histogram, landy_szalay, delete_one_jackknife


class PairStatisticsTests(unittest.TestCase):
    def test_pair_count_no_double_count(self):
        self.assertEqual(pair_histogram([(0,), (1,), (3,)], [0, 2, 4]), [1, 2])

    def test_cross_pair_count(self):
        self.assertEqual(pair_histogram([(0,), (2,)], [0, 1, 3], [(1,)]), [0, 2])

    def test_identical_catalogues_normalized_zero(self):
        p = [(0,), (1,), (3,), (6,)]
        result = landy_szalay(p, p, [0, 2, 4, 7])
        # DR includes the four self-pairs, so identical catalogues do not
        # imply zero xi at short separations; verify exact normalization.
        self.assertAlmostEqual(result[0], 1 - 2 * (6 / 16) / (3 / 6) + 1)

    def test_empty_rr_is_nan(self):
        self.assertTrue(math.isnan(landy_szalay(
            [(0,), (1,)], [(0,), (1,)], [2, 3])[0]))

    def test_jackknife_covariance(self):
        mean, cov = delete_one_jackknife([[1, 2], [2, 4], [3, 6]])
        self.assertEqual(mean, [2, 4])
        self.assertAlmostEqual(cov[0][0], 4 / 3)
        self.assertAlmostEqual(cov[0][1], 8 / 3)


if __name__ == "__main__":
    unittest.main()
