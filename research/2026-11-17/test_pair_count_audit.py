import numpy as np
import pytest
from pair_count_audit import landy_szalay, leave_region_out, jackknife_covariance, covariance_spectrum

def test_zero_correlation():
    assert landy_szalay([1, 2], [1, 2], [1, 2]).tolist() == pytest.approx([0, 0])

def test_pair_count_example():
    assert landy_szalay(1.2, 1.0, 1.0) == pytest.approx(.2)

def test_region_replicates():
    x = leave_region_out([[1.2], [1.0], [.8]], [[1.], [1.], [1.]], [[1.], [1.], [1.]])
    assert x[:, 0] == pytest.approx([.2, 0, -.2])

def test_jackknife_covariance():
    x = np.array([[1., 2.], [2., 4.], [3., 6.]])
    assert np.allclose(jackknife_covariance(x), [[4/3, 8/3], [8/3, 16/3]])

def test_rank_deficiency():
    x = np.array([[1., 11.], [2., 22.], [2., 22.]])
    spectrum = covariance_spectrum(x)
    assert spectrum["rank"] == 1
    assert np.isinf(spectrum["condition"])

def test_full_rank_condition():
    x = np.array([[-2., 0.], [2., 0.], [0., -1.], [0., 1.]])
    spectrum = covariance_spectrum(x)
    assert spectrum["rank"] == 2
    assert spectrum["condition"] == pytest.approx(4)

def test_invalid_counts():
    with pytest.raises(ValueError):
        landy_szalay(1, 1, 0)
    with pytest.raises(ValueError):
        jackknife_covariance([[1], [2]])
