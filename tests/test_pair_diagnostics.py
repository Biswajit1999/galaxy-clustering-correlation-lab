import numpy as np
import pytest
from research.pair_diagnostics import landy_szalay, spatial_jackknife, random_catalogue_sensitivity


def catalogues(seed=2):
    rng = np.random.default_rng(seed)
    return rng.uniform(0, 100, (160, 3)), rng.uniform(0, 100, (500, 3))


def test_exact_pair_normalization_one_full_bin():
    d, r = catalogues()
    res = landy_szalay(d, r, [0, 200])
    assert res.dd[0] == 160 * 159 // 2
    assert res.dr[0] == 160 * 500
    assert res.rr[0] == 500 * 499 // 2
    assert res.xi[0] == pytest.approx(0, abs=1e-14)


def test_empty_random_bin_explicitly_nan():
    d = np.array([[0,0,0],[1,0,0]], float)
    r = np.array([[0,0,0],[100,0,0]], float)
    assert np.isnan(landy_szalay(d, r, [0,2]).xi[0])


def test_permutation_invariance():
    d,r = catalogues()
    a=landy_szalay(d,r,[0,15,30,45,60,90])
    b=landy_szalay(d[::-1],r[::-1],[0,15,30,45,60,90])
    np.testing.assert_array_equal(a.dd,b.dd)
    np.testing.assert_array_equal(a.dr,b.dr)
    np.testing.assert_array_equal(a.rr,b.rr)


def test_periodic_minimum_image_pair():
    d=np.array([[0.1,1,1],[9.9,1,1]])
    r=np.array([[0.1,2,2],[9.9,2,2]])
    assert landy_szalay(d,r,[0,0.5,3],box_size=10).dd[0]==1
    assert landy_szalay(d,r,[0,0.5,3]).dd[0]==0


def test_jackknife_covariance_symmetric_psd():
    d,r=catalogues()
    samples,cov=spatial_jackknife(d,r,[0,30,60,90],box_size=100)
    assert samples.shape==(8,3)
    np.testing.assert_allclose(cov,cov.T,atol=1e-14)
    assert np.linalg.eigvalsh(cov).min()>=-1e-12


def test_jackknife_zero_for_complete_pair_bin():
    d,r=catalogues()
    _,cov=spatial_jackknife(d,r,[0,200],box_size=100)
    assert abs(cov[0,0])<1e-25


def test_random_sensitivity_reproducible():
    d,_=catalogues()
    a=random_catalogue_sensitivity(d,[0,30,60,90],box_size=100,random_sizes=(250,),repeats=4,seed=6)
    b=random_catalogue_sensitivity(d,[0,30,60,90],box_size=100,random_sizes=(250,),repeats=4,seed=6)
    np.testing.assert_array_equal(a[250]['mean'],b[250]['mean'])
    assert np.all(a[250]['finite_fraction']==1)


@pytest.mark.parametrize('bad',[[0,0,1],[-1,2],[0,np.inf]])
def test_invalid_bins(bad):
    d,r=catalogues()
    with pytest.raises(ValueError): landy_szalay(d,r,bad)


def test_invalid_coordinates():
    d,r=catalogues()
    d[0,0]=np.nan
    with pytest.raises(ValueError): landy_szalay(d,r,[0,30])


def test_invalid_periodic_box():
    d,r=catalogues()
    with pytest.raises(ValueError): landy_szalay(d,r,[0,30],box_size=50)
