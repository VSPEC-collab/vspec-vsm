"""
Tests for the vspec_vsm.flare module.
"""
import pytest
from astropy import units as u
import numpy as np

from vspec_vsm.flares import StellarFlare, FlareGenerator, FlareCollection


def get_generator():
    return FlareGenerator(
        dist_teff_mean=1000*u.K,
        dist_teff_sigma=100*u.K,
        dist_fwhm_mean=3*u.hr,
        dist_fwhm_logsigma=0.2,
        alpha=-0.829,
        beta=26.87,
        min_energy=1e33*u.erg,
        cluster_size=2
    )

@pytest.fixture
def generator():
    return get_generator()

@pytest.fixture
def collection():
    gen = get_generator()
    flares = gen.generate_flare_series(10*u.day)
    return FlareCollection(flares)

def test_generate_flare_series(generator: FlareGenerator):
    time = 10 * u.day
    flares = generator.generate_flare_series(time)
    assert isinstance(flares, list)
    assert all(isinstance(flare, StellarFlare) for flare in flares)

def test_get_visible_flares_in_timeperiod(collection: FlareCollection):
    tstart = 0 * u.day
    tfinish = 10 * u.day
    visible_flares = collection.get_visible_flares_in_timeperiod(tstart, tfinish)
    assert isinstance(visible_flares, list)
    assert all(isinstance(flare, StellarFlare) for flare in visible_flares)

def test_gen_fwhm(generator: FlareGenerator):
    n_flares = 10
    fwhm = generator.gen_fwhm(n_flares)
    assert isinstance(fwhm, u.Quantity)
    assert len(fwhm) == n_flares

def test_get_flares_in_timeperiod(collection: FlareCollection):
    tstart = 0 * u.day
    tfinish = 10 * u.day
    flares = collection.get_flares_in_timeperiod(tstart, tfinish)
    assert isinstance(flares, list)
    assert all(isinstance(flare, StellarFlare) for flare in flares)

def test_get_peaks(generator: FlareGenerator):
    n_flares = 10
    tpeaks = generator.get_peaks(n_flares)
    assert isinstance(tpeaks, np.ndarray), f'Got {type(tpeaks)}'
    assert len(tpeaks) == n_flares, f'Got {len(tpeaks)} peaks, expected {n_flares}'

def test_gen_teffs(generator: FlareGenerator):
    n_flares = 10
    teffs = generator.gen_teffs(n_flares)
    assert isinstance(teffs, u.Quantity)
    assert len(teffs) == n_flares



