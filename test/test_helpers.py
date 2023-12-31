"""
Tests for helper functions
"""
from astropy import units as u
import pytest
import numpy as np

import vspec_vsm.helpers as helpers

def test_round_teff():
    """
    Test `VSPEC.helpers.round_teff`
    """
    teff = 100.3*u.K
    assert helpers.round_teff(teff) == 100*u.K
    
def test_get_angle_between():
    cases = [
        {'args': [0, 0, 0, 0]*u.deg, 'pred':0*u.deg},
        {'args': [0, 0, 90, 0]*u.deg, 'pred':90*u.deg},
        {'args': [0, 0, 40, 0]*u.deg, 'pred':40*u.deg},
        {'args': [0, 0, 0, 120]*u.deg, 'pred':120*u.deg},
        {'args': [0, 0, 0, -120]*u.deg, 'pred':120*u.deg},
        {'args': [90, 0, -90, 0]*u.deg, 'pred':180*u.deg},
    ]
    for case in cases:
        assert helpers.get_angle_between(
            *case['args']).to_value(u.deg) == pytest.approx(case['pred'].to_value(u.deg), abs=1)

@pytest.mark.filterwarnings('error')
def test_proj_ortho():
    lat0 = 30 * u.deg
    lon0 = 45 * u.deg
    lats = np.array([40, 50, 60]) * u.deg
    lons = np.array([60, 70, 80]) * u.deg

    x, y = helpers.proj_ortho(lat0, lon0, lats, lons)

    # Check the length of output arrays
    assert len(x) == len(y) == len(lats)

    # Check individual projected coordinates
    assert x[0] == pytest.approx(0.19825408, rel=0.01)
    assert x[1] == pytest.approx(0.27164737, rel=0.01)
    assert x[2] == pytest.approx(0.28682206, rel=0.01)

    assert y[0] == pytest.approx(0.18671294, rel=0.01)
    assert y[1] == pytest.approx(0.37213692, rel=0.01)
    assert y[2] == pytest.approx(0.54519419, rel=0.01)

    # Check for incorrect input type
    with pytest.raises(TypeError):
        helpers.proj_ortho(lat0.value, lon0.value, lats, lons)
        
    lat0 = 0 * u.deg
    lon0 = 10 * u.deg
    lats = np.array([0, 0, 0]) * u.deg
    lons = np.array([100, 40, 10]) * u.deg
    x, y = helpers.proj_ortho(lat0, lon0, lats, lons)
    assert y[0] == pytest.approx(0., rel=0.01)
    assert y[1] == pytest.approx(0., rel=0.01)
    assert y[2] == pytest.approx(0., rel=0.01)
    
    assert x[0] == pytest.approx(1., rel=0.01)
    assert x[1] == pytest.approx(0.5, rel=0.01)
    assert x[2] == pytest.approx(0., rel=0.01)

def test_circle_intersection():
    cases = [
        {'args': [0, 0, 1], 'pred':1.0},
        {'args': [0, 0, 0.5], 'pred':1.0},
        {'args': [0, 0.5, 0.5], 'pred':1.0},
        {'args': [2, 0, 0.5], 'pred':0.0},
    ]
    for case in cases:
        calc = helpers.calc_circ_fraction_inside_unit_circle(*case['args'])
        pred = case['pred']
        assert calc == pytest.approx(pred, rel=1e-6)
    assert helpers.calc_circ_fraction_inside_unit_circle(1, 0, 1) < 0.5
    assert helpers.calc_circ_fraction_inside_unit_circle(0.6, 0, 0.5) > 0.5

def test_clip_teff():
    """
    Test `VSPEC.helpers.clip_teff`
    """
    low_bound = 2300 * u.K
    high_bound = 3900 * u.K

    teff = 2000 * u.K
    clipped_teff = helpers.clip_teff(teff)
    assert clipped_teff == low_bound

    teff = 6500 * u.K
    clipped_teff = helpers.clip_teff(teff)
    assert clipped_teff == high_bound

    teff = 3000 * u.K
    clipped_teff = helpers.clip_teff(teff)
    assert clipped_teff == teff
