"""
Tests for helper functions
"""
from astropy import units as u

import vspec_vsm.helpers as helpers

def test_round_teff():
    """
    Test `VSPEC.helpers.round_teff`
    """
    teff = 100.3*u.K
    assert helpers.round_teff(teff) == 100*u.K
