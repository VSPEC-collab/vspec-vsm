"""
Plot a map of the stellar surface
=================================

This example initializes a ``Star`` object and plots it.
"""
from astropy import units as u
import numpy as np

from vspec_vsm.spots import SpotCollection, StarSpot
from vspec_vsm.faculae import FaculaCollection, Facula
from vspec_vsm.star import Star
from vspec_vsm.config import MSH

SEED = 10
rng = np.random.default_rng(SEED)

# %%
# Initialize the star
# -------------------
#
# First, let's initialize a ``Star`` object.
#
# It needs to be populated by spots and faculae.

n_spots = 10
n_faculae = 10
spot_area = 1000*MSH
facula_radius = 10000*u.km
facula_depth = 10000*u.km

spots = SpotCollection(
    *[
        StarSpot(
            lat=(rng.random() - 0.5)*120*u.deg,
            lon=rng.random()*360*u.deg,
            area_max=spot_area,
            area_current=spot_area,
            teff_umbra=2700*u.K,
            teff_penumbra=2900*u.K,
            area_over_umbra_area=5.,
            is_growing=False,
            growth_rate=0./u.day,
            decay_rate=0*MSH/u.day
        ) for _ in range(n_spots)
    ]
)

faculae = FaculaCollection(
    *[
        Facula(
            lat=(rng.random() - 0.5)*120*u.deg,
            lon=rng.random()*360*u.deg,
            r_max=facula_radius,
            r_init=facula_radius,
            depth=facula_radius,
            lifetime=5*u.hr,
            floor_teff_slope=0*u.K/u.km,
            floor_teff_min_rad=0.1*facula_radius,
            floor_teff_base_dteff=-500*u.K,
            wall_teff_slope=0*u.K/u.km,
            wall_teff_intercept=300*u.K,
            growing=False,
        ) for _ in range(n_faculae)
    ]
)

star_teff = 3300*u.K
star_radius = 0.15*u.R_sun
star_period = 40*u.day
Nlat = 500
Nlon = 1000
ld_params = dict(u1=0.3, u2=0.1)

star = Star(
    teff=star_teff,
    radius=star_radius,
    period=star_period,
    spots=spots,
    faculae=faculae,
    grid_params=(Nlat, Nlon),
    **ld_params
)
# %%
# Simulate the disk
# -----------------
#
# Now let's decide a viewing angle and get an image of the surface.

lon0 = 0*u.deg
lat0 = 0*u.deg

star.plot_surface(lat0, lon0)

# %%
# Add a transit
# -------------
#
# Let's throw in a transiting planet just for fun.

pl_radius = 1*u.R_earth
pl_orbit = 0.05*u.AU
inclination = 89.8*u.deg
phase = 180.4*u.deg

star.plot_surface(lat0, lon0, None, pl_orbit, pl_radius, phase, inclination)
