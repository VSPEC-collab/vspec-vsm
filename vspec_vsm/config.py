"""
Configurations
"""

from astropy import units as u
import numpy as np

MSH = u.def_unit('msh', 1e-6 * 0.5 * 4*np.pi*u.R_sun**2)
"""
Micro-solar hemisphere

This is a standard unit in heliophysics that
equals one millionth of one half the surface area of the Sun.

:type: astropy.units.Unit
"""

stellar_area_unit = MSH
"""
The standard stellar surface area unit.

This unit is used to represent the surface area of stars in VSPEC.
The micro-solar hemisphere is chosen because most Sun Spot literature uses
this unit.

:type: astropy.units.Unit
"""

starspot_initial_area = 10*MSH
"""
Initial ``StarSpot`` area.

Because spots grow exponentially, they can't start at 0 area.
When they are born they are given this small area.

:type: astropy.units.Quantity

.. todo::
    This should optionaly be set by the user. So that smaller
    star spot area regimes are accessible.
"""

NLAT = 500
"""
The default latitude resolution for the stellar model. This should
be set by finding a balance between noticing small changes in spots/faculae
and computation time.

:type: int
"""

NLON = 1000
"""
The default longitude resolution for the stellar model. This should
be set by finding a balance between noticing small changes in spots/faculae
and computation time.

:type: int
"""

grid_teff_bounds = (2300*u.K, 3900*u.K)
"""
The limits on the effective temperature allowed by the grid.

:type: tuple of astropy.units.Quantity
"""
