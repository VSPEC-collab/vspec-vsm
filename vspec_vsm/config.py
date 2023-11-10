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