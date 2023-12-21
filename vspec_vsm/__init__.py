"""
Most of the interface is defined in by the ``Star`` class.
"""

__version__ = "0.3.2"

from vspec_vsm.star import Star
from vspec_vsm.spots import SpotCollection, StarSpot, SpotGenerator
from vspec_vsm.faculae import FaculaCollection, Facula, FaculaGenerator
from vspec_vsm.flares import StellarFlare, FlareCollection, FlareGenerator
from vspec_vsm.granules import Granulation, GranulationKernel
from vspec_vsm.config import MSH
from vspec_vsm.coordinate_grid import CoordinateGrid
