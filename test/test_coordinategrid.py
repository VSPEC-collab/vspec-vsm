"""
Tests for the CoordinateGrid class.
"""

import pytest
import numpy as np
from vspec_vsm.coordinate_grid import RectangularGrid


def test_coordinategrid():
    """
    Tests for the CoordinateGrid class.
    """
    grid = RectangularGrid(Nlat=100, Nlon=200)
    assert grid.Nlat == 100
    assert grid.Nlon == 200


def test_coodinategrid_typeerror():
    """
    Tests for the CoordinateGrid class initialized with floats.
    """
    with pytest.raises(TypeError):
        _ = RectangularGrid(Nlat=100, Nlon=200.0)
    with pytest.raises(TypeError):
        _ = RectangularGrid(Nlat=100.0, Nlon=200)


def test_coordinategrid_oned():
    """
    Tests for the CoordinateGrid oned method.
    """
    grid = RectangularGrid(Nlat=100, Nlon=200)
    lats, lons = grid.oned()
    assert lats.shape == (100,)
    assert lons.shape == (200,)


def test_coordinategrid_grid():
    """
    Tests for the CoordinateGrid grid method.
    """
    grid = RectangularGrid(Nlat=100, Nlon=200)
    lats, lons = grid.grid()
    assert lats.shape == (200, 100)
    assert lons.shape == (200, 100)


def test_coordinategrid_zeros():
    """
    Tests for the CoordinateGrid zeros method.
    """
    grid = RectangularGrid(Nlat=100, Nlon=200)
    arr = grid.zeros()
    assert arr.shape == (200, 100)
    assert arr.dtype == np.float32
    assert np.all(arr == 0)


def test_coordinategrid_eq():
    """
    Tests for the CoordinateGrid __eq__ method.
    """
    grid1 = RectangularGrid(Nlat=100, Nlon=200)
    grid2 = RectangularGrid(Nlat=100, Nlon=200)
    assert grid1 == grid2

    with pytest.raises(TypeError):
        _ = grid1 == 1

    grid3 = RectangularGrid(Nlat=101, Nlon=200)
    assert grid1 != grid3

    grid4 = RectangularGrid(Nlat=200, Nlon=100)
    assert grid1 != grid4
