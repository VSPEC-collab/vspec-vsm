"""
Tests for the CoordinateGrid class.
"""

import pytest
import numpy as np
from astropy import units as u

from vspec_vsm.coordinate_grid import RectangularGrid, SpiralGrid

def test_rectangulargrid():
    """
    Tests for the CoordinateGrid class.
    """
    grid = RectangularGrid(nlat=100, nlon=200)
    assert grid.nlat == 100
    assert grid.nlon == 200


def test_coodinategrid_typeerror():
    """
    Tests for the CoordinateGrid class initialized with floats.
    """
    with pytest.raises(TypeError):
        _ = RectangularGrid(nlat=100, nlon=200.0)
    with pytest.raises(TypeError):
        _ = RectangularGrid(nlat=100.0, nlon=200)


def test_rectangulargrid_oned():
    """
    Tests for the CoordinateGrid oned method.
    """
    grid = RectangularGrid(nlat=100, nlon=200)
    lats, lons = grid.oned()
    assert lats.shape == (100,)
    assert lons.shape == (200,)


def test_rectangulargrid_grid():
    """
    Tests for the CoordinateGrid grid method.
    """
    grid = RectangularGrid(nlat=100, nlon=200)
    lats, lons = grid.grid()
    assert lats.shape == (200, 100)
    assert lons.shape == (200, 100)


def test_rectangulargrid_zeros():
    """
    Tests for the CoordinateGrid zeros method.
    """
    grid = RectangularGrid(nlat=100, nlon=200)
    arr = grid.zeros()
    assert arr.shape == (200, 100)
    assert arr.dtype == np.float32
    assert np.all(arr == 0)


def test_rectangulargrid_eq():
    """
    Tests for the CoordinateGrid __eq__ method.
    """
    grid1 = RectangularGrid(nlat=100, nlon=200)
    grid2 = RectangularGrid(nlat=100, nlon=200)
    assert grid1 == grid2

    with pytest.raises(TypeError):
        _ = grid1 == 1

    grid3 = RectangularGrid(nlat=101, nlon=200)
    assert grid1 != grid3

    grid4 = RectangularGrid(nlat=200, nlon=100)
    assert grid1 != grid4

def test_rectangulargrid_display():
    """
    Tests for the CoordinateGrid display method.
    """
    grid = RectangularGrid(nlat=100, nlon=200)
    lat,lons = grid.grid()
    data = np.sin(2*lat)*np.cos(2*lons)
    llat,llon,dat = grid.display_grid(100,200,data)
    assert llat.shape == (100,)
    assert llon.shape == (200,)
    assert dat.shape == (200, 100)
    assert isinstance(dat, u.Quantity)
    
    llat,llon,dat = grid.display_grid(100,200,data*u.K)
    assert llat.shape == (100,)
    assert llon.shape == (200,)
    assert dat.shape == (200, 100)
    assert isinstance(dat, u.Quantity)
    assert dat.unit == u.K
    

def test_spiralgrid():
    """
    Tests for the SpiralGrid class.
    """
    grid = SpiralGrid(n_points=1000)
    assert grid.n_points == 1000

def test_spiralgrid_grid():
    """
    Tests for the SpiralGrid grid method.
    """
    grid = SpiralGrid(n_points=1000)
    lats, lons = grid.grid()
    assert lats.shape == (1000,)
    assert lons.shape == (1000,)
    assert np.all(lons.value<2*np.pi)
    assert np.all(lons.value>=0.)
    assert np.all(lats.value<=0.5*np.pi)
    assert np.all(lats.value>=-0.5*np.pi)

def test_spiralgrid_zeros():
    """
    Tests for the SpiralGrid zeros method.
    """
    grid = SpiralGrid(n_points=1000)
    arr = grid.zeros()
    assert arr.shape == (1000,)
    assert arr.dtype == np.float32
    assert np.all(arr == 0)

def test_spiralgrid_eq():
    """
    Tests for the SpiralGrid __eq__ method.
    """
    grid1 = SpiralGrid(n_points=1000)
    grid2 = SpiralGrid(n_points=1000)
    assert grid1 == grid2

    with pytest.raises(TypeError):
        _ = grid1 == 1

def test_spiralgrid_display():
    """
    Tests for the SpiralGrid display method.
    """
    grid = SpiralGrid(n_points=1000)
    lat,lons = grid.grid()
    data = np.sin(2*lat)*np.cos(2*lons)
    nlat = 100
    nlon = 200
    llat,llon,dat = grid.display_grid(nlat,nlon,data)
    assert llat.shape == (nlat,)
    assert llon.shape == (nlon,)
    assert dat.shape == (nlon,nlat)
    assert isinstance(dat, u.Quantity)
    assert dat.unit == u.dimensionless_unscaled
    
    llat,llon,dat = grid.display_grid(nlat,nlon,data*u.K)
    assert llat.shape == (nlat,)
    assert llon.shape == (nlon,)
    assert dat.shape == (nlon,nlat)
    assert isinstance(dat, u.Quantity)
    assert dat.unit == u.K