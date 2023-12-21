"""
Coordinate Grid class
"""
from typing import Tuple, Union
import warnings

import numpy as np
from astropy import units as u


def get_lat_points(n_points: int) -> u.Quantity:
    """
    Get a 1D array of latitudes.

    Parameters
    ----------
    n_points : int
        The number of points.

    Returns
    -------
    lats : astropy.units.Quantity
        Array of latitudes.
    """
    lats = np.linspace(-90, 90, n_points)*u.deg
    return lats


def get_lon_points(n_points: int) -> u.Quantity:
    """
    Get a 1D array of longitudes.

    Parameters
    ----------
    n_points : int
        The number of points.

    Returns
    -------
    lons : astropy.units.Quantity
        Array of longitudes.
    """
    lons = np.linspace(0, 360, n_points, endpoint=False)*u.deg
    return lons


class CoordinateGrid:
    """
    Base class for all coordinate grids.
    """

    def _grid(self) -> Tuple[u.Quantity, u.Quantity]:
        raise NotImplementedError(
            'Attempted to call abstract method _grid() from base class')
    @staticmethod
    def new(grid_params: Union[int,Tuple[int, int]]):
        """
        Create an instance of a coordinate grid subclass
        from a set of parameters.
        
        Parameters
        ----------
        grid_params : int or tuple
            If a tuple, the first element is the number of latitude
            points and the second is the number of longitude points.
            If an int, it is the total number of points.
        
        Returns
        -------
        grid : SpiralGrid or RectangularGrid
            An instance of a coordinate grid subclass. If grid_params
            is a tuple, the ``grid`` is a ``RectangularGrid``. If
            ``grid_params`` is an int, the ``grid`` is a ``SpiralGrid``.
        
        Raises
        ------
        TypeError
            If grid_params is not an int or a tuple.
        """
        if isinstance(grid_params, int):
            return SpiralGrid(grid_params)
        elif isinstance(grid_params, tuple):
            return RectangularGrid(*grid_params)
        else:
            raise TypeError('grid_params must be of type int or tuple')
    
    def grid(self) -> Tuple[u.Quantity, u.Quantity]:
        """
        Get a grid of latitudes and longitudes.

        Returns
        -------
        lat : astropy.units.Quantity
            Array of latitudes.
        lon : astropy.units.Quantity
            Array of longitudes.
        """
        return self._grid()

    def cos_angle_from_disk_center(
        self,
        lat0: u.Quantity,
        lon0: u.Quantity
    ) -> np.ndarray:
        """
        Get the cosine of the angle from disk center.

        Parameters
        ----------
        lat0 : astropy.units.Quantity
            The sub-observer latitude.
        lon0 : astropy.units.Quantity
            The sub-observer longitude

        Returns
        -------
        np.ndarray
            An array of cos(x) where x is
            the angle from disk center.

        Notes
        -----
        Recall
        
        .. math::

            \\mu = \\cos{x}

        Where :math:`x` is the angle from center of the disk.
        """
        latgrid, longrid = self.grid()
        mu: u.Quantity = (np.sin(lat0) * np.sin(latgrid)
                          + np.cos(lat0) * np.cos(latgrid)
                          * np.cos(lon0-longrid))
        return mu.to_value(u.dimensionless_unscaled)

    @property
    def _area(self):
        raise NotImplementedError(
            'Attempted to call abstract method _area() from base class')

    @property
    def area(self) -> np.ndarray:
        """
        Get the area of each point as a fraction of the unit sphere.

        Returns
        -------
        np.ndarray
            The area of each point.

        """
        return self._area

    def _zeros(self, dtype='float32'):
        raise NotImplementedError(
            'Attempted to call abstract method _zeros() from base class')

    def zeros(self, dtype: str = 'float32') -> np.ndarray:
        """
        Get an array of zeros.

        Parameters
        ----------
        dtype : str, default='float32'
            Data type to pass to np.zeros.

        Returns
        -------
        np.ndarray
            An array of zeros.

        """
        return self._zeros(dtype=dtype)

    def __eq__(self, other):
        raise NotImplementedError(
            'Attempted to call abstract method __eq__() from base class')

    def _display_grid(self, nlat: int, nlon: int, data: np.ndarray):
        raise NotImplementedError(
            'Attempted to call abstract method _display_grid() from base class')

    def display_grid(
        self,
        nlat: int,
        nlon: int,
        data: Union[np.ndarray, u.Quantity]
    ) -> Tuple[np.ndarray, np.ndarray, u.Quantity]:
        """
        Resample the data to a rectangular grid.

        Parameters
        ----------
        nlat : int
            Number of latitude points to resample to.
        nlon : int
            Number of longitude points to resample to.
        data : np.ndarray or astropy.units.Quantity
            The data to display.

        Returns
        -------
        lat : np.ndarray, shape=(nlat,)
            Array of latitudes.
        lon : np.ndarray, shape=(nlon,)
            Array of longitudes.
        dat : astropy.units.Quantity, shape=(nlat,nlon)
            The resampled data.

        """
        data = data*u.dimensionless_unscaled
        unit: u.Unit = data.unit
        lat, lon, dat = self._display_grid(
            nlat, nlon, data=data.to_value(unit))
        return lat, lon, dat*unit


class RectangularGrid(CoordinateGrid):
    """
    Class to standardize the creation of latitude and longitude grids.

    This class provides a convenient way to create latitude
    and longitude grids of specified dimensions. It allows
    the creation of both one-dimensional arrays and two-dimensional
    grids of latitude and longitude points.

    Parameters
    ----------
    Nlat : int, optional (default=500)
        Number of latitude points.
    Nlon : int, optional (default=1000)
        Number of longitude points.

    Raises
    ------
    TypeError
        If Nlat or Nlon is not an integer.


    Attributes
    ----------
    Nlat : int
        Number of latitude points.
    Nlon : int
        Number of longitude points.

    Examples
    --------
    >>> grid = CoordinateGrid(Nlat=100, Nlon=200)
    >>> lats, lons = grid.oned()
    >>> print(lats.shape, lons.shape)
    (100,) (200,)
    >>> grid_arr = grid.grid()
    >>> print(grid_arr.shape)
    (100, 200)
    >>> zeros_arr = grid.zeros()
    >>> print(zeros_arr.shape)
    (200, 100)
    >>> other_grid = CoordinateGrid(Nlat=100, Nlon=200)
    >>> print(grid == other_grid)
    True

    """

    def __init__(self, nlat=500, nlon=1000):
        if not isinstance(nlat, int):
            raise TypeError('Nlat must be int')
        if not isinstance(nlon, int):
            raise TypeError('Nlon must be int')
        self.nlat = nlat
        self.nlon = nlon

    def oned(self):
        """
        Create one dimensional arrays of latitude and longitude points.

        Returns
        -------
        lats : astropy.units.Quantity , shape=(Nlat,)
            Array of latitude points.
        lons : astropy.units.Quantity , shape=(Nlon,)
            Array of longitude points.

        """
        lats = get_lat_points(self.nlat)
        lons = get_lon_points(self.nlon)
        return lats, lons

    def _grid(self):
        """
        Create a 2 dimensional grid of latitude and longitude points.

        Returns
        -------
        lats : astropy.units.Quantity , shape=(Nlat,Nlon)
            Array of latitude points.
        lons : astropy.units.Quantity , shape=(Nlat,Nlon)
            Array of longitude points.

        """
        lats, lons = self.oned()
        return np.meshgrid(lats, lons)

    def _zeros(self, dtype='float32'):
        """
        Get a grid of zeros.

        Parameters
        ----------
        dtype : str, default='float32'
            Data type to pass to np.zeros.

        Returns
        -------
        arr : np.ndarray, shape=(Nlon, Nlat)
            Grid of zeros.

        """
        return np.zeros(shape=(self.nlon, self.nlat), dtype=dtype)

    def __eq__(self, other):
        """
        Check to see if two CoordinateGrid objects are equal.

        Parameters
        ----------
        other : CoordinateGrid
            Another CoordinateGrid object.

        Returns
        -------
        bool
            Whether the two objects have equal properties.

        Raises
        ------
        TypeError
            If `other` is not a CoordinateGrid object.

        """
        if not isinstance(other, RectangularGrid):
            raise TypeError('other must be of type RectangularGrid')
        else:
            return (self.nlat == other.nlat) & (self.nlon == other.nlon)

    @property
    def _area(self) -> np.ndarray:
        """
        Get the area of each point.

        Returns
        -------
        np.ndarray
            The area of each point.

        """
        latgrid, _ = self.grid()
        jacobian: u.Quantity = np.sin(latgrid + 90*u.deg)
        norm: u.Quantity = np.sum(jacobian)
        return (jacobian/norm).to_value(u.dimensionless_unscaled)

    @property
    def dlat(self):
        """
        The latitude spacing.

        Returns
        -------
        dlat : astropy.units.Quantity
            The latitude spacing.
        """
        return 180*u.deg/(self.nlat-1)

    @property
    def dlon(self):
        """
        The longitude spacing.

        Returns
        -------
        dlon : astropy.units.Quantity
            The longitude spacing.
        """
        return 360*u.deg/(self.nlon)

    def _display_grid(
        self,
        nlat: int,
        nlon: int,
        data: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Check to make sure the grid is correct. Otherwise
        there is little to do as the grid is already rectangular.
        """
        # checks
        if not data.ndim == 2:
            raise TypeError('data must be 2 dimensional')
        if not data.shape[0] == nlon:
            raise TypeError('data must have shape (nlon, nlat)')
        if not data.shape[1] == nlat:
            raise TypeError('data must have shape (nlon, nlat)')
        lats, lons = self.oned()
        return lats.to_value(u.rad), lons.to_value(u.rad), data


class SpiralGrid(CoordinateGrid):
    """
    A class to generate a grid of points
    using a Fibonacci spiral.
    """
    GOLDEN_RATIO = 0.5*(1+np.sqrt(5))

    def __init__(self, n_points: int):
        self.n_points = n_points

    def _grid(self):
        """
        Produce a 1D grid of latitudes and longitudes.

        Returns
        -------
        lat : astropy.units.Quantity
            Array of latitudes.
        lon : astropy.units.Quantity
            Array of longitudes.

        Notes
        -----
        Algorithm from:
            https://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
        """
        i = np.arange(self.n_points)
        lon = 2*np.pi*i/self.GOLDEN_RATIO
        lon = (lon % (2*np.pi))*u.rad
        colat = np.arccos(1-2*(i+0.5)/self.n_points)
        lat = (np.pi/2 - colat)*u.rad
        return lat, lon

    def _zeros(self, dtype='float32') -> np.ndarray:
        return np.zeros(shape=(self.n_points,), dtype=dtype)

    @property
    def _area(self) -> np.ndarray:
        """
        The area as a fraction of the unit sphere.
        """
        return self.zeros() + 1/self.n_points

    def __eq__(self, other) -> bool:
        if isinstance(other, CoordinateGrid):
            if isinstance(other, SpiralGrid):
                return self.n_points == other.n_points
            else:
                warnings.warn(
                    'Comparing two different subclasses of CoordinateGrid. Evaluating to False.')
                return False
        else:
            raise TypeError('other must be of type CoordinateGrid')

    def _display_grid(
        self,
        nlat: int,
        nlon: int,
        data: np.ndarray
    ):
        """
        Resample to a rectangular grid using a gaussian.
        """
        # checks
        if not data.ndim == 1:
            raise ValueError('data must be a 1D array')
        if not data.shape[0] == self.n_points:
            raise ValueError('data must have length n_points')
        data_lat, data_lon = self.grid()
        rect = RectangularGrid(nlat, nlon)

        numerator = np.zeros((nlon, nlat))
        denominator = np.zeros((nlon, nlat))
        # number density in points per steradian
        num_density = self.n_points/(4*np.pi)
        characteristic_len = 1/np.sqrt(num_density)  # characteristic length

        for value, lat, lon in zip(data, data_lat, data_lon):
            cos_r = rect.cos_angle_from_disk_center(lat, lon)
            r = np.arccos(cos_r)
            weight = np.exp(-(r/characteristic_len)**2)
            numerator += weight*value
            denominator += weight
        resampled_data = numerator/denominator
        lats, lons = rect.oned()
        return lats.to_value(u.rad), lons.to_value(u.rad), resampled_data
