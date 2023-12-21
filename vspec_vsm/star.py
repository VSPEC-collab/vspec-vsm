"""
The `VSPEC` star model is designed modularly to allow
for both simple and complex behaviors. Currently, it
is represented by a rectangular grid of points on the stellar
surface, each assigned an effective temperature. At any given
time, the model computes the surface coverage fractions of
each temperature visible to the observer, accounting for the
spherical geometry, limb darkening, and any occultation
by a transiting planet.

Once the surface coverage is calculated, the fractions of each component are
passed allong to VSPEC to compute the spectrum.

The attributes of the ``Star`` class describe the bulk properties of the star,
including radius, period, and the effective temperature of quiet photosphere.
Herein we refer to this temperature as the photosphere temperature to differentiate
it from the temperature of spots, faculae, or other sources of variability.
"""
from typing import Tuple, Union
import warnings
import numpy as np
from astropy import units as u
from astropy.units.quantity import Quantity

from vspec_vsm.coordinate_grid import CoordinateGrid
from vspec_vsm.helpers import (
    get_angle_between, proj_ortho,
    calc_circ_fraction_inside_unit_circle,
    clip_teff
)
from vspec_vsm.spots import SpotCollection, SpotGenerator
from vspec_vsm.faculae import FaculaCollection, FaculaGenerator, Facula
from vspec_vsm.flares import FlareCollection, FlareGenerator
from vspec_vsm.granules import Granulation

class InsufficientResolutionWarning(Warning):
    """Transit calculated but there is not
    enough grid resolution to do it accurately.
    """

class Star:
    """
    Star object representing a variable star.

    Parameters
    ----------
    Teff : astropy.units.Quantity
        Effective temperature of the stellar photosphere.
    radius : astropy.units.Quantity
        Stellar radius.
    period : astropy.units.Quantity
        Stellar rotational period.
    spots : SpotCollection
        Initial spots on the stellar surface.
    faculae : FaculaCollection
        Initial faculae on the stellar surface.
    distance : astropy.units.Quantity , default=1*u.pc
        Distance to the star.
    Nlat : int, default=500
        The number of latitude points on the stellar sufrace.
    Nlon : int, default=1000
        The number of longitude points on the stellar surface.
    gridmaker : CoordinateGrid, default=None
        A `CoordinateGrid` object to create the stellar sufrace grid.
    flare_generator : FlareGenerator, default=None
        Flare generator object.
    spot_generator : SpotGenerator, default=None
        Spot generator object.
    fac_generator : FaculaGenerator, default=None
        Facula generator object.
    ld_params : list, default=[0, 1, 0]
        Limb-darkening parameters.

    Attributes
    ----------
    Teff : astropy.units.Quantity
        Effective temperature of the stellar photosphere.
    radius : astropy.units.Quantity
        Stellar radius.
    distance : astropy.units.Quantity
        Distance to the star.
    period : astropy.units.Quantity
        Stellar rotational period.
    spots : SpotCollection
        Spots on the stellar surface.
    faculae : FaculaCollection
        Faculae on the stellar surface.
    gridmaker : CoordinateGrid
        Object to create the coordinate grid of the surface.
    map : astropy.units.Quantity
        Pixel map of the stellar surface.
    flare_generator : FlareGenerator
        Flare generator object.
    spot_generator : SpotGenerator
        Spot generator object.
    fac_generator : FaculaGenerator
        Facula generator object.
    u1 : float
        Limb-darkening parameter u1.
    u2 : float
        Limb-darkening parameter u2.
    """

    def __init__(self, teff: u.Quantity,
                 radius: u.Quantity,
                 period: u.Quantity,
                 spots: SpotCollection,
                 faculae: FaculaCollection,
                 grid_params: Union[int,Tuple[int, int]] = 1000,
                 flare_generator: FlareGenerator = None,
                 spot_generator: SpotGenerator = None,
                 fac_generator: FaculaGenerator = None,
                 granulation: Granulation = None,
                 u1: float = 0,
                 u2: float = 0,
                 rng: np.random.Generator = np.random.default_rng()
                 ):
        self.teff = teff
        self.radius = radius
        self.period = period
        self.spots:SpotCollection = spots
        self.faculae:FaculaCollection = faculae
        self.rng = rng
        self.gridmaker = CoordinateGrid.new(grid_params)
        self.grid_params = grid_params
        
        self.faculae.set_gridmaker(self.gridmaker)
        self.spots.set_gridmaker(self.gridmaker)

        if flare_generator is None:
            self.flare_generator = FlareGenerator.off(
                rng=self.rng
            )
        else:
            self.flare_generator = flare_generator
        self.flares = None

        if spot_generator is None:
            self.spot_generator = SpotGenerator.off(
                grid_params=grid_params,
                gridmaker=self.gridmaker,
                rng=self.rng
            )
        else:
            self.spot_generator = spot_generator
            try:
                if self.spot_generator.gridmaker != self.gridmaker:
                    self.spot_generator.gridmaker = self.gridmaker
            except TypeError:
                self.spot_generator.gridmaker = self.gridmaker
            

        if fac_generator is None:
            self.fac_generator = FaculaGenerator.off(
                grid_params=grid_params,
                gridmaker=self.gridmaker,
                rng=self.rng
            )
        else:
            self.fac_generator = fac_generator
            try:
                if self.fac_generator.gridmaker != self.gridmaker:
                    self.fac_generator.gridmaker = self.gridmaker
            except TypeError:
                self.fac_generator.gridmaker = self.gridmaker
        if granulation is None:
            self.granulation = Granulation.off(seed=0)
        else:
            self.granulation = granulation
        self.u1 = u1
        self.u2 = u2
        self.set_spot_grid()
        self.set_fac_grid()

    def set_spot_grid(self):
        """
        Set the gridmaker for each spot to be the same.
        """
        for spot in self.spots.spots:
            spot.gridmaker = self.gridmaker

    def set_fac_grid(self):
        """
        Set the gridmaker for each facula to be the same.
        """
        for fac in self.faculae.faculae:
            fac.gridmaker = self.gridmaker

    @property
    def map(self):
        """
        Create a map of the stellar surface based on spots.

        Returns
        -------
        pixelmap : astropy.units.Quantity , Shape(self.gridmaker.Nlon,self.gridmaker.Nlat)
            Map of stellar surface with effective temperature assigned to each pixel.
        """
        return self.spots.map_pixels(self.radius, self.teff)

    def age(self, time):
        """
        Age the spots and faculae on the stellar surface according
        to their own `age` methods. Remove the spots that have decayed.

        Parameters
        ----------
        time : astropy.units.Quantity 
            Length of time to age the features on the stellar surface.
            For most realistic behavior, `time` should be much less than
            spot or faculae lifetime.
        """
        self.spots.age(time)
        self.faculae.age(time)

    def add_spot(self, spot):
        """
        Add one or more spots to the stellar surface.

        Parameters
        ----------
        spot : StarSpot or sequence of StarSpot
            The `StarSpot` object(s) to add.
        """
        self.spots.add_spot(spot)

    def add_fac(self, facula):
        """
        Add one or more faculae to the stellar surface.

        Parameters
        ----------
        facula : Facula or sequence of Facula
            The Facula object(s) to add.

        """
        self.faculae.add_faculae(facula)

    def get_mu(self, lat0: u.Quantity, lon0: u.Quantity) -> np.ndarray:
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
        mu : np.ndarray
            An array of cos(x) where x is
            the angle from disk center.

        """
        return self.gridmaker.cos_angle_from_disk_center(lat0, lon0)

    def ld_mask(self, mu) -> np.ndarray:
        """
        Get a translucent mask based on limb darkeining parameters.

        Parameters
        ----------
        mu : np.ndarray
            The cosine of the angle from disk center

        Returns
        -------
        mask : np.ndarray
            The limb-darkened mask.

        Notes
        -----
        To account for apparent size effect of points on the
        stellar surface, we add a factor of 1 to `u1`. This way,
        in the (Lambertian) case of no-limb darkening, the user
        can set ``u1 = 0``, ``u2 = 0``
        """
        mask = 1 - (self.u1+1) * (1 - mu) - self.u2 * (1 - mu)**2
        behind_star = mu < 0.
        mask[behind_star] = 0
        return mask

    def ld_mask_for_plotting(self, mu) -> np.ndarray:
        """
        Same as above, but does not add the extra 1
        because that is already accounted for by the projection.
        """
        mask = 1 - (self.u1) * (1 - mu) - self.u2 * (1 - mu)**2
        behind_star = mu < 0.
        mask[behind_star] = 0
        return mask

    def get_jacobian(self) -> np.ndarray:
        """
        Get the relative area of each point.

        Returns
        -------
        jacobian : np.ndarray
            The area of each point
        """
        return self.gridmaker.area

    def add_faculae_to_map(
        self,
        lat0: u.Quantity,
        lon0: u.Quantity
    ):
        """
        Add the faculae to the surface map.

        Parameters
        ----------
        lat0 : astropy.units.Quantity
            The sub-observer latitude.
        lon0 : astropy.units.Quantity
            The sub-observer longitude.

        Returns
        -------
        teffmap : astropy.units.Quantity
            A temperature map of the surface
        """
        map_from_spots = self.map
        mu = self.get_mu(lat0, lon0)
        faculae: Tuple[Facula] = self.faculae.faculae
        for facula in faculae:
            angle = get_angle_between(lat0, lon0, facula.lat, facula.lon)
            inside_fac = facula.map_pixels(self.radius)
            if not np.any(inside_fac):  # the facula is too small
                pass
            else:
                fracs = facula.fractional_effective_area(angle)
                dteff_wall, dteff_floor = fracs.keys()
                frac = fracs[dteff_wall].value
                mu_of_fac_pix = mu[inside_fac]
                border_mu = np.percentile(mu_of_fac_pix, 100*frac)
                wall_pix = inside_fac & (mu <= border_mu)
                floor_pix = inside_fac & (mu > border_mu)
                teff_wall = clip_teff(dteff_wall + self.teff)
                teff_floor = clip_teff(dteff_floor + self.teff)
                map_from_spots[wall_pix] = teff_wall
                map_from_spots[floor_pix] = teff_floor
        return map_from_spots

    def get_pl_frac(
        self,
        angle_past_midtransit: u.Quantity,
        orbit_radius: u.Quantity,
        planet_radius: u.Quantity,
        inclination: u.Quantity
    ):
        """
        Get planet fraction

        Parameters
        ----------
        angle_past_midtransit : astropy.units.Quantity
            The phase of the planet past the 180 degree mid transit point.
        orbit_radius : astropy.units.Quantity
            The radius of the planet's orbit.
        radius : astropy.units.Quantity
            The radius of the planet.

        inclination : astropy.units.Quantity
            The inclination of the planet. 90 degrees is transiting.
        """
        x = (orbit_radius/self.radius * np.sin(angle_past_midtransit)
             ).to_value(u.dimensionless_unscaled)
        y = (orbit_radius/self.radius * np.cos(angle_past_midtransit)
             * np.cos(inclination)).to_value(u.dimensionless_unscaled)
        rad = (planet_radius/self.radius).to_value(u.dimensionless_unscaled)
        return 1-calc_circ_fraction_inside_unit_circle(x, y, rad)

    def get_transit_mask(
        self,
        lat0: u.Quantity,
        lon0: u.Quantity,
        orbit_radius: u.Quantity,
        radius: u.Quantity,
        phase: u.Quantity,
        inclination: u.Quantity
    ):
        """
        Get a mask describing which pixels are covered by a transiting planet.

        Parameters
        ----------
        lat0 : astropy.units.Quantity
            The sub-observer latitude.
        lon0 : astropy.units.Quantity
            The sub-observer longitude.
        orbit_radius : astropy.units.Quantity
            The radius of the planet's orbit.
        radius : astropy.units.Quantity
            The radius of the planet.
        phase : astropy.units.Quantity
            The phase of the planet. 180 degrees is mid transit.
        inclination : astropy.units.Quantity
            The inclination of the planet. 90 degrees is transiting.

        Returns
        -------
        mask : np.ndarray
            The fraction of each pixel that is covered by the planet.
        pl_frac : float
            The fraction of the planet that is visible to the observer.
        """
        eclipse = False
        if np.cos(phase) > 0:
            eclipse = True
        angle_past_midtransit: u.Quantity = phase - 180*u.deg
        x: float = (orbit_radius/self.radius * np.sin(angle_past_midtransit)
             ).to_value(u.dimensionless_unscaled)
        y: float = (orbit_radius/self.radius * np.cos(angle_past_midtransit)
             * np.cos(inclination)).to_value(u.dimensionless_unscaled)
        rp_rs: float = (radius/self.radius).to_value(u.dimensionless_unscaled)
        if np.sqrt(x**2 + y**2) > 1 + 2*rp_rs:  # no transit
            return self.gridmaker.zeros().astype('bool'), 1.0
        elif eclipse:
            planet_fraction = self.get_pl_frac(
                angle_past_midtransit, orbit_radius, radius, inclination)
            return self.gridmaker.zeros().astype('bool'), planet_fraction
        else:
            llat, llon = self.gridmaker.grid()
            xcoord, ycoord = proj_ortho(lat0, lon0, llat, llon)
            mu = self.gridmaker.cos_angle_from_disk_center(lat0, lon0)
            area = self.gridmaker.area  # area in units of 4pi steradians ie adds to 1
            point_radii = 2*np.sqrt(area)  # radius of each pixel in radians. The 2
                                           # comes from the 4 in 4pi steradians
            proj_radii = point_radii*mu  # radius of each pixel in projected coords
            # distances in projected coords
            rad_map = np.sqrt((xcoord-x)**2 + (ycoord-y)**2)
            # case 1: Point is completely outside transit radius
            case1 = (rad_map > rp_rs + 2*proj_radii) | np.isnan(rad_map)

            covered_value = np.where(~case1, 1, 0).astype('float')
            if np.any(np.isnan(covered_value)):
                raise ValueError('NaN in covered_value')
            indicies = np.argwhere(~case1)
            if len(indicies)>0:
                relevent_radii = proj_radii[~case1]
                rad_mean = np.mean(relevent_radii)
                if rad_mean > 0.5*rp_rs:
                    area_mean = np.pi*rad_mean**2
                    area_pl = np.pi*rp_rs**2
                    target_area = area_pl/4
                    factor = area_mean/target_area
                    
                    warnings.warn(
                        f'Pixel resolution too low. Increase by factor of {factor:.2f}',
                        InsufficientResolutionWarning
                    )
                    
            for index in indicies:
                s = index if index is int else tuple(index)
                dist_from_transit_center: float = rad_map[s]
                # gauss_sigma = proj_radii[s]
                sigma_x = point_radii[s]*mu[s]
                sigma_y = point_radii[s]
                if rp_rs > dist_from_transit_center+2*sigma_x:
                    covered_value[s] = 1.
                else:
                    x = np.linspace(-3*sigma_x, 3*sigma_x, 100)
                    xx, yy = np.meshgrid(x, x)
                    # zz = 1/(2*np.pi*gauss_sigma**2)*np.exp(-(xx**2 + yy**2)
                    #                                        / (2*gauss_sigma**2))
                    zz = 1/(2*np.pi*sigma_x*sigma_y)*np.exp(-0.5*(xx/sigma_x)**2 - 0.5*(yy/sigma_y)**2)
                    dist = np.sqrt((xx-dist_from_transit_center)**2 + (yy)**2)
                    overlap = np.where(dist < rp_rs, zz, 0)
                    overlap = np.trapz(overlap, x, axis=1)
                    overlap = np.trapz(overlap, x, axis=0)
                    covered_value[s] = overlap
            if np.any(np.isnan(covered_value)):
                raise ValueError('NaN in covered_value')
            return covered_value, 1.0

    def calc_coverage(
        self,
        sub_obs_coords: dict,
        granulation_fraction: float = 0.0,
        orbit_radius: u.Quantity = 1*u.AU,
        planet_radius: u.Quantity = 1*u.R_earth,
        phase: u.Quantity = 90*u.deg,
        inclination: u.Quantity = 0*u.deg

    ):
        """
        Calculate coverage

        Calculate coverage fractions of various Teffs on stellar surface
        given coordinates of the sub-observation point.

        Parameters
        ----------
        sub_obs_coord : dict
            A dictionary giving coordinates of the sub-observation point.
            This is the point that is at the center of the stellar disk from the view of
            an observer. Format: {'lat':lat,'lon':lon} where lat and lon are
            `astropy.units.Quantity` objects.
        granulation_fraction : float
            The fraction of the quiet photosphere that has a lower Teff due to granulation

        Returns
        -------
        total_data : dict
            Dictionary with Keys as Teff quantities and Values as surface fraction floats.
        covered_data : dict
            Dictionary with Keys as Teff quantities and Values as surface fraction floats covered
            by a transiting planet.
        pl_frac : float
            The fraction of the planet that is visble. This is in case of an eclipse.
        """
        cos_c = self.get_mu(sub_obs_coords['lat'], sub_obs_coords['lon'])
        ld = self.ld_mask(cos_c)
        jacobian = self.get_jacobian()

        surface_map = self.add_faculae_to_map(
            sub_obs_coords['lat'], sub_obs_coords['lon'])
        covered, pl_frac = self.get_transit_mask(
            sub_obs_coords['lat'], sub_obs_coords['lon'],
            orbit_radius=orbit_radius,
            radius=planet_radius,
            phase=phase,
            inclination=inclination
        )

        teffs = np.unique(surface_map)
        total_data = {}
        covered_data = {}
        total_area = np.sum(ld*jacobian)
        for teff in teffs:
            pix_has_teff = np.where(surface_map == teff, 1, 0)
            nominal_area = np.sum(pix_has_teff*ld*jacobian)
            covered_area = np.sum(pix_has_teff*ld*jacobian*(covered))
            total_data[teff] = nominal_area/total_area
            covered_data[teff] = covered_area/total_area
        granulation_teff = self.teff - self.granulation.dteff
        # initialize. This way it's okay if there's something else with that Teff too.
        if granulation_teff not in teffs:
            total_data[granulation_teff] = 0
            covered_data[granulation_teff] = 0

        phot_frac = total_data[self.teff]
        total_data[self.teff] = phot_frac * (1-granulation_fraction)
        total_data[granulation_teff] += phot_frac * granulation_fraction
        phot_frac = covered_data[self.teff]
        covered_data[self.teff] = phot_frac * (1-granulation_fraction)
        covered_data[granulation_teff] += phot_frac * granulation_fraction

        return total_data, covered_data, pl_frac

    def birth_spots(self, time):
        """
        Create new spots from a spot generator.

        Parameters
        ----------
        time : astropy.units.Quantity 
            Time over which these spots should be created.

        """
        self.spots.add_spot(self.spot_generator.birth_spots(time, self.radius))

    def birth_faculae(self, time):
        """
        Create new faculae from a facula generator.

        Parameters
        ----------
        time : astropy.units.Quantity 
            Time over which these faculae should be created.


        """
        self.faculae.add_faculae(
            self.fac_generator.birth_faculae(time, self.radius))

    def average_teff(self, sub_obs_coords):
        """
        Calculate the average Teff of the star given a sub-observation point
        using the Stephan-Boltzman law. This can approximate a lightcurve for testing.

        Parameters
        ----------
        sub_obs_coords : dict
            A dictionary containing coordinates of the sub-observation point. This is the 
            point that is at the center of the stellar disk from the view of an observer. 
            Format: {'lat':lat,'lon':lon} where lat and lon are `astropy.units.Quantity` objects.

        Returns
        -------
        astropy.units.Quantity 
            Bolometric average Teff of stellar disk.

        """
        dat, _, _ = self.calc_coverage(sub_obs_coords)
        num = 0
        den = 0
        for teff, value in dat.items():
            num += teff**4 * value
            den += value
        return ((num/den)**(0.25)).to(u.K)

    def plot_surface(
        self,
        lat0: u.Quantity,
        lon0: u.Quantity,
        ax = None,
        orbit_radius: u.Quantity = 1*u.AU,
        radius: u.Quantity = 1*u.R_earth,
        phase: u.Quantity = 90*u.deg,
        inclination: u.Quantity = 0*u.deg,
        nlon: int = 1000,
        nlat: int = 500
    ):
        """
        Add the transit to the surface map and plot.
        """
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        from cartopy.mpl.geoaxes import GeoAxes
        ax: GeoAxes
        proj = ccrs.Orthographic(
            central_latitude=lat0.to_value(u.deg),
            central_longitude=lon0.to_value(u.deg)
        )
        if ax is None:
            _, ax = plt.subplots(subplot_kw={'projection': proj})
            ax: GeoAxes = ax
        elif ax.projection != proj:
            ax.projection = (proj)
        covered, pl_frac = self.get_transit_mask(
            lat0, lon0,
            orbit_radius=orbit_radius,
            radius=radius,
            phase=phase,
            inclination=inclination
        )
        map_with_faculae = self.add_faculae_to_map(lat0, lon0)

        lats, lons, data = self.gridmaker.display_grid(
            nlat, nlon, map_with_faculae)

        lats = (lats*u.rad).to_value(u.deg)
        lons = (lons*u.rad).to_value(u.deg)
        data = data.to_value(u.K)
        im = ax.pcolormesh(lons, lats, data.T,
                           transform=ccrs.PlateCarree())
        plt.colorbar(im, ax=ax, label=r'$T_{\rm eff}$ (K)')

        _, _, data = self.gridmaker.display_grid(nlat, nlon, covered)

        transit_mask = np.where(data > 0.5, 1, np.nan)
        zorder = 100 if pl_frac == 1. else -100
        ax.contourf(lons, lats, transit_mask.T, colors='k', alpha=1,
                    transform=ccrs.PlateCarree(), zorder=zorder)
        mu = self.get_mu(lat0, lon0)
        ld = self.ld_mask_for_plotting(mu)
        alpha = 1-ld.T/np.max(ld.T)
        ax.imshow(np.ones_like(ld), extent=(lons.min(), lons.max(), lats.min(), lats.max()),
                  transform=ccrs.PlateCarree(), origin='lower', alpha=alpha, cmap=plt.cm.get_cmap('gray'), zorder=100)

    def get_flares_over_observation(self, time_duration: Quantity):
        """
        Generate a collection of flares over a specified observation period.

        Parameters
        ----------
        time_duration: astropy.units.Quantity 
            The duration of the observation period.

        Notes
        -----
        This method uses the `FlareGenerator` attribute of the `Star` object to generate
        a distribution of flare energies, and then generates a series of flares over the
        specified observation period using these energies. The resulting collection of
        flares is stored in the `Star`'s `flares` attribute.
        """
        flares = self.flare_generator.generate_flare_series(time_duration)
        self.flares = FlareCollection(flares)

    def get_flare_int_over_timeperiod(self, tstart: Quantity, tfinish: Quantity, sub_obs_coords):
        """
        Compute the total flare integral over a specified time period and sub-observer point.

        Parameters
        ----------
        tstart: astropy.units.Quantity 
            The start time of the period.
        tfinish: astropy.units.Quantity 
            The end time of the period.
        sub_obs_coords : dict
            A dictionary containing coordinates of the sub-observation point. This is the 
            point that is at the center of the stellar disk from the view of an observer. 
            Format: {'lat':lat,'lon':lon} where lat and lon are `astropy.units.Quantity` objects.


        Returns
        -------
        flare_timeareas: list of dict
            List of dictionaries containing flare temperatures and integrated
            time-areas. In the format [{'Teff':9000*u.K,'timearea'=3000*u.Unit('km2 hr)},...]

        Notes
        -----
        This method computes the total flare integral over each flare in the `flares` 
        attribute of the `Star` object that falls within the specified time period and is visible
        from the sub-observer point defined by `sub_obs_coords`. The result is returned
        as a list of dictionaries representing the teff and total flare integral of each flare.
        """
        flare_timeareas = self.flares.get_flare_integral_in_timeperiod(
            tstart, tfinish, sub_obs_coords)
        return flare_timeareas

    def generate_mature_spots(self, coverage: float):
        """
        Generate new mature spots with a specified coverage.

        Parameters
        ----------
        coverage: float
            The coverage of the new spots.

        Notes
        -----
        This method uses the `SpotGenerator` attribute of the current object to generate a
        set of new mature spots with the specified coverage. The new spots are added to 
        the object's `spots` attribute and the pixel map is updated using the new spots.
        """
        new_spots = self.spot_generator.generate_mature_spots(
            coverage, self.radius)
        self.spots.add_spot(new_spots)

    def get_granulation_coverage(self, time: u.Quantity) -> np.ndarray:
        """
        Calculate the coverage by granulation at each point in `time`.

        Parameters
        ----------
        time : astropy.units.Quantity
            The points on the time axis.

        Returns
        -------
        np.ndarray
            The coverage corresponding to each `time` point.
        """
        if self.granulation is None:
            return np.zeros(shape=time.shape)
        else:
            coverage = self.granulation.get_coverage(time)
            return np.where(np.isnan(coverage), 0, coverage)
