import pytest
import astropy.units as u
import numpy as np

from vspec_vsm.star import Star
from vspec_vsm.spots import SpotCollection, SpotGenerator
from vspec_vsm.faculae import FaculaCollection, FaculaGenerator
from vspec_vsm.flares import FlareGenerator

from vspec_vsm.granules import Granulation
from vspec_vsm.coordinate_grid import CoordinateGrid
from vspec_vsm.config import MSH

@pytest.fixture
def star():
    # Create an instance of the Star class for testing
    Teff = 3000 * u.K
    radius = 0.15 * u.Rsun
    period = 10 * u.day
    spots = SpotCollection()
    faculae = FaculaCollection()
    return Star(Teff, radius, period, spots, faculae,grid_params=(500,1000))

@pytest.fixture
def star_with_spots():
    # Create an instance of the Star class for testing
    # This one has spots
    Teff = 3000 * u.K
    radius = 0.15 * u.Rsun
    period = 10 * u.day
    faculae = FaculaCollection()
    sgen = SpotGenerator(500*MSH,0.2,2700*u.K,2600*u.K)
    spots = SpotCollection(*sgen.birth_spots(10*u.day,radius))
    return Star(Teff, radius, period, spots, faculae,
                spot_generator=sgen,grid_params=1000)
@pytest.fixture
def star_with_spots_and_fac():
    # Create an instance of the Star class for testing
    # This one has spots
    Teff = 3000 * u.K
    radius = 0.15 * u.Rsun
    period = 10 * u.day
    
    sgen = SpotGenerator(500*MSH,0.2,2700*u.K,2600*u.K)
    spots = SpotCollection(*sgen.birth_spots(10*u.day,radius))
    fgen = FaculaGenerator(
            dist_r_peak=100*u.km,
            dist_r_logsigma=0.2,
            depth=100*u.km,
            dist_life_peak=6*u.hr,
            dist_life_logsigma=0.2,
            floor_teff_slope=0*u.K/u.km,
            floor_teff_min_rad=20*u.km,
            floor_teff_base_dteff=-100*u.K,
            wall_teff_slope=0*u.K/u.km,
            wall_teff_intercept=100*u.K,
            coverage=0.01,
            grid_params=(300, 600),
            gridmaker=CoordinateGrid.new((300, 600)),
            rng=np.random.default_rng()
        )
    faculae = FaculaCollection(*fgen.birth_faculae(10*u.hr,radius))
    return Star(Teff, radius, period, spots, faculae,
                spot_generator=sgen,fac_generator=fgen)


def test_star_initialization(star:Star):
    assert star.teff == 3000 * u.K
    assert star.radius == 0.15 * u.Rsun
    assert star.period == 10 * u.day
    assert isinstance(star.spots, SpotCollection)
    assert isinstance(star.faculae, FaculaCollection)
    assert isinstance(star.gridmaker, CoordinateGrid)
    assert isinstance(star.map, u.Quantity)
    assert isinstance(star.flare_generator, FlareGenerator)
    assert isinstance(star.spot_generator, SpotGenerator)
    assert isinstance(star.fac_generator, FaculaGenerator)
    assert isinstance(star.granulation,Granulation)
    assert star.u1 == 0
    assert star.u2 == 0


def test_get_pixelmap(star:Star):
    pixelmap = star.map

    assert isinstance(pixelmap, u.Quantity)
    assert pixelmap.shape == (star.gridmaker.nlon, star.gridmaker.nlat)
    assert pixelmap.unit == u.K

def test_star_with_granulation():
    Teff = 5000 * u.K
    radius = 1 * u.Rsun
    period = 10 * u.day
    spots = SpotCollection()
    faculae = FaculaCollection()
    granulation = Granulation(0.2,0.01,5*u.day,200*u.K)
    star = Star(Teff, radius, period, spots, faculae, granulation=granulation)

    assert isinstance(star.granulation, Granulation)


def test_star_with_custom_gridmaker():
    Teff = 5000 * u.K
    radius = 1 * u.Rsun
    period = 10 * u.day
    spots = SpotCollection()
    faculae = FaculaCollection()
    star = Star(Teff, radius, period, spots, faculae,grid_params=1000)

    assert isinstance(star.gridmaker, CoordinateGrid)


def test_age_method_updates_pixelmap(star_with_spots:Star):
    # Store the initial pixelmap
    initial_pixelmap = star_with_spots.map
    

    # Age the spots and faculae
    time = 4 * u.day
    star_with_spots.age(time)

    # Check if the pixelmap has been updated
    assert not np.array_equal(star_with_spots.map, initial_pixelmap)

    # Make additional assertions based on the expected behavior of the method



def test_age_method_with_zero_time(star_with_spots:Star):
    # Add some spots and faculae to the star
    initial_num_spots = len(star_with_spots.spots.spots)
    initial_num_faculae = len(star_with_spots.faculae.faculae)

    # Age the spots and faculae with zero time
    time = 0 * u.day
    star_with_spots.age(time)

    # Check if the spots and faculae remain unchanged
    assert len(star_with_spots.spots.spots) == initial_num_spots
    assert len(star_with_spots.faculae.faculae) == initial_num_faculae

def test_transit_mask():
    Teff = 3000 * u.K
    radius = 0.15 * u.Rsun
    period = 10 * u.day
    spots = SpotCollection()
    faculae = FaculaCollection()
    star = Star(Teff, radius, period, spots, faculae,grid_params=(500,1000))
    
    mask, val = star.get_transit_mask(
        lat0=0*u.deg,
        lon0=0*u.deg,
        orbit_radius=0.05*u.AU,
        radius = 1*u.R_earth,
        phase=180*u.deg,
        inclination=90*u.deg
    )
    assert val == 1
    assert np.any(mask == 0)
    assert np.max(mask) == pytest.approx(1,abs=1e-2)
    assert np.any((mask<1) & (mask>0))
    
    
    # star-sized planet
    # mask, val = star.get_transit_mask(
    #     lat0=0*u.deg,
    #     lon0=0*u.deg,
    #     orbit_radius=1*u.AU,
    #     radius = star.radius,
    #     phase=180*u.deg,
    #     inclination=90*u.deg
    # )
    # jac = star.get_jacobian().to_value(u.dimensionless_unscaled)
    # mu = star.get_mu(lat0=0*u.deg,lon0=0*u.deg).to_value(u.dimensionless_unscaled)
    # facing = mu > 0
    # area = jac*mu
    # area_facing = np.sum(area[facing])
    # area_blocked = np.sum((mask[facing])*area[facing])
    # assert area_blocked == pytest.approx(area_facing,rel=1e-6)
    
    # rp/rs = 0.5
    # x = [8,32,128]
    # y = []
    # phases = (np.linspace(0,1,20)+180)*u.deg
    # t0 = datetime.now()
    # for rs_rp in x:
    #     # rs_rp = 4
    #     lat0 = 0*u.deg
    #     lon0=0*u.deg
    #     a = 0.05*u.AU
    #     rad = star.radius/rs_rp
    #     phase = 180.*u.deg
    #     i = 90*u.deg
    #     yy = []
    #     for phase in phases:
    #         mask, val = star.get_transit_mask(
    #             lat0=lat0,
    #             lon0=lon0,
    #             orbit_radius=a,
    #             radius = rad,
    #             phase=phase,
    #             inclination=i
    #         )
    #         jac = star.get_jacobian().to_value(u.dimensionless_unscaled)
    #         mu = star.get_mu(lat0=0*u.deg,lon0=0*u.deg).to_value(u.dimensionless_unscaled)
    #         facing = mu > 0
    #         area = jac*mu
    #         area_facing = np.sum(area[facing])
    #         area_blocked = np.sum((mask[facing])*area[facing])
    #         expected_area_blocked = area_facing/rs_rp**2
    #         extra_area_blocked = area_blocked - expected_area_blocked
    #         extra_area_blocked_div_exp = extra_area_blocked/expected_area_blocked
    #         yy.append(extra_area_blocked_div_exp*100)
    #         # unity = rs_rp**2*area_blocked/area_facing
    #         # y.append(100*(unity-1))
    #     plt.scatter(phases,yy,label=f'{rs_rp}')
    # # plt.scatter(x,y)
    # plt.xlabel('phase')
    # plt.ylabel('Extra Area blocked (%)')
    # # plt.ylabel('Divergence from expected (%)')
    # # plt.xscale('log')
    # # plt.xticks(x,labels=x)
    # plt.title(datetime.now()-t0)
    # plt.savefig('transit_bad.png',facecolor='w')
    
        
def test_calc_coverage():
    Teff = 3000 * u.K
    radius = 0.15 * u.R_sun
    period = 10 * u.day
    spots = SpotCollection()
    faculae = FaculaCollection()
    star = Star(Teff, radius, period, spots, faculae,grid_params=(500,1000),
                u1=0.1)
    
    rs_rp = 10
    lat0 = 0*u.deg
    lon0=0*u.deg
    sub_obs = {'lat':lat0,'lon':lon0}
    a = 0.05*u.AU
    rad = star.radius/rs_rp
    phase = 180.*u.deg
    i = 90*u.deg

    dat_tot, dat_cov, pl_frac = star.calc_coverage(
        sub_obs_coords=sub_obs,
        granulation_fraction=0.0
    )
    assert dat_tot[f'{star.teff:.2f}'] == 1.0, 'Total coverage must be 1.0'
    assert dat_cov[f'{star.teff:.2f}'] == 0.0, 'Occulted coverage must be 0.0'
    assert pl_frac == 1
    
    dat_tot, dat_cov, pl_frac = star.calc_coverage(
        sub_obs_coords=sub_obs,
        granulation_fraction=0.0,
        orbit_radius=a,
        planet_radius=rad,
        phase=phase,
        inclination=i
    )
    assert dat_tot[f'{star.teff:.2f}'] == 1.0, 'Total coverage must be 1.0'
    # assert dat_cov[star.Teff] == pytest.approx(1/rs_rp**2,rel=1e-3), 'Planet occultation incorrect for phase=90deg'
    assert pl_frac == 1, 'All of planet must be visible.'
    
    dat_tot, dat_cov, pl_frac = star.calc_coverage(
        sub_obs_coords=sub_obs,
        granulation_fraction=0.0,
        orbit_radius=a,
        planet_radius=rad,
        phase=phase+0.3*u.deg,
        inclination=i
    )
    assert dat_tot[f'{star.teff:.2f}'] == 1.0, 'Total coverage must be 1.0'
    # assert dat_cov[star.Teff] == pytest.approx(1/rs_rp**2,rel=1e-3), 'Planet occultation incorrect for phase=90.3deg'
    assert pl_frac == 1, 'All of planet must be visible.'
    
def test_transit_coverage():
    rp_rs = np.logspace(-2.5,-1,10)
    Teff = 3000 * u.K
    radius = 0.15 * u.R_sun
    period = 10 * u.day
    spots = SpotCollection()
    faculae = FaculaCollection()
    star = Star(Teff, radius, period, spots, faculae,grid_params=1000000)
    
    lat0 = 0*u.deg
    lon0=0*u.deg
    sub_obs = {'lat':lat0,'lon':lon0}
    a = 0.05*u.AU
    phase = 180.*u.deg
    i = 90*u.deg
    
    for r in rp_rs:
        dat_tot, dat_cov, pl_frac = star.calc_coverage(
            sub_obs_coords=sub_obs,
            granulation_fraction=0.0,
            orbit_radius=a,
            planet_radius=r*star.radius,
            phase=phase,
            inclination=i
        )
        assert dat_tot[f'{star.teff:.2f}'] == 1.0, 'Total coverage must be 1.0'
        assert dat_cov[f'{star.teff:.2f}'] == pytest.approx(r**2,rel=5e-2), f'Planet occultation incorrect for rp/rs = {r}'
        assert pl_frac == 1, 'All of planet must be visible.'


# def test_add_fac_to_map():
#     Teff = 3000 * u.K
#     radius = 0.15 * u.Rsun
#     period = 10 * u.day
#     sgen = SpotGenerator(5000*MSH,0.2,2900*u.K,2800*u.K,coverage=0.4,init_area=500*MSH)
#     spots = SpotCollection(*sgen.birth_spots(10*u.day,radius))
#     fgen = FaculaGenerator(dist_r_peak = 0.1*u.R_sun,coverage=0.4)
#     # facula = Facula(0*u.deg,40*u.deg,0.1*u.R_sun,0.02*u.R_sun,2700*u.K,3700*u.K,1*u.hr,True,Zw=10000*u.km)
#     faculae = FaculaCollection(*fgen.birth_faculae(50*u.hr,radius,Teff))
#     # faculae = FaculaCollection(facula)
#     star = Star(Teff, radius, period, spots, faculae)
    
#     # fac:Facula = star.faculae.faculae[0]
#     # inside_fac = fac.map_pixels(radius)
#     lon0,lat0 = 0*u.deg,0*u.deg
#     # angle = get_angle_between(lat0,lon0,fac.lat,fac.lon)
#     # fracs = fac.fractional_effective_area(angle)
#     # N_fac_pix = np.sum(inside_fac)
#     # N_wall = N_fac_pix * list(fracs.values())[0]
#     # N_floor = N_fac_pix * list(fracs.values())[1]

#     # mu = star.get_mu(lat0,lon0)
#     # mu_of_fac_pix = mu[inside_fac]
#     # border_mu = np.percentile(mu_of_fac_pix,100*list(fracs.values())[1].value)
#     # wall_pix = inside_fac & (mu >= border_mu)
#     # floor_pix = inside_fac & (mu < border_mu)
#     # wall_teff,floor_teff = fracs.keys()
#     # map = star.map
#     # map[wall_pix] = wall_teff
#     # map[floor_pix] = floor_teff

   

#     lats,lons = star.gridmaker.oned()
#     proj = ccrs.Orthographic(
#                 central_longitude=0,central_latitude=0)
#     fig = plt.figure()
#     gs = fig.add_gridspec(1,1)
#     ax = fig.add_subplot(gs[0,0],projection=proj, fc="r")
#     im = ax.pcolormesh(lons.value,lats.value,star.add_faculae_to_map(lat0,lon0).T,transform=ccrs.PlateCarree(),cmap='viridis')
#     # im = ax.pcolormesh(lons.value,lats.value,map.T.value,transform=ccrs.PlateCarree(),cmap='viridis')
#     fig.colorbar(im,ax=ax)
#     fig.show()
#     0


# def test_transit_mask():
#     Teff = 3000 * u.K
#     radius = 0.15 * u.R_sun
#     period = 10 * u.day
#     sgen = SpotGenerator(1000*MSH,0.2,2900*u.K,2800*u.K,coverage=0.4,init_area=500*MSH,distribution='iso')
#     spots = SpotCollection(*sgen.birth_spots(100*u.day,radius))
#     fgen = FaculaGenerator(dist_r_peak = 0.1*u.R_sun,coverage=0.4)
#     # facula = Facula(0*u.deg,40*u.deg,0.1*u.R_sun,0.02*u.R_sun,2700*u.K,3700*u.K,1*u.hr,True,Zw=10000*u.km)
#     faculae = FaculaCollection(*fgen.birth_faculae(0*u.hr,radius,Teff))
#     # faculae = FaculaCollection(facula)
#     star = Star(Teff, radius, period, spots, faculae,Nlon=1000,Nlat=500)
#     pmap = star.map.value

#     lon0,lat0 = 0*u.deg,0*u.deg
    
#     inclination = 90.*u.deg
#     phase = 0.8*u.deg#+180*u.deg
#     planet_radius = 0.1*radius
#     semimajor_axis = 0.05*u.AU
#     lat,lon = star.gridmaker.oned()
#     # llat,llon = star.gridmaker.grid()
#     # llon = llon-180*u.deg
#     # x,y = proj_ortho(lat0,lon0,llat,llon)

#     # x_pl = (semimajor_axis/radius * np.sin(phase)).to_value(u.dimensionless_unscaled)
#     # y_pl = (semimajor_axis/radius * np.cos(phase)* np.cos(inclination)).to_value(u.dimensionless_unscaled)
#     # rad_pl = (planet_radius/radius).to_value(u.dimensionless_unscaled)
#     # rad_map = np.sqrt((x-x_pl)**2 + (y-y_pl)**2)
#     transited,pl_frac = star.get_transit_mask(lat0,lon0,semimajor_axis,planet_radius,phase,inclination)
#     visible = np.where(~transited,np.nan,0)



#     # plt.scatter(x,y,c=pmap,s=10,cmap='viridis')
#     # plt.colorbar()

#     fig = plt.figure()
#     gs = fig.add_gridspec(1,2)
#     proj = ccrs.Orthographic(
#                 central_longitude=lon0.to_value(u.deg),central_latitude=lat0.to_value(u.deg))
#     ax = fig.add_subplot(gs[0,0],projection=proj, fc="r")
#     im = ax.pcolormesh(lon.value,lat.value,pmap.T,transform=ccrs.PlateCarree(),cmap='viridis')
#     zorder = 100 if pl_frac == 1. else -100
#     ax.contourf(lon.value,lat.value,visible.T,transform=ccrs.PlateCarree(),colors='k',alpha=1,zorder=zorder)
#     fig.colorbar(im,ax=ax)
#     # ax2 = fig.add_subplot(gs[0,1])
#     # ax2.set_aspect('equal')
#     # im = ax2.scatter(x,y,c=pmap,s=5,cmap='viridis')
#     # fig.colorbar(im,ax=ax2)
#     fig.show()

# def test_calc_coverage():
#     Teff = 3000 * u.K
#     radius = 0.15 * u.R_sun
#     period = 10 * u.day
#     sgen = SpotGenerator(1000*MSH,0.2,2900*u.K,2800*u.K,coverage=0.4,init_area=500*MSH,distribution='iso')
#     spots = SpotCollection(*sgen.birth_spots(100*u.day,radius))
#     fgen = FaculaGenerator(dist_r_peak = 0.1*u.R_sun,coverage=0.4)
#     # facula = Facula(0*u.deg,40*u.deg,0.1*u.R_sun,0.02*u.R_sun,2700*u.K,3700*u.K,1*u.hr,True,Zw=10000*u.km)
#     faculae = FaculaCollection(*fgen.birth_faculae(0*u.hr,radius,Teff))
#     # faculae = FaculaCollection(facula)
#     star = Star(Teff, radius, period, spots, faculae,Nlon=1000,Nlat=500,granulation=Granulation(0,0,5*u.day,200*u.K))
    
#     inclination = 90.*u.deg
#     phase = 0.*u.deg+180*u.deg
#     planet_radius = 0.1*radius
#     semimajor_axis = 0.05*u.AU
#     lon0,lat0 = 0*u.deg,0*u.deg
#     total,covered,pl_frac = star.calc_coverage(
#         {'lat':lat0,'lon':lon0},0.0,semimajor_axis,planet_radius,phase,inclination
#     )
#     0

# def test_plot_surface():
#     Teff = 3000 * u.K
#     radius = 0.15 * u.R_sun
#     period = 10 * u.day
#     sgen = SpotGenerator(1000*MSH,0.2,2900*u.K,2800*u.K,coverage=0.4,init_area=500*MSH,distribution='iso')
#     spots = SpotCollection(*sgen.birth_spots(100*u.day,radius))
#     fgen = FaculaGenerator(dist_r_peak = 0.2*u.R_sun,coverage=0.4)
#     # facula = Facula(0*u.deg,40*u.deg,0.1*u.R_sun,0.02*u.R_sun,2700*u.K,3700*u.K,1*u.hr,True,Zw=10000*u.km)
#     faculae = FaculaCollection(*fgen.birth_faculae(10*u.hr,radius,Teff))
#     for fac in faculae.faculae:
#         fac.Zw = 0.01*u.R_sun
#     # faculae = FaculaCollection(facula)
#     star = Star(Teff, radius, period, spots, faculae,Nlon=1000,Nlat=500,u1=0.4,u2=-1)
#     lon0,lat0 = 0*u.deg,0*u.deg
#     inclination = 90.*u.deg
#     phase = -0.0*u.deg+180*u.deg
#     planet_radius = 0.1*radius
#     semimajor_axis = 0.05*u.AU
#     proj = ccrs.Orthographic(
#             central_latitude=lat0.to_value(u.deg),
#             central_longitude=lon0.to_value(u.deg)
#     )
#     fig,axes = plt.subplots(1,1,subplot_kw={'projection': proj },figsize=(5,4))
#     # for i,ax in enumerate(axes):
#     star.plot_surface(lat0,lon0,axes,semimajor_axis,planet_radius,
#                         phase,inclination)
#     0





if __name__ in '__main__':
    pytest.main(args=[__file__])