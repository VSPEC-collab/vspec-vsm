"""
Profiling for Faculae
"""

from line_profiler import LineProfiler
from astropy import units as u
import numpy as np

from VSPEC import params

from vspec_vsm.star import Star
from vspec_vsm.faculae import Facula

def init_star():
    star_params = params.StarParameters(
        psg_star_template='M',
        teff=3300*u.K,
        mass=0.15*u.Msun,
        radius=0.15*u.Rsun,
        period=10*u.day,
        misalignment=0*u.deg,
        misalignment_dir=0*u.deg,
        ld = params.LimbDarkeningParameters(1.0,0.0),
        spots=params.SpotParameters.none(),
        faculae=params.FaculaParameters.std(),
        flares=params.FlareParameters.none(),
        granulation=params.GranulationParameters.none(),
        Nlat=500,Nlon=1000
    )
    return Star.from_params(star_params,rng=np.random.default_rng(10),seed=10)
    
def age_star(star:Star,dtime:u.Quantity,n_steps:int):
    for _ in range(n_steps):
        star.age(dtime)

def calc_coverage(star:Star):
    return star.calc_coverage(sub_obs_coords={'lat':0*u.deg,'lon':0*u.deg})

def sim_observation():
    star = init_star()
    n_steps = 3
    dtime = 1*u.hr
    for _ in range(n_steps):
        _=calc_coverage(star)
        star.birth_faculae(dtime)
        star.age(dtime)


if __name__ in '__main__':
    profiler = LineProfiler()
    profiler.add_function(sim_observation)
    profiler.add_function(Facula.map_pixels)
    
    profiler_wrapper = profiler(sim_observation)
    
    profiler_wrapper()
    
    profiler.print_stats(sort='time')
    