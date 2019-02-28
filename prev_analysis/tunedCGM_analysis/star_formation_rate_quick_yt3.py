from mpi4py import MPI
import numpy as na
import sys

#from yt.funcs import *
#from yt.mods import *

from yt.config import ytcfg

#from yt.utilities.cosmology import EnzoCosmology
from yt.frontends.enzo.simulation_handling import EnzoCosmology

import yt

import pdb

ytcfg.set('yt','serialize','False')
ytcfg.set('yt','StoreParameterFiles','False')

rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size

dataset = sys.argv[1]
output_file = "sfrd.dat"

pf = yt.load(dataset)

rho_crit_now = pf.quan(2.77533874e11, "Msun * Mpc**-3 * h**-2") # Msun Mpc^-3 h^-2

if pf.cosmological_simulation:
    rho_crit_now *= (pf['CosmologyHubbleConstantNow']**2)
    
f_ej = pf.get_parameter('StarMassEjectionFraction',type=float)

if pf.cosmological_simulation:
    pf.get_parameter('CosmologyOmegaLambdaNow',type=float)
    pf.get_parameter('InitialTime',type=float)
    pf.get_parameter('InitialCycleNumber',type=int)
    ec = EnzoCosmology(HubbleConstantNow=100*pf.parameters['CosmologyHubbleConstantNow'],
                       OmegaMatterNow=pf.parameters['CosmologyOmegaMatterNow'],
                       OmegaLambdaNow=pf.parameters['CosmologyOmegaLambdaNow'],
                       InitialRedshift=pf.parameters['CosmologyInitialRedshift'])
    
time_points = pf.parameters['InitialCycleNumber']

if pf.cosmological_simulation:
    initial_time_cu = ec.ComputeTimeFromRedshift(pf.parameters['CosmologyInitialRedshift']) / ec.TimeUnits
else:
    initial_time_cu = 0.0
initial_time_cu = pf.quan(initial_time_cu, "Myr")

final_time_cu = pf.parameters['InitialTime']
final_time_cu = pf.quan(final_time_cu, "Myr")
time_step = (final_time_cu - initial_time_cu) / (time_points - 1)

if pf.cosmological_simulation:
    hubble_time_now = ec.ComputeTimeFromRedshift(0.0)

dt = time_step

star_creation_time = pf.arr([], "Myr")
star_dynamical_time = pf.arr([], "Myr")
star_particle_mass = pf.arr([], "Myr")
star_metal_fraction = pf.arr([], "Myr")

my_grids = na.arange(rank, pf.index.num_grids, size)

pbar = yt.get_pbar('Loading star data:',len(my_grids))

for q, grid in enumerate(my_grids):
    pbar.update(q)
    creation_time = pf.index.grids[grid]['creation_time'].in_units("Myr")
    stars = creation_time > 0
    if na.where(stars)[0].size > 0:
        star_creation_time = pf.arr(na.concatenate([star_creation_time, creation_time[stars].in_units("Myr")]), "Myr")
        star_dynamical_time = pf.arr(na.concatenate([star_dynamical_time,
                                                     pf.index.grids[grid]['dynamical_time'][stars].in_units("Myr")]), "Myr")
        star_particle_mass = pf.arr(na.concatenate([star_particle_mass,
                                                    pf.index.grids[grid]['particle_mass'][stars].in_units("Msun")]), "Msun")
        star_metal_fraction = na.concatenate([star_metal_fraction, pf.index.grids[grid]['metallicity_fraction'][stars]])
    del stars
    pf.index.grids[grid].clear_data()

pbar.finish()

del my_grids

if rank == 0:
    f = open(output_file,'w')
    f.write("#z\tLookback time\tComoving SFRD\tStars\tRho_star\tOmega_star\tComoving MFRD\tRho_star_metal\tOmega_star_metal\n")
    f.close()

# get initial star masses
xv1 = (final_time_cu - star_creation_time) / star_dynamical_time
mass_initial = star_particle_mass / (1 - f_ej * (1 - ((1 + xv1)*na.exp(-xv1))))
too_big = mass_initial > (star_particle_mass / (1 - f_ej))
mass_initial[too_big] = (star_particle_mass[too_big] / (1 - f_ej))
del too_big

current_time_cu = initial_time_cu
if pf.cosmological_simulation:
    all_mass = pf.quan(0.0, "Msun/Mpccm**3")
    all_metals = pf.quan(0.0, "Msun/Mpccm**3")
    
else:
    all_mass = pf.quan(0.0, "Msun/Mpc**3")
    all_metals = pf.quan(0.0, "Msun/Mpc**3")

print(all_mass)
    
while current_time_cu <= final_time_cu:

    if pf.cosmological_simulation:
        current_redshift = ec.ComputeRedshiftFromTime(current_time_cu * ec.TimeUnits)
        lookback_time = hubble_time_now - current_time_cu * ec.TimeUnits
    else:
        current_redshift = 0.0
        # This could probably be computed better somehow?
        lookback_time = current_time_cu

    if rank == 0: print("Calculating sfrd for z = %.4e." % (current_redshift))

    MPI.COMM_WORLD.Barrier()
    formed_stars = star_creation_time < current_time_cu
    total_stars = na.where(formed_stars)[0].size
    sfrd = pf.quan(0.0, "Msun")
    mfrd = pf.quan(0.0, "Msun")
    if total_stars > 0:
        xv1 = (current_time_cu - star_creation_time[formed_stars]) / star_dynamical_time[formed_stars]
        xv2 = (current_time_cu + dt - star_creation_time[formed_stars]) / star_dynamical_time[formed_stars]
        young = xv2 < 12.0
        if na.where(young)[0].size > 0:
            xv1 = xv1[young]
            xv2 = xv2[young]
            mass_form = mass_initial[formed_stars][young] * ((1 + xv1)*na.exp(-xv1) - (1 + xv2)*na.exp(-xv2))
            mask = mass_form > star_particle_mass[formed_stars][young]
            mass_form[mask] = star_particle_mass[formed_stars][young][mask]
            mass_form[mass_form < 0] = 0.0
            star_metal = star_metal_fraction[formed_stars][young]
            sfrd = mass_form.sum()
            mfrd = (mass_form * star_metal).sum()
    if pf.cosmological_simulation:
        print("Proc %d: stars: %d, sfrd: %e." % (rank, total_stars, 
                                                 (sfrd / (dt * (pf.length_unit.in_units('Mpccm')**3)))))
    else:
        print("Proc %d: stars: %d, sfrd: %e." % (rank, total_stars, 
                                                 (sfrd / (dt * (pf.length_unit.in_units('Mpc')**3)))))
        
    total_stars = MPI.COMM_WORLD.allreduce(total_stars, op=MPI.SUM)
    sfrd = MPI.COMM_WORLD.allreduce(sfrd, op=MPI.SUM)
    mfrd = MPI.COMM_WORLD.allreduce(mfrd, op=MPI.SUM)

    print(current_time_cu)
    print(sfrd)
    print(all_mass)
    
    if pf.cosmological_simulation:
        all_mass = all_mass + (sfrd / (pf.units['mpccm']**3))
        all_metals += (mfrd / (pf.units['mpccm']**3))
        omega_star = all_mass / rho_crit_now
        omega_star_metals = all_metals / rho_crit_now
        sfrd /=  (dt * pf.time_units['years'] * (pf.units['mpccm']**3))
        mfrd /=  (dt * pf.time_units['years'] * (pf.units['mpccm']**3))
    else:
        all_mass = all_mass + (sfrd / (pf.length_unit.in_units('Mpc')**3))
        all_metals += (mfrd / (pf.length_unit.in_units('Mpc')**3))
        omega_star = all_mass / rho_crit_now
        omega_star_metals = all_metals / rho_crit_now
        sfrd /=  dt.in_units("yr") #* (pf.length_unit.in_units('Mpc')**3))
        mfrd /=  (dt * (pf.length_unit.in_units('Mpc')**3))
        
    MPI.COMM_WORLD.Barrier() 
    if rank == 0:
        f = open(output_file,'a')
        f.write("%.6e %.6e %.6e %d %.6e %.6e %.6e %.6e %.6e\n" % (current_redshift, lookback_time, sfrd, total_stars, 
                                                                  all_mass, omega_star, mfrd, all_metals, omega_star_metals))
        f.close()
    current_time_cu += time_step
