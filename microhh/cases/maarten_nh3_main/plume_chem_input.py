import matplotlib.pyplot as pl
import numpy as np
import pandas as pd
import netCDF4 as nc

import microhh_tools as mht
import helpers as hlp
from constants import *
from lsm_input import LSM_input

"""
Settings.
"""
float_type = "f8"

xsize = 6400
ysize = 3200
zsize = 4000

itot = 64
jtot = 32
ktot = 96

start_hour = 6
end_hour = 12

# Enable resolved plume rise:
sw_plume_rise = False

# Enable non-linear KPP chemistry:
sw_chemistry = True

# Enable land-surface model and more detailled deposition.
sw_land_surface = True


"""
Read base .ini file for case settings.
"""
ini = mht.Read_namelist("plume_chem.ini.base")


"""
Create case input.
"""
if sw_chemistry:

    # Read TUV output table.
    # There is a total of 24 hour available, generated for the
    # Jaenschwalde power plant on 23/05/2022.
    columns = [
        "time",
        "sza",
        "jo31d",
        "jh2o2",
        "jno2",
        "jno3",
        "jn2o5",
        "jch2or",
        "jch2om",
        "jch3o2h",
    ]
    tuv = pd.read_table(
        "plume_chem_tuv_output.txt",
        sep="\\s+",
        skiprows=12,
        skipfooter=1,
        engine="python",
        names=columns,
        index_col="time",
    )

    # NOTE: `.loc` is value based, not on index.
    tuv = tuv.loc[start_hour:end_hour]

    # Convert to seconds, and subtract starting time.
    tuv.index *= 3600
    tuv.index -= tuv.index.values[0]

    # Emissions (?)
    emi_no = np.zeros(tuv.index.size)
    emi_isop = np.zeros(tuv.index.size)

    # Concentrations, for now constant with height.
    species = {"nh3": 1e-9}

    deposition_species = ["nh3"]
else:
    species = {}


"""
Define vertical grid and initial thl/qt/co2 profiles.
"""
dz = zsize / ktot
z = np.arange(0.5 * dz, zsize, dz)


def profile(zi, v_bulk, dv, gamma_v, clip_at_zero=False):
    """
    Create well mixed profile with jump and constant lapse rate above.
    """

    k_zi = np.abs(z - zi).argmin()

    profile = np.zeros(ktot)
    profile[:k_zi] = v_bulk
    profile[k_zi:] = v_bulk + dv + gamma_v * (z[k_zi:] - zi)

    if clip_at_zero:
        profile[profile < 0] = 0.0

    return profile


# Vertical profiles.
thl = profile(zi=1000, v_bulk=290, dv=2, gamma_v=0.006)
# qt  = profile(zi=1000, v_bulk=10e-3, dv=-2e-3, gamma_v=-0.003e-3, clip_at_zero=True)
qt = profile(zi=1000, v_bulk=6e-3, dv=-2e-3, gamma_v=-0.002e-3, clip_at_zero=True)
u = np.ones(ktot) * 5

# Surface fluxes.
t0 = start_hour * 3600
t1 = end_hour * 3600
td = 12 * 3600
time = np.linspace(t0, t1, 32)

# wthl = 0.15 * np.sin(np.pi * (time-t0) / td)
# wqt  = 8e-5 * np.sin(np.pi * (time-t0) / td)
wthl = 0.20 * np.sin(np.pi * (time - t0) / td)
wqt = 6e-5 * np.sin(np.pi * (time - t0) / td)

# Surface radiation (only used with land-surface enabled).
sw_flux_dn = 600 * np.sin(np.pi * (time - t0) / td)
sw_flux_dn[sw_flux_dn < 0] = 0

sw_flux_up = 0.2 * sw_flux_dn

lw_flux_dn = np.ones_like(sw_flux_dn) * 340
lw_flux_up = np.ones_like(sw_flux_dn) * 400

pl.figure(figsize=(10, 5))

pl.subplot(131)
pl.plot(time / 3600, wthl)
pl.xlabel("time (h)")
pl.ylabel("w`thl` (K m s-1)")

pl.subplot(132)
pl.plot(time / 3600, wqt * 1000)
pl.xlabel("time (h)")
pl.ylabel("w`qt` (g kg-1 m s-1)")

pl.subplot(133)
pl.plot(time / 3600, sw_flux_dn, label="sw_flux_dn")
pl.plot(time / 3600, sw_flux_up, label="sw_flux_up")
pl.plot(time / 3600, lw_flux_dn, label="lw_flux_dn")
pl.plot(time / 3600, lw_flux_up, label="lw_flux_up")
pl.xlabel("time (h)")
pl.ylabel("sw_flux_dn` (W m-2)")
pl.legend()

pl.tight_layout()

"""
Write input NetCDF file.
"""


def add_nc_var(name, dims, nc, data):
    if dims is None:
        var = nc.createVariable(name, np.float64)
    else:
        var = nc.createVariable(name, np.float64, dims)
    var[:] = data


nc_file = nc.Dataset("plume_chem_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", ktot)
add_nc_var("z", ("z"), nc_file, z)

# Atmospheric input.
nc_group_init = nc_file.createGroup("init")

add_nc_var("u", ("z"), nc_group_init, u)
add_nc_var("thl", ("z"), nc_group_init, thl)
add_nc_var("qt", ("z"), nc_group_init, qt)
# add_nc_var('co2_inflow', ('z'), nc_group_init, co2)

nc_tdep = nc_file.createGroup("timedep")
nc_tdep.createDimension("time_surface", time.size)

add_nc_var("time_surface", ("time_surface"), nc_tdep, time - time[0])
add_nc_var("thl_sbot", ("time_surface"), nc_tdep, wthl)
add_nc_var("qt_sbot", ("time_surface"), nc_tdep, wqt)

if sw_chemistry:
    # Chemistry input.
    nc_chem = nc_file.createGroup("timedep_chem")
    nc_chem.createDimension("time_chem", tuv.index.size)

    add_nc_var("time_chem", ("time_chem"), nc_chem, tuv.index)
    add_nc_var("jo31d", ("time_chem"), nc_chem, tuv.jo31d)
    add_nc_var("jh2o2", ("time_chem"), nc_chem, tuv.jh2o2)
    add_nc_var("jno2", ("time_chem"), nc_chem, tuv.jno2)
    add_nc_var("jno3", ("time_chem"), nc_chem, tuv.jno3)
    add_nc_var("jn2o5", ("time_chem"), nc_chem, tuv.jn2o5)
    add_nc_var("jch2or", ("time_chem"), nc_chem, tuv.jch2or)
    add_nc_var("jch2om", ("time_chem"), nc_chem, tuv.jch2om)
    add_nc_var("jch3o2h", ("time_chem"), nc_chem, tuv.jch3o2h)
    add_nc_var("emi_isop", ("time_chem"), nc_chem, emi_isop)
    add_nc_var("emi_no", ("time_chem"), nc_chem, emi_no)

    for name, value in species.items():
        profile = np.ones(ktot, dtype=np.float64) * value
        add_nc_var(name, ("z"), nc_group_init, profile)
        add_nc_var("{}_inflow".format(name), ("z"), nc_group_init, profile)

if sw_land_surface:
    nc_soil = nc_file.createGroup("soil")
    nc_soil.createDimension("z", 4)
    add_nc_var("z", ("z"), nc_soil, np.array([-1.945, -0.64, -0.175, -0.035]))

    add_nc_var("theta_soil", ("z"), nc_soil, np.array([0.34, 0.25, 0.21, 0.18]))
    add_nc_var("t_soil", ("z"), nc_soil, np.array([282, 287, 290, 286]))
    add_nc_var("index_soil", ("z"), nc_soil, np.ones(4) * 2)
    add_nc_var("root_frac", ("z"), nc_soil, np.array([0.05, 0.3, 0.4, 0.25]))

    # Add idealized (prescribed) radiation.
    add_nc_var("sw_flux_dn", ("time_surface"), nc_tdep, sw_flux_dn)
    add_nc_var("sw_flux_up", ("time_surface"), nc_tdep, sw_flux_up)
    add_nc_var("lw_flux_dn", ("time_surface"), nc_tdep, lw_flux_dn)
    add_nc_var("lw_flux_up", ("time_surface"), nc_tdep, lw_flux_up)

nc_file.close()


"""
Define emissions.
"""
# Coordinates of central cooling tower (m):
x0 = 200
y0 = ysize / 2.0

# Std-dev of plume widths:
sigma_x = 50
sigma_y = 50
sigma_z = 21
z0 = 0

## if sw_plume_rise:
##     # Emissions from tower height.
##     z0 = 120
##     sigma_z = 25
## else:
##     # The heights and sigma are from Dominik's .csv profiles, curve fitted with Python.
##     z0 = 299.68  # of 599.69 for high
##     sigma_z = 122.37
## 
## # x,y spacing towers:
## dx = 290
## dy = 120
## ddx = 40

## Strength of plumes, from the CoCO2 simulation protocol:
#strength_co2 = 732.5 / 9.0 / MCO2  # kmol(CO2) s-1
#strength_no2 = 0.0289 / 9.0 / MNO2  # kmol(NO2) s-1
#strength_no = 0.359 / 9.0 / MNO  # kmol(NO) s-1
#strength_co = 0.223 / 9.0 / MCO  # kmol(CO) s-1
# strength_nh3  = 1.0  / 9. / MNH3    # kmol(NH3) s-1
# strength_nh3  = 1.48950934102587E-06  # kmol(NH3) s-1 (equivalent to one barn of 80 cows)
strength_nh3 = 0.0

# Emission of heat and moisture. Numbers are from:
# Effective pollutant emission heights for atmospheric transport modelling based on real-world information
# Pregger and Friedrich, 2009, 10.1016/j.envpol.2008.09.027
Tp = 50 + T0  # Flue gass temperature (K)
Mp = 790  # Volume-flux (m3 s-1)

# This is not very accurate...:
pp = 1e5
rhop = pp / (Rd * Tp)
rhp = 1.0
qp = rhp * hlp.calc_qsat(Tp, pp)

strength_q = np.round(Mp * rhop * qp, decimals=2)
strength_T = np.round(Mp * rhop * Tp, decimals=2)

# Emission input model.
emi = hlp.Emissions()

# Add emission from individual towers.
# NOTE: Jaenschwalde has nine cooling towers, in three groups.
#       In reality, only two groups are active at the same time.
#       Here, we distribute the emissions over all nine towers......
##for j in range(-1, 2):
##    for i in range(-1, 2):
##        x = x0 + i * dx + j * ddx
##        y = y0 + j * dy
##        z = z0
##


x = x0
y = y0
z = z0

##        # emi.add('co2', strength_co2, True, x, y, z, sigma_x, sigma_y, sigma_z)

if sw_chemistry:
    emi.add("nh3", strength_nh3, True, x, y, z, sigma_x, sigma_y, sigma_z)

if sw_plume_rise:
    emi.add("thl", strength_T, False, x, y, z, sigma_x, sigma_y, sigma_z)
    emi.add("qt", strength_q, False, x, y, z, sigma_x, sigma_y, sigma_z)

"""
Create heterogeneous land-surface with three adjacent zones
"""
if (sw_land_surface):
    blocksize_i = 4
    blocksize_j = 4
    
    # Create three separate masks
    mask1 = np.zeros((jtot, itot), dtype=bool)
    mask2 = np.zeros((jtot, itot), dtype=bool)
    mask3 = np.zeros((jtot, itot), dtype=bool)
    
    # Divide domain into thirds in x-direction
    third_size = itot // 3
    
    for j in range(jtot):
        for i in range(itot):
            if i < third_size:
                mask1[j,i] = True
            elif i < 2*third_size:
                mask2[j,i] = True
            else:
                mask3[j,i] = True
    
    ls = LSM_input(itot, jtot, 4, sw_water=True, TF=float_type, debug=True, exclude_fields=['z0m', 'z0h'])
    
    # Set properties for each mask zone
    # First zone (mask1)
    ls['c_veg'][mask1] = 0.0
    ls['lai'][mask1] = 0.0
    ls['water_mask'][mask1] = 0.0
    
    # Second zone (mask2)
    ls['c_veg'][mask2] = 1.0
    ls['lai'][mask2] = 6.0
    ls['water_mask'][mask2] = 0.0
    
    # Third zone (mask3)
    ls['c_veg'][mask3] = 1.0
    ls['lai'][mask3] = 0.1
    ls['water_mask'][mask3] = 0
    
    # Constant properties
    ls['gD'][:,:] = 0.6
    ls['rs_veg_min'][:,:] = 300
    ls['rs_soil_min'][:,:] = 200
    ls['lambda_stable'][:,:] = 20
    ls['lambda_unstable'][:,:] = 20
    ls['cs_veg'][:,:] = 0.5
    ls['t_bot_water'][:,:] = 298
    
    ls['t_soil'][:,:,:] = 295
    ls['theta_soil'][:,:,:] = 0.1
    ls['index_soil'][:,:,:] = 0
    ls['root_frac'][:,:,:] = 0.15
    
    # Check if all values are set
    ls.check()
    
    # Save to binary and NetCDF
    ls.save_binaries(allow_overwrite=True)
    ls.save_netcdf('lsm_input.nc', allow_overwrite=True)

# """
# Create heterogeneous land-surface, here with simple block pattern to look at deposition differences.
# """
# if (sw_land_surface):
##
# blocksize_i = 8
# blocksize_j = 8
##
# mask = np.zeros((jtot, itot), dtype=bool)
# mask[:] = False
##
# for j in range(jtot):
# for i in range(itot):
# patch_i = i // blocksize_i % 2 == 0
# patch_j = j // blocksize_j % 2 == 0
##
# if (patch_i and patch_j) or (not patch_i and not patch_j):
# mask[j,i] = True
##
# ls = LSM_input(itot, jtot, 4, sw_water=True, TF=float_type, debug=True, exclude_fields=['z0m', 'z0h'])
##
# def set_value(variable, masked, non_masked):
# ls[variable][ mask] = masked
# ls[variable][~mask] = non_masked
##
# Spatially varying properties.
# set_value('c_veg',      masked=0.8, non_masked=0)
# set_value('lai'  ,      masked=2,   non_masked=0)
# set_value('lai'  ,      masked=0.2,   non_masked=0)
# set_value('water_mask', masked=0,   non_masked=0)
# set_value('water_mask', masked=0,   non_masked=1)
##
# Constant properties.
# ls['gD'][:,:] = 0.5
# ls['rs_veg_min'][:,:] = 100
# ls['rs_soil_min'][:,:] = 50
# ls['lambda_stable'][:,:] = 10
# ls['lambda_unstable'][:,:] = 10
# ls['cs_veg'][:,:] = 0
# ls['t_bot_water'][:,:] = 295
##
# ls['t_soil'][:,:,:] = 290
# ls['theta_soil'][:,:,:] = 0.1
# ls['index_soil'][:,:,:] = 0
# ls['root_frac'][:,:,:] = 0.25
##
# Check if all values are set.
# ls.check()
##
# Save to binary (used by MicroHH) and NetCDF (just for checking).
# ls.save_binaries(allow_overwrite=True)
# ls.save_netcdf('lsm_input.nc', allow_overwrite=True)


"""
Add settings to .ini file.
"""

ini["grid"]["itot"] = itot
ini["grid"]["jtot"] = jtot
ini["grid"]["ktot"] = ktot

ini["grid"]["xsize"] = xsize
ini["grid"]["ysize"] = ysize
ini["grid"]["zsize"] = zsize

# scalars = ['co2'] + list(species.keys())
scalars = list(species.keys())
ini["advec"]["fluxlimit_list"] = scalars
ini["limiter"]["limitlist"] = scalars
ini["fields"]["slist"] = scalars
ini["boundary"]["scalar_outflow"] = scalars

ini["time"]["endtime"] = (end_hour - start_hour) * 3600

ini["source"]["sourcelist"] = emi.source_list

ini["chemistry"]["swchemistry"] = sw_chemistry

crosslist = ["thl", "qt", "u", "v", "w", "thl_fluxbot", "qt_fluxbot"]

if sw_chemistry:
    # Add chemicial species and their vertical integrals.
    crosslist += list(species.keys())
    crosslist += [f"{x}_path" for x in species.keys()]
    crosslist += [f"vd{x}" for x in deposition_species]

if sw_land_surface and sw_chemistry:
    # Add deposition for each land-surface tile.
    for s in deposition_species:
        for t in ["soil", "wet", "veg"]:
            crosslist.append(f"vd{s}_{t}")

if sw_land_surface:
    ini["boundary"]["swboundary"] = "surface_lsm"
    ini["boundary"]["sbcbot"] = "flux"
    ini["boundary"]["sbot"] = "0"
    ini["boundary"]["thl"] = "dirichlet"
    ini["boundary"]["qt"] = "dirichlet"
    ini["boundary"]["swtimedep"] = False
    ini["boundary"]["timedeplist"] = "empty"

    ini["radiation"]["swradiation"] = "prescribed"
else:
    ini["boundary"]["swboundary"] = "surface"
    ini["boundary"]["sbcbot"] = "flux"
    ini["boundary"]["swtimedep"] = True
    ini["boundary"]["timedeplist"] = ["thl_sbot", "qt_sbot"]

    ini["radiation"]["swradiation"] = False

if sw_chemistry and sw_land_surface:
    ini["deposition"]["swdeposition"] = True
else:
    ini["deposition"]["swdeposition"] = False

ini["cross"]["crosslist"] = crosslist
ini["cross"]["xz"] = ysize / 2

ini["source"]["source_x0"] = emi.x0
ini["source"]["source_y0"] = emi.y0
ini["source"]["source_z0"] = emi.z0

ini["source"]["sigma_x"] = emi.sigma_x
ini["source"]["sigma_y"] = emi.sigma_y
ini["source"]["sigma_z"] = emi.sigma_z

ini["source"]["strength"] = emi.strength
ini["source"]["swvmr"] = emi.sw_vmr

ini["source"]["line_x"] = emi.line_x
ini["source"]["line_y"] = emi.line_y
ini["source"]["line_z"] = emi.line_z

ini.save("plume_chem.ini", allow_overwrite=True)

