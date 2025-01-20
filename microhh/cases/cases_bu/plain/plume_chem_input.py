# Scenario: Forest on the right side of the Grassland 
# Domain Size: ... km (x) × ... km (y) × ... km (z)
# Species: NH3 (Background Concentration: 1 mol/mol)
# Point Source: None

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

# EDITED: Updated domain sizes as requested
xsize = 3200   # Changed from 3200
ysize = 3200   # Changed from 3200
zsize = 4000   # Changed from 3200

# EDITED: Updated grid points as requested
itot = 32     # Changed from 32
jtot = 32      # Changed from 32
ktot = 96     # Changed from 32

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


##########################################################
# Surface radiation (only used with land-surface enabled).
##########################################################
sw_flux_dn = 600 * np.sin(np.pi * (time-t0) / td)
sw_flux_dn[sw_flux_dn < 0] = 0

sw_flux_up = 0.2 * sw_flux_dn

lw_flux_dn = np.ones_like(sw_flux_dn) * 340
lw_flux_up = np.ones_like(sw_flux_dn) * 400
##########################################################

##pl.figure(figsize=(10,5))
##
##pl.subplot(131)
##pl.plot(time/3600, wthl)
##pl.xlabel("time (h)")
##pl.ylabel("w`thl` (K m s-1)")
##
##pl.subplot(132)
##pl.plot(time/3600, wqt*1000)
##pl.xlabel("time (h)")
##pl.ylabel("w`qt` (g kg-1 m s-1)")
##
##pl.subplot(133)
##pl.plot(time/3600, sw_flux_dn, label="sw_flux_dn")
##pl.plot(time/3600, sw_flux_up, label="sw_flux_up")
##pl.plot(time/3600, lw_flux_dn, label="lw_flux_dn")
##pl.plot(time/3600, lw_flux_up, label="lw_flux_up")
##pl.xlabel("time (h)")
##pl.ylabel("sw_flux_dn` (W m-2)")
##pl.legend()
##
##pl.tight_layout()

################################
#Write input NetCDF file.
################################

# Defines helper function add_nc_var to add variables to NetCDF file:
def add_nc_var(name, dims, nc, data):
    if dims is None:
        var = nc.createVariable(name, np.float64)
    else:
        var = nc.createVariable(name, np.float64, dims)
    var[:] = data

# Creates "plume_chem_input.nc" file
nc_file = nc.Dataset("plume_chem_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

###############################
# Sets up dimensions and groups
###############################

# "z" dimension for vertical levels:
nc_file.createDimension("z", ktot)
add_nc_var("z", ("z"), nc_file, z)

# Atmospheric input.
# ("init" group for initial atmospheric conditions)
nc_group_init = nc_file.createGroup("init")

add_nc_var("u", ("z"), nc_group_init, u)
add_nc_var("thl", ("z"), nc_group_init, thl)
add_nc_var("qt", ("z"), nc_group_init, qt)
#add_nc_var("co2", ("z"), nc_group_init, co2)
#add_nc_var("co2_inflow", ("z"), nc_group_init, co2)

#("timedep" group for time-dependent surface conditions)
nc_tdep = nc_file.createGroup("timedep");
nc_tdep.createDimension("time_surface", time.size)

add_nc_var("time_surface", ("time_surface"), nc_tdep, time-time[0])
add_nc_var("thl_sbot", ("time_surface"), nc_tdep, wthl)
add_nc_var("qt_sbot", ("time_surface"), nc_tdep, wqt)

if (sw_chemistry):
    # Chemistry input.
    # "timedep_chem" group (if chemistry enabled)
    nc_chem = nc_file.createGroup("timedep_chem");
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
        profile = np.ones(ktot, dtype=np.float64)*value
        add_nc_var(name, ("z"), nc_group_init, profile)
        add_nc_var("{}_inflow".format(name), ("z"), nc_group_init, profile)

    # Add flux_nh3 variable
    add_nc_var("flux_nh3", ("z"), nc_group_init, np.zeros(ktot, dtype=np.float64))

if (sw_land_surface):
    # "soil" group (if land surface enabled)
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
y0 = ysize/2.

# # Std-dev of plume widths:
# sigma_x = 25
# sigma_y = 25
#sigma_x = (xsize/itot)*0.5
#sigma_y = (ysize/jtot)*0.5
#sigma_z = (zsize/ktot)*0.5
#z0 = (zsize/ktot)*0.5
sigma_x = 50
sigma_y = 50
sigma_z = 21
z0 = 0

# # Handles plume rise conditions:
# # If enabled: uses tower height (120m)
# # If disabled: uses fitted heights from CSV profiles
# 
# if sw_plume_rise:
#     # Emissions from tower height.
#     z0 = 120
#     sigma_z = 25
# else:
#     # The heights and sigma are from Dominik"s .csv profiles, curve fitted with Python.
#     z0 = 299.68    # of 599.69 for high
#     sigma_z = 122.37

# # x,y spacing towers:
# dx = 290
# dy = 120
# ddx = 40

# # Strength of plumes, from the CoCO2 simulation protocol:
# strength_co2 = 0.0  / 9. / MCO2   # kmol(CO2) s-1
# strength_no2 = 0.0 / 9. / MNO2   # kmol(NO2) s-1
# strength_no  = 0.0  / 9. / MNO    # kmol(NO) s-1
# strength_co  = 0.0  / 9. / MCO    # kmol(CO) s-1
# strength_nh3  = 1.0  / 9. / MNH3    # kmol(NH3) s-1
# strength_nh3  = 1.48950934102587E-06  # kmol(NH3) s-1 (equivalent to one barn of 80 cows)
strength_nh3  = 0.0

# Emission of heat and moisture. Numbers are from:
# Effective pollutant emission heights for atmospheric transport modelling based on real-world information
# Pregger and Friedrich, 2009, 10.1016/j.envpol.2008.09.027
Tp = 50+T0    # Flue gass temperature (K)
Mp = 790      # Volume-flux (m3 s-1)

# This is not very accurate...:
pp = 1e5
rhop = pp/(Rd*Tp)
rhp = 1.
qp = rhp * hlp.calc_qsat(Tp, pp)

strength_q = np.round(Mp * rhop * qp, decimals=2)
strength_T = np.round(Mp * rhop * Tp, decimals=2)

# Emission input model.
emi = hlp.Emissions()

x = x0
y = y0
z = z0

##        # emi.add("co2", strength_co2, True, x, y, z, sigma_x, sigma_y, sigma_z)

if sw_chemistry:
    emi.add("nh3", strength_nh3, True, x, y, z, sigma_x, sigma_y, sigma_z)

if sw_plume_rise:
    emi.add("thl", strength_T, False, x, y, z, sigma_x, sigma_y, sigma_z)
    emi.add("qt", strength_q, False, x, y, z, sigma_x, sigma_y, sigma_z)

##############################################################################################
##Create heterogeneous land-surface with three distinct regions (grass_in, forest, grass_out).
##############################################################################################


if (sw_land_surface):
    # Define block sizes for x and y directions
    blocksize_i = 4  # Size of block in x direction (in grid points)
    blocksize_j = 4  # Size of block in y direction (in grid points)
    
    # Calculate size of each block in total grid points
    block_points = blocksize_i * blocksize_j  # Total grid points in one block
    
    # Calculate total number of possible regions in each direction
    region_sizex = itot // blocksize_i  # Number of possible regions in x direction
    region_sizey = jtot // blocksize_j  # Number of possible regions in y direction
    
    # Print information about block structure
    print(f"\nBlock structure:")
    print(f"Each block is {blocksize_i}×{blocksize_j} grid points")
    print(f"Maximum possible regions: {region_sizex} (x) × {region_sizey} (y)")

    # EDITED: Define boundaries for the three regions (in grid points)
    # Each grid point represents 58m (11136/192 = 58)
    # Adjusting boundaries to align with block size
    grass_in_end = (10 // blocksize_i) * blocksize_i      # Align to block boundary (near 1508m)
    forest_end = (21 // blocksize_i) * blocksize_i       # Align to block boundary (near 6496m)
    # grass_out continues to end (192 * 58 = 11136m)
    
    # Print information about domain setup
    print(f"\nDomain setup (aligned to block boundaries):")
    print(f"Grass_in region: 0 to {grass_in_end*58}m (0 to {grass_in_end} points)")
    print(f"Forest region: {grass_in_end*58}m to {forest_end*58}m ({grass_in_end} to {forest_end} points)")
    print(f"Grass_out region: {forest_end*58}m to {itot*58}m ({forest_end} to {itot} points)")
    
    ##########################################################
    # Create and save mask files according to MicroHH format
    ##########################################################
    
    # Initialize masks for three regions
    mask_grass_in = np.zeros((jtot, itot), dtype=float_type)
    mask_forest = np.zeros((jtot, itot), dtype=float_type)
    mask_grass_out = np.zeros((jtot, itot), dtype=float_type)
    
    # Set the regions according to specified boundaries using block-wise assignment
    for j in range(0, jtot, blocksize_j):
        j_slice = slice(j, j + blocksize_j)
        
        # Assign grass_in region blocks
        for i in range(0, grass_in_end, blocksize_i):
            i_slice = slice(i, i + blocksize_i)
            mask_grass_in[j_slice, i_slice] = 1
        
        # Assign forest region blocks
        for i in range(grass_in_end, forest_end, blocksize_i):
            i_slice = slice(i, i + blocksize_i)
            mask_forest[j_slice, i_slice] = 1
        
        # Assign grass_out region blocks
        for i in range(forest_end, itot, blocksize_i):
            i_slice = slice(i, i + blocksize_i)
            mask_grass_out[j_slice, i_slice] = 1
    
    # Save binary mask files according to MicroHH format
    mask_grass_in.tofile("grass_in.0000000")
    mask_forest.tofile("forest.0000000")
    mask_grass_out.tofile("grass_out.0000000")
    
    ##########################################################
    # Initialize and setup Land Surface Model
    ##########################################################
    
    # Initialize LSM
    ls = LSM_input(itot, jtot, 4, sw_water=True, TF=float_type, debug=True, exclude_fields=["z0m", "z0h"])
    
    # Create boolean masks for LSM setup (aligned with block boundaries)
    mask_forest_bool = np.zeros((jtot, itot), dtype=bool)
    mask_forest_bool[:, grass_in_end:forest_end] = True
    
    def set_value(variable, forest, grass):
        ls[variable][mask_forest_bool] = forest      # Forest values
        ls[variable][~mask_forest_bool] = grass      # Grass values (both in and out)
    
    # Set surface properties for forest and grass regions
    set_value("c_veg", forest=1.0, grass=1.0)      # Vegetation fraction
    set_value("lai", forest=6.0, grass=0.1)        # Leaf Area Index
    set_value("water_mask", forest=0, grass=0)     # No water surfaces

    # c_veg (Vegetation Cover Fraction)

    
    ##########################################################
    # Set constant properties for all surfaces
    ##########################################################
    
    ### Surface parameters
    ##ls["gD"][:,:] = 0                 # Vegetation water stress parameter
    ##ls["rs_veg_min"][:,:] = 100       # Minimum vegetation surface resistance
    ##ls["rs_soil_min"][:,:] = 50       # Minimum soil surface resistance
    ##ls["lambda_stable"][:,:] = 10     # Skin conductivity for stable conditions (W/m²/K) 
    ##ls["lambda_unstable"][:,:] = 10   # Skin conductivity for unstable conditions (W/m²/K) 
    ##ls["cs_veg"][:,:] = 0             # Vegetation heat capacity
    ##ls["t_bot_water"][:,:] = 295     # Bottom water temperature (K)
    ##
    ### Soil properties (4 layers)
    ##ls["t_soil"][:,:,:] = 290        # Soil temperature
    ##ls["theta_soil"][:,:,:] = 0.2    # Volumetric soil moisture content (m³/m³)
    ##ls["index_soil"][:,:,:] = 0      # Soil type index
    ##ls["root_frac"][:,:,:] = 0.25    # Root fraction distribution in each soil layer
    

    # Surface parameters
    ls["gD"][:,:] = 0.6               # Vegetation water stress parameter
    ls["rs_veg_min"][:,:] =300        # Minimum vegetation surface resistance
    ls["rs_soil_min"][:,:] = 200      # Minimum soil surface resistance
    ls["lambda_stable"][:,:] = 20     # Skin conductivity for stable conditions (W/m²/K) 
    ls["lambda_unstable"][:,:] = 20   # Skin conductivity for unstable conditions (W/m²/K) 
    ls["cs_veg"][:,:] = 0.5           # Vegetation heat capacity
    ls["t_bot_water"][:,:] = 298     # Bottom water temperature (K)
    
    # Soil properties (4 layers)
    ls["t_soil"][:,:,:] = 295        # Soil temperature
    ls["theta_soil"][:,:,:] = 0.1    # Volumetric soil moisture content (m³/m³)
    ls["index_soil"][:,:,:] = 0      # Soil type index
    ls["root_frac"][:,:,:] = 0.15    # Root fraction distribution in each soil layer

    # Check if all values are set
    ls.check()
    
    # Save LSM setup
    ls.save_binaries(allow_overwrite=True)
    ls.save_netcdf("lsm_input.nc", allow_overwrite=True)

    ###########################
    #Add settings to .ini file.
    ###########################
    
    # Sets grid parameters:
    ini["grid"]["itot"] = itot
    ini["grid"]["jtot"] = jtot
    ini["grid"]["ktot"] = ktot
    
    ini["grid"]["xsize"] = xsize
    ini["grid"]["ysize"] = ysize
    ini["grid"]["zsize"] = zsize
    
    
    # Add statistics settings for masks
    ini["stats"]["xymasklist"] = "grass_in,forest,grass_out"
    

    # Handles scalar variables:
    scalars = list(species.keys())
    ini["advec"]["fluxlimit_list"] = scalars
    ini["limiter"]["limitlist"] = scalars
    ini["fields"]["slist"] = scalars
    ini["boundary"]["scalar_outflow"] = scalars
    
    ini["time"]["endtime"] = (end_hour - start_hour) * 3600
    
    # Configures sources and chemistry:
    ini["source"]["sourcelist"] = emi.source_list
    
    ini["chemistry"]["swchemistry"] = sw_chemistry
    
    crosslist = ["thl", "qt", "u", "v", "w", "thl_fluxbot", "qt_fluxbot","flux_nh3"]
    
    if (sw_chemistry):
        # Add chemicial species and their vertical integrals.
        crosslist += list(species.keys())
        crosslist += [f"{x}_path" for x in species.keys()]
        crosslist += [f"vd{x}" for x in deposition_species]
    
    if (sw_land_surface and sw_chemistry):
        # Add deposition for each land-surface tile.
        for s in deposition_species:
            for t in ["soil", "wet", "veg"]:
                crosslist.append(f"vd{s}_{t}")

    ## Configures boundary conditions based on land-surface flag:
    # With land-surface:
    if (sw_land_surface):
        ini["boundary"]["swboundary"] = "surface_lsm"
        ini["boundary"]["sbcbot"] = "flux"
        ini["boundary"]["sbot"] = "0"
        ini["boundary"]["thl"] = "dirichlet"
        ini["boundary"]["qt"] = "dirichlet"
        ini["boundary"]["swtimedep"] = False
        ini["boundary"]["timedeplist"] = "empty"
    
        ini["radiation"]["swradiation"] = "prescribed"
    
    # Without land-surface (Radiation source is off):
    else:
        ini["boundary"]["swboundary"] = "surface"
        ini["boundary"]["sbcbot"] = "flux"
        ini["boundary"]["swtimedep"] = True
        ini["boundary"]["timedeplist"] = ["thl_sbot", "qt_sbot"]
    
        ini["radiation"]["swradiation"] = False
    
    
    ## Sets deposition settings:
    #  if BOTH chemistry AND land surface are enabled, turns ON deposition in the model!
    if (sw_chemistry and sw_land_surface):
        ini["deposition"]["swdeposition"] = True
    
    # If either one or both are disabled, turns OFF deposition in the model!
    else:
        ini["deposition"]["swdeposition"] = False
    
    # Configures cross-section output and source parameters:
    ini["cross"]["crosslist"] = crosslist
    ini["cross"]["xz"] = ysize/2
    
    # Adds emission source locations and parameters:
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
    
    #  saves the configuration:
    ini.save("plume_chem.ini", allow_overwrite=True)

###########################################
## "gD"
##
## According to the IFS documentation, gD is directly related to vegetation water stress 
## through the stomatal conductance equation. 
## 
## The stomatal conductance (inverse of resistance) in the IFS model is calculated using 
## a coupling factor f that depends on atmospheric moisture deficit Ds:
##     f = f0(1 - Ds/Dmax) + fmin(Ds/Dmax)
## where:
##     Ds    - Atmospheric moisture deficit (vapor pressure deficit)
##     Dmax  - Maximum moisture deficit parameter
##     f0    - Value of f when Ds = 0 (unstressed condition)
##     fmin  - Minimum value of f
## 
## The parameter gD, mentioned in Table 8.1 of the IFS documentation, determines the sensitivity 
## of the stomatal response to vapor pressure deficit. It"s used in calculating the coupling factor f 
## through the relationship:
##     1/f3(Da) = exp(-gD * Da)
## where:
##     f3    - Stress function
##     Da    - Atmospheric humidity deficit
##     gD    - Vapor pressure deficit sensitivity parameter
## 
## From Table 8.1:
##     - High vegetation (trees) typically has gD ≈ 0.03
##     - Low vegetation (grass, crops) has gD = 0
## 
## This indicates that trees are more sensitive to atmospheric drying than grasses.
## 
## In conclusion, gD directly affects how plants respond to atmospheric water stress by 
## modifying their stomatal conductance. A higher gD implies that vegetation reduces its 
## conductance more strongly as the air becomes drier, making it more conservative with 
## water use under stress conditions.
#
###########################################
## Vapor Pressure Deficit (VPD)
## VPD is a crucial concept in plant physiology and agriculture. 
## It measures the difference between the amount of moisture in the air 
## and the amount of moisture the air can hold when saturated. 
## This metric plays a vital role in understanding and managing plant growth, 
## transpiration, and overall health.
##
## Understanding VPD:
## VPD is calculated by subtracting the actual vapor pressure in the air 
## from the saturation vapor pressure at a given temperature. 
## It is typically expressed in units of pressure, such as kilopascals (kPa). 
## The formula for VPD can be simplified as:
## VPD = Vapor Pressure (saturation) – Vapor Pressure (air)
##
## Importance in Plant Growth:
## VPD directly affects plant transpiration rates. 
## It makes VPD a proactive tool for plant empowerment. 
## - Transpiration Control: VPD drives the transpiration process in plants. 
##   When VPD is high (dry air), transpiration rates increase. 
##   When VPD is low, transpiration decreases or stops.
## - Nutrient Uptake: As water evaporates from leaves, 
##   it creates a "pull" effect that draws water and nutrients through the roots. 
##   This process is essential for plant metabolism and development.
## - Growth Regulation: Optimal VPD levels promote healthy growth. 
##   Extreme levels can impede development. 
##   Low VPD can slow metabolism and increase disease susceptibility. 
##   Very high VPD can lead to water stress and impede CO2 intake.
##
## Optimal VPD Range:
## For most plants, the ideal VPD range in a greenhouse setting 
## is between 0.45 kPa to 1.25 kPa. 
## The optimal point is around 0.85 kPa. 
## Specific plants may have different requirements. 
## The ideal VPD can vary depending on the growth stage 
## and environmental conditions.
#
###########################################
##
## "rs_veg_min"
##
## The minimum canopy resistance (rs_veg_min) appears in both the A-gs model and the photosynthesis schemes. Here"s the detailed explanation:
## The total canopy resistance rc (reciprocal of conductance) is calculated as:
## rc = rs,min/(LAI * f1(Rs) * f2(θ̄) * f3(Da))
## where:
## rs,min is the minimum stomatal resistance (rs_veg_min in the code)
## LAI is the leaf area index
## f1(Rs) is a hyperbolic function of downward shortwave radiation:
## 1/f1(Rs) = min[1, (bRs + c)/(a(bRs + 1))]
## where a = 0.81, b = 0.004 W⁻¹m², c = 0.05
## f2(θ̄) is a soil moisture stress function:
## For θ < θpwp: f2 = 0
## For θpwp ≤ θ ≤ θcap: f2 = (θ - θpwp)/(θcap - θpwp)
## For θ > θcap: f2 = 1
## where θpwp is permanent wilting point and θcap is field capacity
## f3(Da) is the vapor pressure deficit function discussed earlier:
## 1/f3(Da) = exp(-gD*Da)
## 
## Typical values for rs_veg_min from Table 8.1:
## Crops/mixed farming: 100 s/m
## Short grass: 100 s/m
## Evergreen needleleaf trees: 250 s/m
## Deciduous broadleaf trees: 175 s/m
## Tall grass: 100 s/m
## Tundra: 80 s/m
## Irrigated crops: 180 s/m
## Bogs and marshes: 240 s/m
## 
## This parameter represents the minimum resistance to water vapor transfer when:
## Light conditions are optimal (f1 ≈ 1)
## Soil moisture is adequate (f2 ≈ 1)
## Vapor pressure deficit is low (f3 ≈ 1)
###########################################
## For Dutch Grassland on a typical sunny day:
###########################################
## rs_veg_min = 100 s/m
## Reasoning:
## Dutch grasslands are typically well-watered due to regular rainfall.
## The temperate climate means moderate vapor pressure deficits.
## Grass species in the Netherlands are generally lush and productive.
## Regular precipitation and moderate temperatures allow for lower resistance values.
## Most Dutch grasslands are managed for agriculture (dairy farming), so they"re usually in good condition.
## 
###########################################
## For Dutch coniferous forest on a typical sunny day:
###########################################
## rs_veg_min = 250 s/m
## 
## Reasoning from the IFS documentation and coniferous forest characteristics:
##
## Species-specific factors:
## - Dutch coniferous forests are mainly composed of Scots pine (Pinus sylvestris) and Douglas fir.
## - These species have evolved more conservative water use strategies.
## - Needle structure is designed to minimize water loss.
## - Table 8.1 shows evergreen needleleaf trees have the highest rs_veg_min of 250 s/m.
##
## Climate considerations:
## - Despite the Netherlands" maritime climate with regular rainfall, 
##   coniferous trees maintain a conservative water use strategy even in well-watered conditions.
## - Their evergreen nature means they need a year-round water management strategy.
## - Higher resistance helps protect against winter water stress when soil might be frozen.
##
## This value of 250 s/m is appropriate for:
## - A typical sunny summer day in the Netherlands.
## - Normal soil moisture conditions.
## - An established coniferous forest.
## - Non-drought conditions.
##
## The value might need to be increased:
## - During prolonged dry periods.
## - Under drought conditions.
## - During winter when soil water might be less accessible.
## - In sandy soil areas (common in Dutch coniferous forests).
###########################################
## For Dutch Bare Soil Areas (based on IFS documentation):
###########################################
##
## rs_soil_min (Minimum soil resistance):
## - Recommended value: 50 s/m
##
## Reasoning:
## - Dutch soils are typically well-watered.
## - Common soil types in the Netherlands include:
##   - Sandy soils (coastal areas and eastern regions).
##   - Clay soils (river areas and polders).
##   - Peat soils (western and northern regions).
## - The maritime climate keeps soil moisture relatively high.
## - Regular precipitation reduces the need for high resistance values.
##
## Adjustments to rs_soil_min:
## - Higher values (100-150 s/m) for:
##   - Sandy soils with low water retention.
##   - During dry periods.
##   - Urban or compacted soils.
## - Lower values (30-40 s/m) for:
##   - Clay or peat soils with high water retention.
##   - After rainfall.
##   - Areas with a high water table.
##
## Other Relevant Parameters for Dutch Bare Soil:
## - gD = N/A (no vegetation).
## - lambda_stable = 15 W/m²/K (thermal conductivity for stable conditions).
## - lambda_unstable = 15 W/m²/K (thermal conductivity for unstable conditions).
## - cs_veg = 0 (no vegetation).
## - theta_soil ≈ 0.2-0.3 m³/m³ (typical soil moisture content for Dutch soils under normal conditions).
## - index_soil depends on location:
##   - 0 (Coarse) for sandy areas.
##   - 2-3 (Medium-fine to Fine) for clay areas.
##   - 5 (Organic) for peat areas.


