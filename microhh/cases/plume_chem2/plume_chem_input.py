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
float_type = 'f8'

xsize = 3200
ysize = 3200
zsize = 3200

itot = 32
jtot = 32
ktot = 32

start_hour = 6
end_hour = 18

# Enable resolved plume rise:
sw_plume_rise = False

# Enable non-linear KPP chemistry:
sw_chemistry = True

# Enable land-surface model and more detailled deposition.
sw_land_surface = True


"""
Read base .ini file for case settings.
"""
ini = mht.Read_namelist('plume_chem.ini.base')


"""
Create case input.
"""
if (sw_chemistry):

    # Read TUV output table.
    # There is a total of 24 hour available, generated for the
    # Jaenschwalde power plant on 23/05/2022.

    # columns for time, solar zenith angle, and various photolysis rates:
    columns = ['time', 'sza', 'jo31d', 'jh2o2', 'jno2', 'jno3', 'jn2o5', 'jch2or', 'jch2om', 'jch3o2h']
    tuv = pd.read_table(
            'plume_chem_tuv_output.txt',
            sep='\\s+',
            skiprows=12,
            skipfooter=1,
            engine='python',
            names=columns,
            index_col='time')

    # NOTE: `.loc` is value based, not on index.
    # TUV (Tropospheric Ultraviolet and Visible)
    tuv = tuv.loc[start_hour:end_hour]

    # Convert to seconds, and subtract starting time.
    tuv.index *= 3600
    tuv.index -= tuv.index.values[0]

    # Emissions (?)
    # Sets up emissions arrays for NO and ISOP (isoprene)
    emi_no = np.zeros(tuv.index.size)
    emi_isop = np.zeros(tuv.index.size)

    # Concentrations, for now constant with height.
    species = {
        'nh3': 1.0}

    deposition_species = ['nh3']
else:
    species = {}


"""
Define vertical grid and initial thl/qt/co2 profiles.
"""
dz = zsize / ktot
z = np.arange(0.5*dz, zsize, dz)

def profile(zi, v_bulk, dv, gamma_v, clip_at_zero=False):
    """
    Create well mixed profile with jump and constant lapse rate above.
    """

    k_zi = np.abs(z - zi).argmin()

    profile = np.zeros(ktot)
    profile[:k_zi] = v_bulk
    profile[k_zi:] = v_bulk + dv + gamma_v * (z[k_zi:] - zi)

    # Option to clip negative values:
    if clip_at_zero:
        profile[profile < 0] = 0.

    return profile


# Vertical profiles.
thl = profile(zi=1000, v_bulk=290,   dv=2,     gamma_v=0.006)
qt  = profile(zi=1000, v_bulk=10e-3, dv=-2e-3, gamma_v=-0.003e-3, clip_at_zero=True)
u   = np.ones(ktot) * 5
co2 = np.zeros(ktot)

# Surface fluxes.
t0 = start_hour*3600
t1 = end_hour*3600
td = 12*3600
time = np.linspace(t0, t1, 32)


wthl = 0.15
wqt  = 8e-5 

# wthl = 0.15 * np.sin(np.pi * (time-t0) / td)
# wqt  = 8e-5 * np.sin(np.pi * (time-t0) / td)

# ##########################################################
# # Surface radiation (only used with land-surface enabled).
# ##########################################################
# sw_flux_dn = 600 * np.sin(np.pi * (time-t0) / td)
# sw_flux_dn[sw_flux_dn < 0] = 0
# 
# sw_flux_up = 0.2 * sw_flux_dn
# 
# lw_flux_dn = np.ones_like(sw_flux_dn) * 340
# lw_flux_up = np.ones_like(sw_flux_dn) * 400
# 
# ##########################################################

"""
pl.figure(figsize=(10,5))

pl.subplot(131)
pl.plot(time/3600, wthl)
pl.xlabel('time (h)')
pl.ylabel('w`thl` (K m s-1)')

pl.subplot(132)
pl.plot(time/3600, wqt*1000)
pl.xlabel('time (h)')
pl.ylabel('w`qt` (g kg-1 m s-1)')

pl.subplot(133)
pl.plot(time/3600, sw_flux_dn, label='sw_flux_dn')
pl.plot(time/3600, sw_flux_up, label='sw_flux_up')
pl.plot(time/3600, lw_flux_dn, label='lw_flux_dn')
pl.plot(time/3600, lw_flux_up, label='lw_flux_up')
pl.xlabel('time (h)')
pl.ylabel('sw_flux_dn` (W m-2)')
pl.legend()

pl.tight_layout()
"""

"""
Write input NetCDF file.
"""

# Defines helper function add_nc_var to add variables to NetCDF file:
def add_nc_var(name, dims, nc, data):
    if dims is None:
        var = nc.createVariable(name, np.float64)
    else:
        var = nc.createVariable(name, np.float64, dims)
    var[:] = data

# Creates 'plume_chem_input.nc' file
nc_file = nc.Dataset('plume_chem_input.nc', mode='w', datamodel='NETCDF4', clobber=True)

###############################
# Sets up dimensions and groups
###############################

# 'z' dimension for vertical levels:
nc_file.createDimension('z', ktot)
add_nc_var('z', ('z'), nc_file, z)

# Atmospheric input.
# ('init' group for initial atmospheric conditions)
nc_group_init = nc_file.createGroup('init');

add_nc_var('u', ('z'), nc_group_init, u)
add_nc_var('thl', ('z'), nc_group_init, thl)
add_nc_var('qt', ('z'), nc_group_init, qt)
#add_nc_var('co2', ('z'), nc_group_init, co2)
#add_nc_var('co2_inflow', ('z'), nc_group_init, co2)

#('timedep' group for time-dependent surface conditions)
nc_tdep = nc_file.createGroup('timedep');
nc_tdep.createDimension("time_surface", time.size)

add_nc_var('time_surface', ('time_surface'), nc_tdep, time-time[0])
add_nc_var('thl_sbot', ('time_surface'), nc_tdep, wthl)
add_nc_var('qt_sbot', ('time_surface'), nc_tdep, wqt)

if (sw_chemistry):
    # Chemistry input.
    # 'timedep_chem' group (if chemistry enabled)
    nc_chem = nc_file.createGroup('timedep_chem');
    nc_chem.createDimension("time_chem", tuv.index.size)

    add_nc_var("time_chem", ('time_chem'), nc_chem, tuv.index)
    add_nc_var("jo31d", ('time_chem'), nc_chem, tuv.jo31d)
    add_nc_var("jh2o2", ('time_chem'), nc_chem, tuv.jh2o2)
    add_nc_var("jno2", ('time_chem'), nc_chem, tuv.jno2)
    add_nc_var("jno3", ('time_chem'), nc_chem, tuv.jno3)
    add_nc_var("jn2o5", ('time_chem'), nc_chem, tuv.jn2o5)
    add_nc_var("jch2or", ('time_chem'), nc_chem, tuv.jch2or)
    add_nc_var("jch2om", ('time_chem'), nc_chem, tuv.jch2om)
    add_nc_var("jch3o2h", ('time_chem'), nc_chem, tuv.jch3o2h)
    add_nc_var("emi_isop", ('time_chem'), nc_chem, emi_isop)
    add_nc_var("emi_no", ('time_chem'), nc_chem, emi_no)

    for name, value in species.items():
        profile = np.ones(ktot, dtype=np.float64)*value
        add_nc_var(name, ('z'), nc_group_init, profile)
        add_nc_var('{}_inflow'.format(name), ('z'), nc_group_init, profile)


    # Add flux_nh3 variable
    add_nc_var('flux_nh3', ('z'), nc_group_init, np.zeros(ktot, dtype=np.float64))

if (sw_land_surface):
    # 'soil' group (if land surface enabled)
    nc_soil = nc_file.createGroup('soil')
    nc_soil.createDimension('z', 4)
    add_nc_var('z', ('z'), nc_soil, np.array([-1.945, -0.64, -0.175, -0.035]))

    add_nc_var('theta_soil', ('z'), nc_soil, np.array([0.34, 0.25, 0.21, 0.18]))
    add_nc_var('t_soil', ('z'), nc_soil, np.array([282, 287, 290, 286]))
    add_nc_var('index_soil', ('z'), nc_soil, np.ones(4) * 2)
    add_nc_var('root_frac', ('z'), nc_soil, np.array([0.05, 0.3, 0.4, 0.25]))

#    # Add idealized (prescribed) radiation.
#    add_nc_var('sw_flux_dn', ('time_surface'), nc_tdep, sw_flux_dn)
#    add_nc_var('sw_flux_up', ('time_surface'), nc_tdep, sw_flux_up)
#    add_nc_var('lw_flux_dn', ('time_surface'), nc_tdep, lw_flux_dn)
#    add_nc_var('lw_flux_up', ('time_surface'), nc_tdep, lw_flux_up)

nc_file.close()


"""
Define emissions.
"""
# Coordinates of central cooling tower (m):
x0 = 200.0
y0 = ysize/2.

# # Std-dev of plume widths:
# sigma_x = 25
# sigma_y = 25
sigma_x = 9
sigma_y = 9
sigma_z = 9
z0 = 9.0

# # Handles plume rise conditions:
# # If enabled: uses tower height (120m)
# # If disabled: uses fitted heights from CSV profiles
# 
# if sw_plume_rise:
#     # Emissions from tower height.
#     z0 = 120
#     sigma_z = 25
# else:
#     # The heights and sigma are from Dominik's .csv profiles, curve fitted with Python.
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

# Add emission from individual towers.
# NOTE: Jaenschwalde has nine cooling towers, in three groups.
#       In reality, only two groups are active at the same time.
#       Here, we distribute the emissions over all nine towers......
# (Creates emissions for nine cooling towers in three groups)
# for j in range(-1,2):
##    for i in range(-1,2):
##        x = x0 + i*dx + j*ddx
##        y = y0 + j*dy
##        z = z0


x = x0
y = y0
z = z0

        #emi.add('co2', strength_co2, True, x, y, z, sigma_x, sigma_y, sigma_z)

if (sw_chemistry):
    emi.add('nh3', strength_nh3, True, x, y, z, sigma_x, sigma_y, sigma_z)

if sw_plume_rise:
    emi.add('thl', strength_T, False, x, y, z, sigma_x, sigma_y, sigma_z)
    emi.add('qt',  strength_q, False, x, y, z, sigma_x, sigma_y, sigma_z)

"""
Create heterogeneous land-surface, here with simple block pattern to look at deposition differences.
"""

if (sw_land_surface):
    # Define block sizes for x and y directions
    blocksize_i = 1  # Size of block in x direction (in grid points)
    blocksize_j = 1  # Size of block in y direction (in grid points)
    
    # Calculate size of each block in total grid points
    block_points = blocksize_i * blocksize_j  # Total grid points in one block
    
    # Calculate total number of possible regions in each direction
    region_sizex = itot // blocksize_i  # Number of possible regions in x direction 
    region_sizey = jtot // blocksize_j  # Number of possible regions in y direction
    
    # Print information about domain setup
    print(f"Each block is {blocksize_i}×{blocksize_j} grid points")
    print(f"Maximum possible splits in x-direction: {region_sizex}")
    print(f"Maximum possible splits in y-direction: {region_sizey}")
    
    # Define number of splits desired in each direction
    n_splits_x = 2  # Number of different surface type regions in x direction
    n_splits_y = 1  # Keep y direction as one region for now
    
    # Check if requested splits are valid
    if n_splits_x > region_sizex:
       raise ValueError(f"Number of splits ({n_splits_x}) cannot be larger than total regions in x ({region_sizex})")
    
    ##########################################################
    # Calculate region sizes and positions
    ##########################################################
    
    # Calculate how many regions belong to each split
    regions_per_split_x = region_sizex // n_splits_x  # Width of each split in x
    regions_per_split_y = region_sizey // n_splits_y  # Width of each split in y
    
    # Handle both odd and even numbers of splits
    if n_splits_x % 2 == 0:  # Even number of splits
        # For even numbers like 4, you can choose position 2 or 3
        # middle_split = n_splits_x // 2 - 1  # Left-middle (position 2 in case of 4 splits)
        middle_split = n_splits_x // 2      # Right-middle (position 3 in case of 4 splits)
    else:  # Odd number of splits
        middle_split = (n_splits_x - 1) // 2  # True middle position
    
    start_forest = middle_split * regions_per_split_x * blocksize_i
    end_forest = (middle_split + 1) * regions_per_split_x * blocksize_i
    
    print(f"Total splits: {n_splits_x}")
    print(f"Forest split position: {middle_split + 1} out of {n_splits_x}")
    print(f"Each split width: {regions_per_split_x} grid points")
    
    ##########################################################
    # Create and save mask files according to MicroHH format
    ##########################################################
    
    # Initialize masks for forest and grass
    mask_forest = np.zeros((jtot, itot), dtype=float_type)
    mask_grass = np.zeros((jtot, itot), dtype=float_type)
    
    # Set the middle section as forest, rest as grass
    mask_forest[:, start_forest:end_forest] = 1  # Sets middle section to forest
    mask_grass[:, :start_forest] = 1            # Sets left sections to grass
    mask_grass[:, end_forest:] = 1              # Sets right sections to grass
    
    # Save binary mask files according to MicroHH format
    mask_forest.tofile('forest.0000000')
    mask_grass.tofile('grass.0000000')
    
    print(f"Created masks with dimensions: {jtot}x{itot}")
    print(f"Forest region: columns {start_forest} to {end_forest}")
    
    ##########################################################
    # Initialize and setup Land Surface Model
    ##########################################################
    
    # Initialize LSM
    ls = LSM_input(itot, jtot, 4, sw_water=True, TF=float_type, debug=True, exclude_fields=['z0m', 'z0h'])
    
    # Create boolean mask for LSM setup
    mask = np.zeros((jtot, itot), dtype=bool)
    mask[:, start_forest:end_forest] = True  # True for forest section
    
    def set_value(variable, forest, grass):
        ls[variable][ mask] = forest    # Forest values
        ls[variable][~mask] = grass     # Grass values
    
    # Set surface properties for forest and grass regions
    set_value('c_veg', forest=1.0, grass=0.9)      # Vegetation fraction
    set_value('lai', forest=5.0, grass=2.0)        # Leaf Area Index
    set_value('water_mask', forest=0, grass=0)     # No water surfaces
    
    ##########################################################
    # Set constant properties for all surfaces
    ##########################################################
    
    # Surface parameters
    ls['gD'][:,:] = 0                 # Vegetation water stress parameter
    ls['rs_veg_min'][:,:] = 100       # Minimum vegetation surface resistance
    ls['rs_soil_min'][:,:] = 50       # Minimum soil surface resistance
    ls['lambda_stable'][:,:] = 10     # Stability parameter
    ls['lambda_unstable'][:,:] = 10   # Stability parameter
    ls['cs_veg'][:,:] = 0            # Surface heat capacity
    ls['t_bot_water'][:,:] = 295     # Bottom water temperature
    
    # Soil properties (4 layers)
    ls['t_soil'][:,:,:] = 290        # Soil temperature
    ls['theta_soil'][:,:,:] = 0.2    # Soil moisture
    ls['index_soil'][:,:,:] = 0      # Soil type index
    ls['root_frac'][:,:,:] = 0.25    # Root distribution
    
    # Check if all values are set
    ls.check()
    
    # Save LSM setup
    ls.save_binaries(allow_overwrite=True)
    ls.save_netcdf('lsm_input.nc', allow_overwrite=True)

    """
    Add settings to .ini file.
    """
    
    # Sets grid parameters:
    ini['grid']['itot'] = itot
    ini['grid']['jtot'] = jtot
    ini['grid']['ktot'] = ktot
    
    ini['grid']['xsize'] = xsize
    ini['grid']['ysize'] = ysize
    ini['grid']['zsize'] = zsize
    
    
    # Add statistics settings for masks
    ini['stats']['xymasklist'] = 'forest,grass'
    
    # Handles scalar variables:
    
    # scalars = ['co2'] + list(species.keys())
    scalars = list(species.keys())
    ini['advec']['fluxlimit_list'] = scalars
    ini['limiter']['limitlist'] = scalars
    ini['fields']['slist'] = scalars
    ini['boundary']['scalar_outflow'] = scalars
    
    ini['time']['endtime'] = (end_hour - start_hour) * 3600
    
    
    # Configures sources and chemistry:
    ini['source']['sourcelist'] = emi.source_list
    
    ini['chemistry']['swchemistry'] = sw_chemistry
    
    ## Sets up crosslist (variables to output):
    # If chemistry enabled: adds species and their paths
    # If land-surface and chemistry enabled: adds deposition for each surface type
    crosslist = ['u', 'v', 'w', 'flux_nh3']
    
    if (sw_chemistry):
        # Add chemicial species and their vertical integrals.
        crosslist += list(species.keys())
        crosslist += [f'{x}_path' for x in species.keys()]
        crosslist += [f'vd{x}' for x in deposition_species]
    
    if (sw_land_surface and sw_chemistry):
        # Add deposition for each land-surface tile.
        for s in deposition_species:
            for t in ['soil', 'wet', 'veg']:
                crosslist.append(f'vd{s}_{t}')
    
    
    
    ## Configures boundary conditions based on land-surface flag:
    # With land-surface:
    if (sw_land_surface):
        ini['boundary']['swboundary'] = 'surface_lsm'
        ini['boundary']['sbcbot'] = 'flux'
        ini['boundary']['sbot'] = '0'
        ini['boundary']['thl'] = 'dirichlet'
        ini['boundary']['qt'] = 'dirichlet'
        ini['boundary']['swtimedep'] = False
        ini['boundary']['timedeplist'] = 'empty'
    
        ini['radiation']['swradiation'] = 'prescribed'
    
    # Without land-surface (Radiation source is off):
    else:
        ini['boundary']['swboundary'] = 'surface'
        ini['boundary']['sbcbot'] = 'flux'
        ini['boundary']['swtimedep'] = True
        ini['boundary']['timedeplist'] = ['thl_sbot', 'qt_sbot']
    
        ini['radiation']['swradiation'] = False
    
    
    ## Sets deposition settings:
    #  if BOTH chemistry AND land surface are enabled, turns ON deposition in the model!
    if (sw_chemistry and sw_land_surface):
        ini['deposition']['swdeposition'] = True
    
    # If either one or both are disabled, turns OFF deposition in the model!
    else:
        ini['deposition']['swdeposition'] = False
    
    # Configures cross-section output and source parameters:
    ini['cross']['crosslist'] = crosslist
    ini['cross']['xz'] = ysize/2
    
    # Adds emission source locations and parameters:
    ini['source']['source_x0'] = emi.x0
    ini['source']['source_y0'] = emi.y0
    ini['source']['source_z0'] = emi.z0
    
    ini['source']['sigma_x'] = emi.sigma_x
    ini['source']['sigma_y'] = emi.sigma_y
    ini['source']['sigma_z'] = emi.sigma_z
    
    ini['source']['strength'] = emi.strength
    ini['source']['swvmr'] = emi.sw_vmr
    
    ini['source']['line_x'] = emi.line_x
    ini['source']['line_y'] = emi.line_y
    ini['source']['line_z'] = emi.line_z
    
    #  saves the configuration:
    ini.save('plume_chem.ini', allow_overwrite=True)
