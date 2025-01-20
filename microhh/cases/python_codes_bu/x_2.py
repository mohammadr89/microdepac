import netCDF4 as nc    # Library for reading NetCDF files
import numpy as np      # Library for numerical operations

# Define domain specifications
domain_size = 3200  # Physical size of domain in meters (0 to 3200m in x, y, and z)
npoints = 32        # Number of grid points in each direction (32x32x32 grid)

# Calculate grid spacing
# Since domain is 3200m with 32 points, each grid cell is 3200/32 = 100m
dy = domain_size / npoints  # Grid spacing in y-direction (meters)
dz = domain_size / npoints  # Grid spacing in z-direction (meters)

# Open NetCDF files in read mode ('r')
# These files contain our simulation data
default_nc = nc.Dataset('plume_chem.default.0010800.nc', 'r')  # Contains density data
u_yz_nc = nc.Dataset('u.yz.nc', 'r')                          # Contains velocity data

# Get density profile (rhoref) from the thermo group in default_nc
# rhoref varies with height (z) only
rhoref = default_nc.groups['thermo'].variables['rhoref'][:]  # [:] loads all data

# Set target x-location for flux calculation
x_target = 2200  # We want to calculate flux at x = 2200m

# Find the index in the x-grid that corresponds to x=1000m
# xh[:] reads all x-coordinates, np.where finds where x equals target value
x_index = np.where(u_yz_nc.variables['xh'][:] == x_target)[0][0]

# Extract velocity data (u) at x=1000m
# u_data dimensions: [time, z, y, x]
# [0, :, :, x_index] means:
#   - time = 0 (first timestep)
#   - : means all points in z direction
#   - : means all points in y direction
#   - x_index is our target x-location
u_data = u_yz_nc.variables['u'][0, :, :, x_index]

# Prepare density for multiplication with velocity
# rhoref is 1D (varies with z only)
# We need to make it 2D to multiply with u_data (which is z,y)
# [:, np.newaxis] adds a new axis, making rhoref a 2D array
rho_2d = rhoref[:, np.newaxis]  # Shape becomes (32, 1)

# Multiply density and velocity
# Through broadcasting, rho_2d (32,1) * u_data (32,32) becomes (32,32)
flux_density = rho_2d * u_data

# Perform double integration using simple rectangular method
# np.sum adds up all values in the array
# Multiply by dy*dz to account for grid cell area
mass_flux = np.sum(flux_density) * dy * dz

# Print results with detailed information
print(f"Mass flux calculation for x = {x_target}m:")
print(f"X-index in grid: {x_index}")
print(f"Available x coordinates: {u_yz_nc.variables['xh'][:]}")
print(f"Mass flux value = {mass_flux:.2f} kg/s")
print(f"\nGrid information:")
print(f"Domain size: {domain_size}m x {domain_size}m x {domain_size}m")
print(f"Grid points: {npoints} x {npoints} x {npoints}")
print(f"dy = {dy}m")
print(f"dz = {dz}m")
print(f"Density range: {rhoref.min():.3f} to {rhoref.max():.3f} kg/m³")

# Close NetCDF files to free up system resources
default_nc.close()
u_yz_nc.close()
