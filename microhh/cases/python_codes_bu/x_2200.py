import netCDF4 as nc    # Library for reading NetCDF files
import numpy as np      # Library for numerical operations

# Define domain specifications
domain_size = 3200  # Physical size of domain in meters (0 to 3200m in x, y, and z)
npoints = 32        # Number of grid points in each direction (32x32x32 grid)

# Calculate grid spacing
dy = domain_size / npoints  # Grid spacing in y-direction (meters)
dz = domain_size / npoints  # Grid spacing in z-direction (meters)

# Open NetCDF files in read mode ('r')
default_nc = nc.Dataset('plume_chem.default.0010800.nc', 'r')  
u_yz_nc = nc.Dataset('u.yz.nc', 'r')                          
nh3_yz_nc = nc.Dataset('nh3.yz.nc', 'r')                     

# Print x coordinates for all fields
print("Available x coordinates in each field:")
print(f"U field x coordinates: {u_yz_nc.variables['xh'][:]}")
print(f"NH3 field x coordinates: {nh3_yz_nc.variables['x'][:]}")

# Get density profile
rho = default_nc.groups['thermo'].variables['rho'][:]  

# Target x-location and NH3 interpolation points
x_target = 2200  # Target x location

# Get x coordinates for NH3
nh3_x = nh3_yz_nc.variables['x'][:]

# Find NH3 interpolation points (2150 and 2250)
x1 = 2150  # First NH3 x point
x2 = 2250  # Second NH3 x point
idx1 = np.where(nh3_x == x1)[0][0]
idx2 = np.where(nh3_x == x2)[0][0]

# Get velocity data index
u_x_index = np.where(u_yz_nc.variables['xh'][:] == x_target)[0][0]

print(f"Calculation for x = {x_target}m:")
print(f"U field: using exact value at x = {x_target}m (index {u_x_index})")
print(f"NH3 field: interpolating between:")
print(f"  - x = {x1}m (index {idx1})")
print(f"  - x = {x2}m (index {idx2})")

# Get NH3 data at both x positions
nh3_1 = nh3_yz_nc.variables['nh3'][0, :, :, idx1]  
nh3_2 = nh3_yz_nc.variables['nh3'][0, :, :, idx2]  

# Linear interpolation for NH3 at x=2200
nh3_interp = nh3_1 + (x_target - x1)/(x2 - x1) * (nh3_2 - nh3_1)

# Get velocity data
u_data = u_yz_nc.variables['u'][0, :, :, u_x_index]

# Prepare density for multiplication (make 2D)
rho_2d = rho[:, np.newaxis]  

# Calculate triple product ρ.u.nh3
flux_density = rho_2d * u_data * nh3_interp

# Perform double integration
mass_flux = np.sum(flux_density) * dy * dz

print(f"\nResults:")
print(f"NH3 mass flux at x = {x_target}m: {mass_flux:.6e} kg/s")
print(f"\nGrid information:")
print(f"dy = {dy}m")
print(f"dz = {dz}m")
print(f"Density range: {rho.min():.3f} to {rho.max():.3f} kg/m³")
print(f"NH3 range at x={x1}m: {np.min(nh3_1):.6e} to {np.max(nh3_1):.6e}")
print(f"NH3 range at x={x2}m: {np.min(nh3_2):.6e} to {np.max(nh3_2):.6e}")
print(f"NH3 interpolated range: {np.min(nh3_interp):.6e} to {np.max(nh3_interp):.6e}")
print(f"U range at x={x_target}m: {np.min(u_data):.6e} to {np.max(u_data):.6e}")

# Close NetCDF files
default_nc.close()
u_yz_nc.close()
nh3_yz_nc.close()
