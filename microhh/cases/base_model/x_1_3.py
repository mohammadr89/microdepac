import netCDF4 as nc    
import numpy as np      

# Define domain specifications
domain_size = 3200  
npoints = 32        
dy = domain_size / npoints  
dz = domain_size / npoints  

# Open NetCDF files
default_nc = nc.Dataset('plume_chem.default.0000000.nc', 'r')  
u_yz_nc = nc.Dataset('u.yz.nc', 'r')                          
nh3_yz_nc = nc.Dataset('nh3.yz.nc', 'r')                     

# Print x coordinates for all fields
print("Available x coordinates in each field:")
print(f"U field x coordinates: {u_yz_nc.variables['xh'][:]}")
print(f"NH3 field x coordinates: {nh3_yz_nc.variables['x'][:]}")
print(f"Note: Density (rho) varies only with z\n")

# Get density profile
rhoref = default_nc.groups['thermo'].variables['rhoref'][:]  

def calculate_flux(x_target, nh3_x):
    """Calculate mass flux for given x location and corresponding NH3 x value"""
    
    # Get indices
    u_x_index = np.where(u_yz_nc.variables['xh'][:] == x_target)[0][0]
    nh3_x_index = np.where(nh3_yz_nc.variables['x'][:] == nh3_x)[0][0]
    
    print(f"\nCalculation for x = {x_target}m:")
    print(f"U field: using exact value at x = {x_target}m (index {u_x_index})")
    print(f"NH3 field: using value at x = {nh3_x}m (index {nh3_x_index})")
    print(f"Density: using z-dependent profile")
    
    # Get data
    u_data = u_yz_nc.variables['u'][0, :, :, u_x_index]
    nh3_data = nh3_yz_nc.variables['nh3'][0, :, :, nh3_x_index]
    rho_2d = rhoref[:, np.newaxis]
    
    # Calculate flux
    flux_density = rho_2d * u_data * nh3_data
    mass_flux = np.sum(flux_density) * dy * dz
    
    # Print results
    print(f"\nResults for x = {x_target}m:")
    print(f"Mass flux: {mass_flux:.6e} kg/s")
    print(f"Using NH3 at x = {nh3_x}m")
    print(f"NH3 range: {np.min(nh3_data):.6e} to {np.max(nh3_data):.6e}")
    print(f"U range: {np.min(u_data):.6e} to {np.max(u_data):.6e}")
    return mass_flux

# Calculate for x=1000 using NH3 at x=1050
flux_1000 = calculate_flux(1000, 1050)

# Calculate for x=1200 using NH3 at x=1150
flux_1200 = calculate_flux(1200, 1150)

print(f"\nSummary:")
print(f"Mass flux at x=1000m: {flux_1000:.6e} kg/s (using NH3 at x=1050m)")
print(f"Mass flux at x=1200m: {flux_1200:.6e} kg/s (using NH3 at x=1150m)")
print(f"\nGrid spacing:")
print(f"dy = {dy}m")
print(f"dz = {dz}m")
print(f"Density range: {rhoref.min():.3f} to {rhoref.max():.3f} kg/m³")

# Close NetCDF files
default_nc.close()
u_yz_nc.close()
nh3_yz_nc.close()
