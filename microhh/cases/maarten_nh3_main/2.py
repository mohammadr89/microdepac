import netCDF4 as nc    
import numpy as np      

# Open files
nc_def = nc.Dataset('plume_chem.default.0010800.nc', 'r')  
nc_u = nc.Dataset('u.yz.nc', 'r')                          
nc_nh3 = nc.Dataset('nh3.yz.nc', 'r')                     

# Grid spacing (3200m/32points)
dy = dz = 100  
uflux1100 = []
uflux2200 = []
flux1100 = []
flux2200 = []

# Calculate the index corresponding to z=1200m (0-1200m = first 13 points since we start at 0)
z_max_index = 32  # 0,100,200,...,1200 = 13 points

for i in range(360):
    # Get density (slice to match dimensions)
    rho = nc_def.groups['thermo'].variables['rho'][i,:z_max_index]  # Only up to 1200m
    
    # Calculate for x=1100m forest entry
    u_1100 = nc_u.variables['u'][i, :z_max_index, :, 1]  
    nh3_1050 = nc_nh3.variables['nh3'][i, :z_max_index, :, 0]  
    nh3_1150 = nc_nh3.variables['nh3'][i, :z_max_index, :, 1]  
    nh3_1100 = 0.5 * (nh3_1050 + nh3_1150)  
    
    # Add np.newaxis to rho after slicing
    flux_1100 = np.sum(rho[:, np.newaxis] * u_1100 * nh3_1100) * dy * dz * 17.031 / 28.
    uflux_1100 = np.sum(rho[:, np.newaxis] * u_1100) * dy * dz * 17.031 / 28.
    
    # Calculate for x=2200m forest exit
    u_2200 = nc_u.variables['u'][i, :z_max_index, :, 4]  
    nh3_2150 = nc_nh3.variables['nh3'][i, :z_max_index, :, 3]  
    nh3_2250 = nc_nh3.variables['nh3'][i, :z_max_index, :, 4]  
    nh3_2200 = 0.5 * (nh3_2150 + nh3_2250)  
    
    flux_2200 = np.sum(rho[:, np.newaxis] * u_2200 * nh3_2200) * dy * dz * 17.031 / 28.
    uflux_2200 = np.sum(rho[:, np.newaxis] * u_2200) * dy * dz * 17.031 / 28.
    
    flux1100.append(flux_1100)
    flux2200.append(flux_2200)
    uflux1100.append(uflux_1100)
    uflux2200.append(uflux_2200)
 
# Calculate statistics
uavg_diff = np.mean(np.array(uflux2200) - np.array(uflux1100))
upercent_diff = np.mean((np.array(uflux2200) - np.array(uflux1100)) / np.array(uflux1100) * 100)

avg_diff = np.mean(np.array(flux2200) - np.array(flux1100))
percent_diff = np.mean((np.array(flux2200) - np.array(flux1100)) / np.array(flux1100) * 100)
print(f"\numass flux at x = 1100m: {uflux_1100:.12e} kg/s")
print(f"nh3 mass flux at x = 2200m: {uflux_2200:.12e} kg/s")
print(f"\nAverage difference in fluxes (2200m - 1100m): {uavg_diff:.12e} kg/s")
print(f"Average percentage difference: {upercent_diff:.12e}%")


print(f"\nnh3 mass flux at x = 1100m: {flux_1100:.6e} kg[nh3]/s")
print(f"nh3 mass flux at x = 2200m: {flux_2200:.6e} kg[nh3]/s")
print(f"\nAverage difference in fluxes (2200m - 1100m): {avg_diff:.6e} g[NH3]/s")
print(f"Average percentage difference: {percent_diff:.2f}%")

# Close files
nc_def.close()
nc_u.close()
nc_nh3.close()
