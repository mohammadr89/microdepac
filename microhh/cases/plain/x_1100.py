#
import netCDF4 as nc    
import numpy as np      

# Open files
nc_def = nc.Dataset('plume_chem.default.0010800.nc', 'r')  
nc_u = nc.Dataset('u.yz.nc', 'r')                          
nc_nh3 = nc.Dataset('nh3.yz.nc', 'r')                     

# Grid spacing (3200m/32points)
dy = dz = 100  
flux1100 = []
flux2200 = []

# Calculate the index corresponding to z=1200m
z_max_index = 12  # since 12 * 100m = 1200m

# Calculate fluxes at x=1100m for all time steps
for i in range(180):
    rho = nc_def.groups['thermo'].variables['rho'][i,:z_max_index]
    u_1100 = nc_u.variables['u'][i, :z_max_index, :, 1]
    nh3_1050 = nc_nh3.variables['nh3'][i, :z_max_index, :, 0]
    nh3_1150 = nc_nh3.variables['nh3'][i, :z_max_index, :, 1]
    nh3_1100 = 0.5 * (nh3_1050 + nh3_1150)
    
    flux_1100 = np.sum(rho[:, np.newaxis] * u_1100 * nh3_1100) * dy * dz * 17.031 / 28.
    flux1100.append(flux_1100)

# Calculate average flux at x=1100m
avg_flux_1100 = np.mean(flux1100)

# Calculate fluxes at x=2200m for all time steps
for i in range(180):
    rho = nc_def.groups['thermo'].variables['rho'][i,:z_max_index]
    u_2200 = nc_u.variables['u'][i, :z_max_index, :, 4]
    nh3_2150 = nc_nh3.variables['nh3'][i, :z_max_index, :, 3]
    nh3_2250 = nc_nh3.variables['nh3'][i, :z_max_index, :, 4]
    nh3_2200 = 0.5 * (nh3_2150 + nh3_2250)
    
    flux_2200 = np.sum(rho[:, np.newaxis] * u_2200 * nh3_2200) * dy * dz * 17.031 / 28.
    flux2200.append(flux_2200)

# Calculate average flux at x=2200m
avg_flux_2200 = np.mean(flux2200)

# Calculate the difference between average fluxes
flux_difference = avg_flux_2200 - avg_flux_1100
percent_diff = (flux_difference / avg_flux_1100) * 100

print(f"\nAverage NH3 mass flux at x = 1100m: {avg_flux_1100:.6e} kg[NH3]/s")
print(f"Average NH3 mass flux at x = 2200m: {avg_flux_2200:.6e} kg[NH3]/s")
print(f"\nDifference in average fluxes (2200m - 1100m): {flux_difference:.6e} kg[NH3]/s")
print(f"Percentage difference: {percent_diff:.2f}%")

# Close files
nc_def.close()
nc_u.close()
nc_nh3.close()

























### Get density
##rho = nc_def.groups['thermo'].variables['rho'][0,:]  # time=0
##
### Calculate for x=1100m forest entry
##u_1100 = nc_u.variables['u'][0, :, :, 1]  # x=1100m
##nh3_1050 = nc_nh3.variables['nh3'][0, :, :, 0]  # NH3 at x=1050
##nh3_1250 = nc_nh3.variables['nh3'][0, :, :, 2]  # NH3 at x=1250
##nh3_1100 = nh3_1050 + (1100-1050)/(1250-1050) * (nh3_1250 - nh3_1050)  # interpolate
##flux_1100 = np.sum(rho[:, np.newaxis] * u_1100 * nh3_1100) * dy * dz
##
### Calculate for x=2200m forest exit
##u_2200 = nc_u.variables['u'][0, :, :, 4]  # x=2200m
##nh3_2150 = nc_nh3.variables['nh3'][0, :, :, 3]  # NH3 at x=2150
##nh3_2350 = nc_nh3.variables['nh3'][0, :, :, 5]  # NH3 at x=2350
##nh3_2200 = nh3_2150 + (2200-2150)/(2350-2150) * (nh3_2350 - nh3_2150)  # interpolate
##flux_2200 = np.sum(rho[:, np.newaxis] * u_2200 * nh3_2200) * dy * dz
##
##print(f"NH3 mass flux at x = 1100m: {flux_1100:.6e} kg/s")
##print(f"NH3 mass flux at x = 2200m: {flux_2200:.6e} kg/s")
##
### Close files
##nc_def.close()
##nc_u.close()
##nc_nh3.close()
