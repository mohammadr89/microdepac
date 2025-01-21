import netCDF4 as nc    
import numpy as np      
import pylab as plt
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

for i in range(180):

	# Get density
	rho = nc_def.groups['thermo'].variables['rho'][i,:z_max_index]  # time=0, rho shape is (32,)
	# Calculate for x=1100m forest entry
	u_1100 = nc_u.variables['u'][i, :z_max_index, :, 1]  # x=1100m, u shape is (32,32)
	nh3_1050 = nc_nh3.variables['nh3'][i, :z_max_index, :, 0]  # NH3 at x=1050, shape is (32,32)
	nh3_1150 = nc_nh3.variables['nh3'][i, :z_max_index, :, 1]  # NH3 at x=1150, shape is (32,32)
	nh3_1100 = 0.5 * (nh3_1050 + nh3_1150)  # simple average
	flux_1100 = np.sum(rho[:, np.newaxis] * u_1100 * nh3_1100) * dy * dz * 17.031 / 28.
	# rho[:,np.newaxis] changes shape from (32,) to (32,1) to match u_1100 and nh3_1100 shapes

	# Calculate for x=2200m forest exit
	u_2200 = nc_u.variables['u'][i, :z_max_index, :, 4]  # x=2200m, u shape is (32,32)
	nh3_2150 = nc_nh3.variables['nh3'][i, :z_max_index, :, 3]  # NH3 at x=2150, shape is (32,32)
	nh3_2250 = nc_nh3.variables['nh3'][i, :z_max_index, :, 4]  # NH3 at x=2250, shape is (32,32)
	nh3_2200 = 0.5 * (nh3_2150 + nh3_2250)  

	flux_2200 = np.sum(rho[:, np.newaxis] * u_2200 * nh3_2200) * dy * dz * 17.031 / 28.

	flux1100.append(flux_1100)
	flux2200.append(flux_2200)

	# rho[:,np.newaxis] changes shape from (32,) to (32,1) to match u_2200 and nh3_2200 shapes
	#print(f"NH3 mass flux at x = 1100m: {flux_1100:.6e} kg[NH3]/s")
	#print(f"NH3 mass flux at x = 2200m: {flux_2200:.6e} kg[NH3]/s")
# Close files
nc_def.close()
nc_u.close()
nc_nh3.close()
f,ax = plt.subplots()
ax.plot(flux1100)
ax.plot(flux2200)

