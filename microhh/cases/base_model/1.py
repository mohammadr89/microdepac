import netCDF4 as nc    
import numpy as np      
import matplotlib.pyplot as plt

# Open files
nc_def = nc.Dataset('plume_chem.default.0010800.nc', 'r')  
nc_u = nc.Dataset('u.yz.nc', 'r')                          
nc_nh3 = nc.Dataset('nh3.yz.nc', 'r')                     

# Grid spacing
dy = dz = 100  
flux1100 = []
flux2200 = []
deposition_flux = []
vd_nh3 = []

for i in range(360):
    # Get density and calculate mass fluxes as before
    rho = nc_def.groups['thermo'].variables['rho'][i,:]
    
    # Calculate horizontal fluxes
    u_1100 = nc_u.variables['u'][i, :, :, 1]  
    nh3_1050 = nc_nh3.variables['nh3'][i, :, :, 0]  
    nh3_1150 = nc_nh3.variables['nh3'][i, :, :, 1]  
    nh3_1100 = 0.5 * (nh3_1050 + nh3_1150)  
    flux_1100 = np.sum(rho[:, np.newaxis] * u_1100 * nh3_1100) * dy * dz * 17.031 / 28.
    
    # Get deposition data
    dep_vel = nc_def.groups['deposition'].variables['vdnh3'][i]
    dep_flux = nc_def.groups['deposition'].variables['flux_nh3'][i]
    
    flux1100.append(flux_1100)
    vd_nh3.append(dep_vel)
    deposition_flux.append(dep_flux)

# Convert to arrays
flux1100 = np.array(flux1100)
vd_nh3 = np.array(vd_nh3)
deposition_flux = np.array(deposition_flux)

# Plot results
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))

# Plot horizontal mass flux
ax1.plot(flux1100, label='Horizontal flux at x=1100m')
ax1.set_ylabel('Mass flux (kg/s)')
ax1.legend()
ax1.grid(True)

# Plot deposition velocity
ax2.plot(vd_nh3, label='NH3 deposition velocity')
ax2.set_ylabel('Velocity (m/s)')
ax2.legend()
ax2.grid(True)

# Plot deposition flux
ax3.plot(deposition_flux, label='NH3 deposition flux')
ax3.set_ylabel('Flux (mol m-2 s-1)')
ax3.set_xlabel('Time step')
ax3.legend()
ax3.grid(True)

plt.tight_layout()
plt.show()

# Print summary statistics
print("Summary Statistics:")
print(f"Mean horizontal flux: {np.mean(flux1100):.2e} kg/s")
print(f"Mean deposition velocity: {np.mean(vd_nh3):.2e} m/s")
print(f"Mean deposition flux: {np.mean(deposition_flux):.2e} mol m-2 s-1")
print(f"\nTotal deposited NH3: {np.sum(deposition_flux):.2e} mol m-2")

# Close files
nc_def.close()
nc_u.close()
nc_nh3.close()
