import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

# Read the NetCDF files
depac_file = 'vdnh3_soil.xy_depac.nc'
ifs_file = 'vdnh3_soil.xy_ifs.nc'

# Open datasets with decode_times=False to handle non-standard time units
ds_depac = xr.open_dataset(depac_file, decode_times=False)
ds_ifs = xr.open_dataset(ifs_file, decode_times=False)

# Calculate difference (DEPAC - IFS)
difference = ds_depac['vdnh3_soil'] - ds_ifs['vdnh3_soil']

# Create a function to plot the difference at any timestep
def plot_difference_at_time(difference, time_index):
    plt.figure(figsize=(10, 8))
    
    # Extract data for the specific time
    data = difference.isel(time=time_index)
    time_value = float(difference.time[time_index])
    
    # Create the plot
    im = plt.pcolormesh(difference.x, difference.y, data, 
                        cmap='RdBu_r', shading='auto')
    
    # Add colorbar and labels
    plt.colorbar(im, label='Vegetation NH3 deposition velocity difference (DEPAC - IFS)')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title(f'Vegetation NH3 deposition velocity difference at time {time_value} s')
    
    plt.show()

# Plot for the last timestep
plot_difference_at_time(difference, -1)

# Calculate some basic statistics
stats = {
    'mean_difference': float(difference.mean()),
    'max_difference': float(difference.max()),
    'min_difference': float(difference.min()),
    'std_difference': float(difference.std())
}

print("\nStatistics of the difference (DEPAC - IFS):")
for key, value in stats.items():
    print(f"{key}: {value:.6f}")

# Save the difference to a new NetCDF file if needed
difference.to_netcdf('vdnh3_soil_difference.nc')
