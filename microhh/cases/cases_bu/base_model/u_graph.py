import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Read NetCDF file
u_yz_nc = nc.Dataset('u.yz.nc', 'r')

# Get dimensions
time = u_yz_nc.variables['time'][:]  # time array
x = u_yz_nc.variables['xh'][:]       # x cross-sections
y = u_yz_nc.variables['y'][:]        # y coordinates
z = u_yz_nc.variables['z'][:]        # z coordinates

# Show available cross-sections
print("Available cross-sections (x locations in meters):")
print(x)

# Get user input
x_target = float(input("Choose an x location from the above list (in meters): "))

# Find the index for chosen x
x_index = np.where(x == x_target)[0][0]

# Select a point in the y-z plane
y_index = 16  # Middle y point (32/2)
z_index = 0  # Middle z point (32/2)

# Extract u values over time at this point
u_time = u_yz_nc.variables['u'][:, z_index, y_index, x_index]

# Create the plot
plt.figure(figsize=(12, 6))
plt.plot(time/3600, u_time, 'b-', linewidth=2)  # Convert time to hours
plt.grid(True)
plt.xlabel('Time (hours)')
plt.ylabel('u velocity (m/s)')
plt.title(f'U velocity vs Time\nCross-section x={x_target}m, y={y[y_index]}m, z={z[z_index]}m')

# Add a horizontal line at u=0 for reference
plt.axhline(y=0, color='r', linestyle='--', alpha=0.5)

plt.show()

# Print information about the data
print(f"\nData Information:")
print(f"Cross-section at x = {x_target}m")
print(f"Location in y-z plane: y={y[y_index]}m, z={z[z_index]}m")
print(f"Time range: {time[0]/3600:.1f} to {time[-1]/3600:.1f} hours")
print(f"U velocity range: {np.min(u_time):.3f} to {np.max(u_time):.3f} m/s")

# Close NetCDF file
u_yz_nc.close()
