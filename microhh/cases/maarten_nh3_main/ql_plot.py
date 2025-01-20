import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

# Read the netCDF file
#ds = nc.Dataset('plume_chem.default.0000000.nc')
ds = nc.Dataset('plume_chem.default.0010800.nc')

# Get the variables
z = ds.variables['z'][:]  # height levels
ql = ds.groups['thermo'].variables['ql'][:]  # liquid water
time = ds.variables['time'][:]  # time steps

# Create figure with a specific size and with space for colorbar
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111)

# Plot for each time step with different colors
num_timesteps = len(time)
colors = plt.cm.rainbow(np.linspace(0, 1, num_timesteps))  # Create color gradient

for t in range(num_timesteps):
    ax.plot(ql[t, :], z, '-', color=colors[t], linewidth=1, alpha=0.7)

# Customize the plot
ax.set_xlabel('ql [kg/kg]')
ax.set_ylabel('z [m]')
ax.set_title('Liquid Water Profile')
ax.grid(True, linestyle='--', alpha=0.5)

# Set axis limits
ax.set_ylim(0, 4000)
min_ql = np.min(ql)
max_ql = np.max(ql)
ax.set_xlim(min_ql - 0.1*(max_ql-min_ql), max_ql + 0.1*(max_ql-min_ql))

# Scientific notation for x-axis since these are typically small values
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

# Add colorbar to show time progression
norm = plt.Normalize(time[0], time[-1])
sm = plt.cm.ScalarMappable(cmap=plt.cm.rainbow, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax)
cbar.set_label('Time [s]')

# Adjust layout
plt.tight_layout()

# Save the plot
plt.savefig('ql_profile_timesteps.png', dpi=300, bbox_inches='tight')

# Print some basic statistics
print(f"Time range: {time[0]:.1f} to {time[-1]:.1f} seconds")
print(f"Height range: {z.min():.2f} to {z.max():.2f} m")
print(f"Ql range: {ql.min():.6f} to {ql.max():.6f} kg/kg")
print(f"Number of time steps: {num_timesteps}")

# Show the plot
plt.show()

ds.close()
