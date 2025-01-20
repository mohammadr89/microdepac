# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 2025
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import os

# Get current working directory
folder = os.getcwd() + os.sep

# Open datasets
ds_lai = xr.open_dataset(folder+'nh3.xy.nc', decode_times=False)
ds_no_lai = xr.open_dataset(folder+'nh3.xy_no_lai.nc', decode_times=False)

# Calculate global min and max values across all timesteps
vmin_global = min(
    ds_lai.nh3.isel(z=0).min(),
    ds_no_lai.nh3.isel(z=0).min()
)
vmax_global = max(
    ds_lai.nh3.isel(z=0).max(),
    ds_no_lai.nh3.isel(z=0).max()
)

print(f"Global range: {vmin_global:.2e} to {vmax_global:.2e}")

# Create figure and axes (vertical arrangement)
fig, ax = plt.subplots(2, 1, figsize=(8, 10), dpi=300)

# Initialize first frame with global limits
# For top panel (with LAI)
img1 = ax[0].imshow(ds_lai.nh3.isel(time=0, z=0).values,
                    extent=[ds_lai.x[0], ds_lai.x[-1], ds_lai.y[-1], ds_lai.y[0]],
                    vmin=vmin_global, vmax=vmax_global)
plt.colorbar(img1, ax=ax[0], orientation='horizontal', label='NH3')

# For bottom panel (without LAI)
img2 = ax[1].imshow(ds_no_lai.nh3.isel(time=0, z=0).values,
                    extent=[ds_no_lai.x[0], ds_no_lai.x[-1], ds_no_lai.y[-1], ds_no_lai.y[0]],
                    vmin=vmin_global, vmax=vmax_global)
plt.colorbar(img2, ax=ax[1], orientation='horizontal', label='NH3')

# Set titles and labels
ax[0].set_title("With LAI: NH3 at z=50m")
ax[1].set_title("No LAI: NH3 at z=50m")

# Add axis labels
for axis in ax:
    axis.set_xlabel('x [m]')
    axis.set_ylabel('y [m]')

# Add time annotation
time_text = fig.text(0.5, 0.95, '', ha='center')

# Update function for animation
def update(frame):
    # Update each view
    img1.set_array(ds_lai.nh3.isel(time=frame, z=0).values)
    img2.set_array(ds_no_lai.nh3.isel(time=frame, z=0).values)
    
    # Update time text
    time_hours = ds_lai.time[frame].values / 3600
    time_text.set_text(f'Time = {time_hours:.2f} hours')
    
    # No need to update color limits as they're fixed
    return [img1, img2, time_text]

# Create and save animation
plt.tight_layout()
n_frames = len(ds_lai.time)
anim = FuncAnimation(fig, update, frames=n_frames, interval=100, blit=True)

# Save as MP4
writer = animation.FFMpegWriter(fps=10, bitrate=3000)
anim.save('nh3_xy_comparison_vertical_animation.mp4', writer=writer)
plt.close()
