# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 2025
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

# Set the folder path
folder = "/home/mohamad/ammonia/microhh/cases/plain/"

# Open datasets
dsxy = xr.open_dataset(folder+'nh3.xy.nc', decode_times=False)
dsxz = xr.open_dataset(folder+'nh3.xz.nc', decode_times=False)

# Create figure and axes
fig, ax = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

# Initialize first frame
# For xy view at z=50m (first z level)
img1 = ax[0].imshow(dsxy.nh3.isel(time=0, z=0).values,
                    extent=[dsxy.x[0], dsxy.x[-1], dsxy.y[-1], dsxy.y[0]])
plt.colorbar(img1, ax=ax[0], orientation='horizontal', label='NH3')

# For xz view
img2 = ax[1].imshow(dsxz.nh3.isel(time=0, y=0).values,
                    extent=[dsxz.x[0], dsxz.x[-1], dsxz.z[0], dsxz.z[-1]])
plt.colorbar(img2, ax=ax[1], orientation='horizontal', label='NH3')

# Set titles and labels
ax[0].set_title("Top view (xy) at z=50m")
ax[1].set_title("Side view (xz)")

ax[0].set_xlabel('x [m]')
ax[0].set_ylabel('y [m]')
ax[1].set_xlabel('x [m]')
ax[1].set_ylabel('z [m]')

# Add time annotation
time_text = fig.text(0.5, 0.95, '', ha='center')

# Update function for animation
def update(frame):
    # Update each view
    img1.set_array(dsxy.nh3.isel(time=frame, z=0).values)
    img2.set_array(dsxz.nh3.isel(time=frame, y=0).values)
    
    # Update time text
    time_hours = dsxy.time[frame].values / 3600
    time_text.set_text(f'Time = {time_hours:.2f} hours')
    
    # Update color scales
    vmin = min(dsxy.nh3.isel(time=frame, z=0).min(), 
              dsxz.nh3.isel(time=frame).min())
    vmax = max(dsxy.nh3.isel(time=frame, z=0).max(), 
              dsxz.nh3.isel(time=frame).max())
    
    img1.set_clim(vmin, vmax)
    img2.set_clim(vmin, vmax)
    
    return [img1, img2, time_text]

# Create and save animation
plt.tight_layout()
n_frames = len(dsxy.time)
anim = FuncAnimation(fig, update, frames=n_frames, interval=100, blit=True)

# Save as MP4
writer = animation.FFMpegWriter(fps=10, bitrate=3000)
anim.save('nh3_animation.mp4', writer=writer)
plt.close()
