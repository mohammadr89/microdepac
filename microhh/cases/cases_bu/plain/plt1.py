# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 2025
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# Set the folder path
folder = "/home/mohamad/ammonia/microhh/cases/plain/"

# Open the datasets
dsxy = xr.open_dataset(folder+'nh3.xy.nc', decode_times=False)
dsxz = xr.open_dataset(folder+'nh3.xz.nc', decode_times=False)

# Create a list of datasets (now only 2)
ds = [dsxy, dsxz]

# Set time to last timestep
time = -1

# Create figure with two panels
fig, ax = plt.subplots(1, 2, dpi=300, figsize=(12, 5))

# For xy view at z=50m
img = ax[0].imshow(dsxy.nh3.isel(time=time, z=0).values,
                   extent=[dsxy.x[0], dsxy.x[-1], dsxy.y[-1], dsxy.y[0]])
plt.colorbar(img, ax=ax[0], orientation='horizontal', label='NH3')

# For xz view
img = ax[1].imshow(dsxz.nh3.isel(time=time, y=0).values,
                   extent=[dsxz.x[0], dsxz.x[-1], dsxz.z[0], dsxz.z[-1]])
plt.colorbar(img, ax=ax[1], orientation='horizontal', label='NH3')

# Set titles and labels
ax[0].set_title("Top view (xy) at z=50m")
ax[1].set_title("Side view (xz)")

# Add axis labels
ax[0].set_xlabel('x [m]')
ax[0].set_ylabel('y [m]')
ax[1].set_xlabel('x [m]')
ax[1].set_ylabel('z [m]')

# Add time in hours
time_hours = dsxy.time[time].values / 3600
fig.suptitle(f"time = {time_hours:.2f} hours")

# Adjust layout and save
fig.tight_layout()
fig.savefig("nh3_snapshot.png")
plt.close()
