import os
import numpy as np
import xarray as xr

# Get the current working directory
current_dir = os.getcwd()

# Define the file name
file_name = "plume_chem.forest.0000000.nc"

# Combine the current directory with the file name
file_path = os.path.join(current_dir, file_name)

# Load the NetCDF file
data = xr.open_dataset(file_path, group="default")

# Extract tke and z variables
tke = data["tke"]  # Shape: (time, z)
z = data["z"]  # Shape: (z,)

# Calculate layer thickness (dz) using z
dz = np.diff(z, prepend=z[0])

# Compute the total TKE by integrating over z
# Multiply tke by dz at each z-level and sum across z
total_tke = (tke * dz).sum(dim="z")

# Convert to a pandas dataframe for visualization (optional)
total_tke_df = total_tke.to_pandas()

# Save or display the results
print(total_tke_df)

# Save the results to CSV in the current directory
output_file = os.path.join(current_dir, "total_tke_timeseries.csv")
total_tke_df.to_csv(output_file)

print(f"Total TKE time series saved to {output_file}")

