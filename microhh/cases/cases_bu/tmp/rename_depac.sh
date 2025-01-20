#!/bin/bash

# Directory containing the files
DIRECTORY="."

# Loop through all files with the pattern *.xy.nc
for file in "$DIRECTORY"/*.xy.nc; do
    # Check if the file exists
    if [[ -f "$file" ]]; then
        # Create the new file name by replacing ".xy.nc" with "_depac.nc"
        new_name="${file%.xy.nc}.xy_depac.nc"
        
        # Rename the file
        mv "$file" "$new_name"
    fi
done

echo "Files renamed successfully!"

