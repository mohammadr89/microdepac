import netCDF4 as nc

# Function to print structure of NetCDF file
def print_ncfile_structure(filename):
    dataset = nc.Dataset(filename)
    print("\nFile:", filename)
    print("\nDimensions:")
    for dim in dataset.dimensions:
        print(f"{dim}: {len(dataset.dimensions[dim])}")
    
    print("\nVariables:")
    for var in dataset.variables:
        print(f"{var}: {dataset.variables[var].dimensions}")
    
    print("\nGroups:")
    for group in dataset.groups:
        print(f"\nGroup: {group}")
        print("Dimensions:")
        for dim in dataset.groups[group].dimensions:
            print(f"{dim}: {len(dataset.groups[group].dimensions[dim])}")
        print("Variables:")
        for var in dataset.groups[group].variables:
            print(f"{var}: {dataset.groups[group].variables[var].dimensions}")
    
    dataset.close()

print_ncfile_structure("moistcblles_input.nc")
