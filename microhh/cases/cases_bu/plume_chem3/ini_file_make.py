import os

# Set the file path relative to the current directory
file_path = os.path.join(os.getcwd(), 'plume_chem.ini')

# Load the ini file content
with open(file_path, 'r') as file:
    ini_content = file.read()

# Create two new ini contents with modified starttime and endtime
ini_content_1 = ini_content.replace('starttime=0', 'starttime=0').replace('endtime=7200', 'endtime=7200')
ini_content_2 = ini_content.replace('starttime=0', 'starttime=7200').replace('endtime=7200', 'endtime=21000')

# Define the new file paths
file_1_path = os.path.join(os.getcwd(), 'plume_chem1.ini')
file_2_path = os.path.join(os.getcwd(), 'plume_chem2.ini')

# Write the new contents to respective files
with open(file_1_path, 'w') as file:
    file.write(ini_content_1)

with open(file_2_path, 'w') as file:
    file.write(ini_content_2)

print(f"Files created:\n1. {file_1_path}\n2. {file_2_path}")

