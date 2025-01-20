import os

# Set the file path relative to the current directory
file_path = os.path.join(os.getcwd(), 'plume_chem.ini')

# Load the ini file content
with open(file_path, 'r') as file:
    ini_content = file.read()

# Create three ini files with different time ranges
# plume_chem1.ini: 0 to 10800
ini_content_1 = ini_content.replace('endtime=21600', 'endtime=10800')

# plume_chem2.ini: 10800 to 21600  
ini_content_2 = ini_content.replace('starttime=0', 'starttime=10800').replace('endtime=21600', 'endtime=21600')


# Define the new file paths 
file_1_path = os.path.join(os.getcwd(), 'plume_chem1.ini')
file_2_path = os.path.join(os.getcwd(), 'plume_chem2.ini')

# Write the modified contents to new files
with open(file_1_path, 'w') as file:
    file.write(ini_content_1)
    
with open(file_2_path, 'w') as file:
    file.write(ini_content_2)
