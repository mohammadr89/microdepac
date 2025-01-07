import configparser
import shutil

def modify_ini_files():
    """
    Creates two modified copies of plume_chem.ini with different time settings.
    The second file's start time matches the first file's end time.
    """
    # First, create copies of the original file
    shutil.copy('plume_chem.ini', 'plume_chem1.ini')
    shutil.copy('plume_chem.ini', 'plume_chem2.ini')
    
    # Process first file
    config1 = configparser.ConfigParser()
    config1.read('plume_chem1.ini')
    
    # Modify time settings for first file
    config1['time'] = {
        'endtime': '10800',
        'dt': '6.0',
        'dtmax': '60.0',
        'savetime': '3600',
        'outputiter': '5',
        'adaptivestep': 'true',
        'starttime': '0',
        'rkorder': '3'
    }
    
    # Save first file
    with open('plume_chem1.ini', 'w') as configfile:
        config1.write(configfile)
    
    # Process second file
    config2 = configparser.ConfigParser()
    config2.read('plume_chem2.ini')
    
    # Modify time settings for second file
    # Note: starttime equals endtime of first file
    config2['time'] = {
        'endtime': '25200',
        'dt': '6.0',
        'dtmax': '60.0',
        'savetime': '3600',
        'outputiter': '5',
        'adaptivestep': 'true',
        'starttime': '10800',  # Matches endtime of first file
        'rkorder': '3'
    }
    
    # Save second file
    with open('plume_chem2.ini', 'w') as configfile:
        config2.write(configfile)
    
    print("Created two modified ini files:")
    print("1. plume_chem1.ini (endtime=10800, starttime=0)")
    print("2. plume_chem2.ini (endtime=25200, starttime=10800)")

if __name__ == "__main__":
    modify_ini_files()
