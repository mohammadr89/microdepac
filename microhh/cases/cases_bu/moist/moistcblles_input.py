import matplotlib.pyplot as pl
import numpy as np
import netCDF4 as nc

float_type = "f8"
# float_type = "f4"

np_dtype = np.float64 if float_type == "f8" else np.float32

# Get number of vertical levels and size from .ini file
with open('moistcblles.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

# Parameters for the well-mixed profile
h = 1000.        # Height of the mixed layer (m)
dthl = 10.       # Temperature jump at inversion (K)
dthz = 100.      # Thickness of the inversion layer (m)
dthetadz = 0.003 # Temperature gradient above inversion (K/m)
dqtdz = -1.e-6   # Simple moisture gradient

# set the height
z   = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
u   = np.ones(np.size(z))*5
v   = np.zeros(np.size(z))
thl = np.zeros(np.size(z))
qt  = np.zeros(np.size(z))

# Create well-mixed profile with jump for temperature only
for k in range(kmax):
    # Temperature profile with jump
    if(z[k] <= h - 0.5*dthz):
        thl[k] = 300.  # Well-mixed layer
    elif(z[k] <= h + 0.5*dthz):
        thl[k] = 300. + dthl/dthz * (z[k]-(h-0.5*dthz))  # Inversion layer
    else:
        thl[k] = 300. + dthl + dthetadz*(z[k]-(h+0.5*dthz))  # Free atmosphere

    # Simple linear moisture profile
    qt[k] = 0.015 + dqtdz*z[k]

nc_file = nc.Dataset("moistcblles_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", kmax)
nc_z   = nc_file.createVariable("z"  , float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_u   = nc_group_init.createVariable("u"  , float_type, ("z"))
nc_v   = nc_group_init.createVariable("v"  , float_type, ("z"))
nc_thl = nc_group_init.createVariable("thl", float_type, ("z"))
nc_qt  = nc_group_init.createVariable("qt" , float_type, ("z"))

nc_z  [:] = z  [:]
nc_u  [:] = u  [:]
nc_v  [:] = v  [:]
nc_thl[:] = thl[:]
nc_qt [:] = qt [:]

nc_file.close()

# Create time varying surface fields with uniform values
endtime = 10800
dt = 3600
nt = int((endtime / dt)+1)

itot = 64
jtot = 32

thl_sbot = np.zeros((nt, jtot, itot), dtype=np_dtype)
qt_sbot = np.zeros((nt, jtot, itot), dtype=np_dtype)

# Set uniform surface values
thl_sbot[:] = 0.1       # Uniform temperature flux
qt_sbot[:] = 5.e-5      # Uniform moisture flux

# Write as binary input files for MicroHH
for t in range(nt):
    thl_sbot[t,:].tofile('thl_bot_in.{0:07d}'.format(t*dt))
    qt_sbot[t,:].tofile('qt_bot_in.{0:07d}'.format(t*dt))
