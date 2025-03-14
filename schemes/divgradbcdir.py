from pylab import *

xsize = 0.5
igc   = 3

def geterror(u, uref, n=2):
  error = sqrt(sum((uref - u)**2.) / u.size)
  return error

def refdata(n):
  x = linspace(0.5*xsize/n, xsize-0.5*xsize/n, n)
  u = sin(2.*pi*x/xsize)
  return x, u

def dgx2nd(x, u):
  istart = igc
  iend   = u.size + igc
  dx     = xsize / u.size
  dgu    = zeros(u.size)
  ucalc  = zeros(u.size + 2*igc)
  ucalc[istart:iend] = u[:]

  # # periodic bcs
  # ucalc[0   :igc     ] = u[u.size-igc:u.size]
  # ucalc[iend:iend+igc] = u[0:igc]

  # non-periodic bc
  ucalc[istart-1] = - ucalc[istart]
  ucalc[iend    ] = - ucalc[iend-1]

  for i in range(istart, iend):
    dgu[i-igc] = (ucalc[i-1] - 2.*ucalc[i] + ucalc[i+1]) / (dx**2.)

  dguref = -4*pi**2./xsize**2. * sin(2.*pi*x/xsize)

  erru = geterror(dgu, dguref)

  return dgu, dguref, erru

def dgx4th(x, u):
  istart = igc
  iend   = u.size + igc
  dx     = xsize / u.size
  dgu    = zeros(u.size)
  ucalc  = zeros(u.size + 2*igc)
  ucalc[istart:iend] = u[:]

  # periodic bcs
  #ucalc[0   :igc     ] = u[u.size-igc:u.size]
  #ucalc[iend:iend+igc] = u[0:igc]

  # ghost cell
  ucalc[istart-1] = (-15.*ucalc[istart] + 5.*ucalc[istart+1] - ucalc[istart+2]) / 5.
  ucalc[iend    ] = (-15.*ucalc[iend-1] + 5.*ucalc[iend-2  ] - ucalc[iend-3  ]) / 5.

  i = istart
  dgu[i-igc] = (550.*ucalc[i-1] - 1047.*ucalc[i] + 416.*ucalc[i+1] + 110.*ucalc[i+2] - 30.*ucalc[i+3] + ucalc[i+4]) / (576.*dx**2.)
  i = istart+1
  dgu[i-igc] = (-50.*ucalc[i-2] + 777.*ucalc[i-1] - 1456.*ucalc[i] + 782.*ucalc[i+1] - 54.*ucalc[i+2] + ucalc[i+3]) / (576.*dx**2.)
  for i in range(istart+2, iend-2):
    dgu[i-igc] = ((ucalc[i-3]+ucalc[i+3]) -54.*(ucalc[i-2]+ucalc[i+2]) + 783.*(ucalc[i-1]+ucalc[i+1]) - 1460.*ucalc[i]) / (576.*dx**2.)
  i = iend-2
  dgu[i-igc] = (-50.*ucalc[i+2] + 777.*ucalc[i+1] - 1456.*ucalc[i] + 782.*ucalc[i-1] - 54.*ucalc[i-2] + ucalc[i-3]) / (576.*dx**2.)
  i = iend-1
  dgu[i-igc] = (550.*ucalc[i+1] - 1047.*ucalc[i] + 416.*ucalc[i-1] + 110.*ucalc[i-2] - 30.*ucalc[i-3] + ucalc[i-4]) / (576.*dx**2.)

  dguref = -4*pi**2./xsize**2. * sin(2.*pi*x/xsize)

  erru = geterror(dgu, dguref)

  return dgu, dguref, erru


x8, u8 = refdata(8)
dx8    = xsize / 8.
dg8_2nd, dgref8_2nd, err8_2nd = dgx2nd(x8, u8)
dg8_4th, dgref8_4th, err8_4th = dgx4th(x8, u8)

x16, u16 = refdata(16)
dx16     = xsize / 16.
dg16_2nd, dgref16_2nd, err16_2nd = dgx2nd(x16, u16)
dg16_4th, dgref16_4th, err16_4th = dgx4th(x16, u16)

x32, u32 = refdata(32)
dx32     = xsize / 32.
dg32_2nd, dgref32_2nd, err32_2nd = dgx2nd(x32, u32)
dg32_4th, dgref32_4th, err32_4th = dgx4th(x32, u32)

x64, u64 = refdata(64)
dx64     = xsize / 64.
dg64_2nd, dgref64_2nd, err64_2nd = dgx2nd(x64, u64)
dg64_4th, dgref64_4th, err64_4th = dgx4th(x64, u64)

x128, u128 = refdata(128)
dx128      = xsize / 128.
dg128_2nd, dgref128_2nd, err128_2nd = dgx2nd(x128, u128)
dg128_4th, dgref128_4th, err128_4th = dgx4th(x128, u128)

dxs  = array([dx8, dx16, dx32, dx64, dx128])
errs_2nd = array([err8_2nd, err16_2nd, err32_2nd, err64_2nd, err128_2nd])
errs_4th = array([err8_4th, err16_4th, err32_4th, err64_4th, err128_4th])

print('convergence 2nd', (log(errs_2nd[-1])-log(errs_2nd[0])) / (log(dxs[-1])-log(dxs[0])) )
print('convergence 4th', (log(errs_4th[-1])-log(errs_4th[0])) / (log(dxs[-1])-log(dxs[0])) )

off2 = 8.
off3 = 4.
off4 = 2.
slope2 = off2*(dxs[:] / dxs[0])**2.
slope3 = off3*(dxs[:] / dxs[0])**3.
slope4 = off4*(dxs[:] / dxs[0])**4.

close('all')

figure()
plot(x8 , dg8_2nd , 'b-',  label="8_2nd" )
plot(x16, dg16_2nd, 'g-',  label="16_2nd")
plot(x32, dg32_2nd, 'r-',  label="32_2nd")
plot(x8 , dg8_4th , 'b--', label="8_4th" )
plot(x16, dg16_4th, 'g--', label="16_4th")
plot(x32, dg32_4th, 'r--', label="32_4th")
#plot(x64, dg64_4th, 'c-', label="64_4th")
plot(x128, dgref128_4th, 'k--', label="ref")
legend(loc=4, frameon=False)

figure()
loglog(dxs, errs_2nd, 'bo-', label="g2nd")
loglog(dxs, errs_4th, 'go-', label="g4th")
loglog(dxs, slope2, 'k--', label="2nd")
loglog(dxs, slope3, 'k-.', label="3rd")
loglog(dxs, slope4, 'k:' , label="4th")
legend(loc=4, frameon=False)
#xlim(0.01, 0.2)

