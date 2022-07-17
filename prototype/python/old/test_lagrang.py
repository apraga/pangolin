import numpy as np
import math

# Input in radian
def get_psi(lat,lon,U_0):
  res = U_0*np.cos(lon*0.5)**2
  res = res*np.cos(lat)**2
  return res

def to_radian(lat,lon):
  lat2 = lat*np.pi/180
  lon2 = lon*np.pi/180
  return lat2,lon2

def to_degree(lat,lon):
  lat2 = int(lat*180/np.pi)
  lon2 = int(lon*180/np.pi)
  return lat2,lon2

# Input in radian
def get_lon_prev(psi,lat_prev,lat_cur,lon_cur,U_0):
  # Compute longitude at t=0 (longitude is always positive)
  delta = np.cos(lon_cur*0.5)*np.cos(lat_cur)/np.cos(lat_prev)
  lon_prev = 2*np.arccos(delta)
  return lon_prev

###########################################

lat0 = 44
lon0 = 1
U_0 = 1./174
lat0r,lon0r = to_radian(lat0,lon0)
psi = get_psi(lat0r,lon0r,U_0)

name = "data.dat"
f = open(name,'w')

for lat in range(43,-44,-1):
  latr = lat*np.pi/180
  lon_prevr = get_lon_prev(psi,latr,lat0r,lon0r,U_0)
  lon_prev = lon_prevr*180./np.pi
  s = "%f %f \n" % (lon_prev,lat)
  f.write(s)

f.close()
