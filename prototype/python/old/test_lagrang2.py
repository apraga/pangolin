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
def get_lon_next(t1,lat_next,lat_prev,lon_prev,U_0):
  # Compute longitude at t=0 (longitude is always positive)
  delta = np.cos(lon_prev*0.5)*np.cos(lat_prev)/np.cos(lat_next)
  margin = 2*np.pi/180
  lon_next = 2*np.arccos(delta)
  lon_next2 = 2*np.arccos(-delta)
  #if (t1 > 0 and abs(lon_prev) < margin and lon_next > 0):
  #if (abs(lon_prev) < margin and abs(lon_next) < margin):
    #lon_next = -lon_next
  #print "test",2*np.arccos(-delta)*180/np.pi
  return lon_next,lon_next2

def get_lat_next(t2,t1,psi,lat_prev,lon_prev,U_0):
  a = 1.
  beta = np.sqrt(1 - psi/(a*U_0))
  alpha = np.sqrt(psi*U_0/(a**3))
  value = np.sin(lat_prev)/beta
  threshold = 1.0
  threshold_value = 0.99999999
  if (value > threshold):
    print value,"> 1"
    value = threshold_value
  if (value < -threshold):
    print value,"< -1"
    value = threshold_value
  arcsin_tmp = np.arcsin(value)
  # Get sign
  sign = 1.
  if (lon_prev < 0):
    sign = -1
  sin_tmp = beta*np.sin(arcsin_tmp - alpha*sign*(t2-t1))
  lat_prev = np.arcsin(sin_tmp)
  return lat_prev


###########################################

lat_prev = 44
lon_prev = 1
U_0 = 1./174
lat_prevr,lon_prevr = to_radian(lat_prev,lon_prev)
psi = get_psi(lat_prevr,lon_prevr,U_0)
t1 = 0
dt = 20.
t2 = dt

name = "data.dat"
f = open(name,'w')

for t in range(45):
  lat_nextr = get_lat_next(t2,t1,psi,lat_prevr,lon_prevr,U_0)
  lon_nextr,lon_nextr2 = get_lon_next(t1,lat_nextr,lat_prevr,lon_prevr,U_0)
  lat_next = lat_nextr*180/np.pi
  lon_next = lon_nextr*180/np.pi
  lon_next2 = lon_nextr2*180/np.pi
  s = "%f %f \n" % (lon_next,lat_next)
  f.write(s)
  s = "%f %f \n" % (lon_next2,lat_next)
  #print lon_next,lat_next
  f.write(s)
  t2 = t2 + dt
  t1 = t1 + dt
  lat_prevr = lat_nextr
  lon_prevr = lon_nextr

f.close()
