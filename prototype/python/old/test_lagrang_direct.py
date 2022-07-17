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
# Sign is the sign of the longitude 
def get_lon_next(sign,lat_next,lat_prev,lon_prev,U_0):
  delta = np.cos(lon_prev*0.5)*np.cos(lat_prev)/np.cos(lat_next)
  lon_next = 2*np.arccos(sign*delta)
  return lon_next

def get_sign(t2,t1,lat_prev,psi,U_0):
  # Compute time for changing zones
  a = 1.
  beta = np.sqrt(1 - psi/(a*U_0))
  alpha = np.sqrt(psi*U_0/(a**3))
  lat_max = np.arccos(np.sqrt(psi/(a*U_0)))
  value = np.sin(lat_max)/beta
  threshold = 1.0
  threshold_value = 0.99999999
  if (value > threshold):
    value = threshold_value
  if (value < -threshold):
    value = threshold_value
  #delta_t1 = np.arcsin(value) - np.arcsin(np.sin(lat_prev)/beta) 
  delta_t = np.arcsin(-value) - np.arcsin(np.sin(lat_prev)/beta) 
  delta_t = -delta_t/alpha
  sign = 1.
  #print "Dt, real",delta_t,t2-t1
  if ((t2-t1) > delta_t):
    sign = -1.
  return sign

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
  sin_tmp = beta*np.sin(arcsin_tmp - alpha*(t2-t1))
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

name = "data.dat"
f = open(name,'w')

for t in range(3,80):
  t2 = t*dt
  lat_nextr = get_lat_next(t2,t1,psi,lat_prevr,lon_prevr,U_0)
  sign = get_sign(t2,t1,lat_prevr,psi,U_0)
  lon_nextr = get_lon_next(sign,lat_nextr,lat_prevr,lon_prevr,U_0)
  lat_next = lat_nextr*180/np.pi
  lon_next = lon_nextr*180/np.pi
  s = "%f %f \n" % (lon_next,lat_next)
  f.write(s)
#  print lon_next,lat_next
f.close()
