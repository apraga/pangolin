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
  #lon_next = 2*np.arccos(sign*delta)
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
  sin_tmp = beta*np.sin(arcsin_tmp - alpha*(t2-t1))
  lat_prev = np.arcsin(sin_tmp)
  return lat_prev

def write_ncl(k,q,nb_lat,nb_lon):
  name = "../plot/data_python/"+str(k)+".dat"
  f = open(name,'w')
  dlat = 1.
  dlon = 1.
  k = 0
  for i in range(nb_lat):
    lat_cur = 90. - i*dlat
    for j in range(nb_lon):
      lon_cur = j*dlon
      s = "%f %f %f \n" % (lat_cur,lon_cur,q[k])
      k = k + 1
      f.write(s)
  f.close()



###########################################

lon_prev = 1
U_0 = 1./174
t1 = 0
dt = 20.
nb_iter = 1

nb_lat = 89
nb_lon = 360
nb_points = nb_lat*nb_lon
q = np.zeros(nb_points)
q2 = np.zeros(nb_points)
init = np.zeros(nb_points)

dlat = 1.
dlon = 1.
k = 0
t2 = nb_iter*dt
# Init
for i in range(nb_lat):
  for j in range(nb_lon):
    q[k] = i
    k = k + 1

k = 0
# Compute new position
for i in range(nb_lat):
  lat_prevr = (90. - i*dlat)*np.pi/180
  for j in range(nb_lon):
    if (i == 90):
      continue 
    lon_prevr = j*dlon*np.pi/180
    psi = get_psi(lat_prevr,lon_prevr,U_0)
    alpha = np.sqrt(psi*U_0)
    lat_nextr = get_lat_next(t2,t1,psi,lat_prevr,lon_prevr,U_0)
    sign = get_sign(t2,t1,lat_prevr,psi,U_0)
    lon_nextr = get_lon_next(sign,lat_nextr,lat_prevr,lon_prevr,U_0)
    lat_next = lat_nextr*180/np.pi
    lon_next = lon_nextr*180/np.pi
    if (math.isnan(lat_nextr)):
      print "nan : lat prev",i,j
    if (math.isnan(lon_nextr)):
      print "nan : lon prev",i,j

    if (lat_next < 0):
      continue

    i_lon = int(lon_next/dlon)
    i_lat = int((90.-lat_next)/dlat)
    ind = int((i_lat-1)*nb_lon + i_lon)
    #print "prev",lat_prevr*180/np.pi,lon_prevr*180/np.pi
    #print "next",lat_next,lon_next
    q2[ind] = q[k]
    k = k + 1

write_ncl(nb_iter,q2,nb_lat,nb_lon)
