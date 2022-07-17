import numpy as np
import math

# Input in radian
def get_psi(lat,lon,U_0):
  res = U_0*np.cos(lon*0.5)**2
  res = res*np.cos(lat)**2
  return res

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

def to_radian(lat,lon):
  lat2 = lat*np.pi/180
  lon2 = lon*np.pi/180
  return lat2,lon2

def to_degree(lat,lon):
  lat2 = int(lat*180/np.pi)
  lon2 = int(lon*180/np.pi)
  return lat2,lon2

# Input in radian
def get_prev(t2,t1,lat_cur,lon_cur,U_0):
  a = 1.
  # Compute latitude at t=0
  if (lat_cur == 0 and lon_cur == 0):    # Equator
    lat_prev = lat_cur 
    lon_prev = 2*np.arccos(1./np.cos(lat_prev))
  else:
    psi = get_psi(lat_cur,lon_cur,U_0)
    beta = np.sqrt(1 - psi/(a*U_0))
    alpha = np.sqrt(psi*U_0/(a**3))
    value = np.sin(lat_cur)/beta
    threshold = 1.0
    threshold_value = 0.99999999
    if (value > threshold):
      print value,"> 1"
      value = threshold_value
    if (value < -threshold):
      print value,"< -1"
      value = threshold_value
    arcsin_tmp = np.arcsin(value)
    sin_tmp = beta*np.sin(arcsin_tmp + alpha*(t2-t1))
    lat_prev = np.arcsin(sin_tmp)

    # Compute longitude at t=0 (longitude is always positive)
    delta = np.cos(lon_cur*0.5)*np.cos(lat_cur)/np.cos(lat_prev)
    lon_prev = 2*np.arccos(delta)

  
  return lat_prev,lon_prev

###########################################

nb_lat = 90
#nb_lon = 180#360
nb_lon = 360
U_0 = 1./174
nb_points = nb_lat*nb_lon
q = np.zeros(nb_points)
q2 = np.zeros(nb_points)
init = np.zeros(nb_points)

dlat = 1.
dlon = 1.
k = 0
# Init
for i in range(nb_lat):
  for j in range(nb_lon):
    q[k] = i
    k = k + 1

# Iteration
dt = 20
nb_iter = 30
t1 = 0
t2 = nb_iter*dt
k = 0
for i in range(nb_lat):
  lat_cur = 89. - i*dlat
  for j in range(nb_lon):
    lon_cur = j*dlon
    # Get previous position
    lat_cur2,lon_cur2 = to_radian(lat_cur,lon_cur)
    lat_prev2,lon_prev2 = get_prev(t2,t1,lat_cur2,lon_cur2,U_0)
    lat_prev,lon_prev = to_degree(lat_prev2,lon_prev2)
    if (lat_prev < 0):
      lat_prev = 0

    # Search the position in the array
    i_lat = int((90.-lat_prev)/dlat)
    #print i_lat,(90.-lat_prev)/dlat
    i_lon = int(lon_prev/dlon)
    k_prev = (i_lat-1)*nb_lon + i_lon
    #print "cur :",lat_cur,lon_cur
    #print "prev :",lat_prev,lon_prev
    #print "ind :",i_lat,i_lon
    #print "lon, lon afetr:",i_lon*dlon,(i_lon+1)*dlon
    psi_cur = get_psi(lat_cur2,lon_cur2,U_0)
    psi_prev = get_psi(lat_prev2,lon_prev2,U_0)
    #print k,k_prev
    #lon_tmp = i_lon*dlon
    #lat_tmp = 90. - i_lat*dlat
    #lat_tmp = i_lat*dlat
    q2[k] = q[k_prev]
    k = k + 1

write_ncl(nb_iter,q2,nb_lat,nb_lon)
