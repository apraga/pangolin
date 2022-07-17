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
  reduced = lon_prev
  translate = False
  #print sign
  if (lon_prev > np.pi):
    translate = True
    reduced = 2*np.pi - lon_prev

  delta = np.cos(reduced*0.5)*np.cos(lat_prev)/np.cos(lat_next)
  lon_next = 2*np.arccos(delta)
  if (sign == -1):
      lon_next = 2*np.pi - lon_next
  return lon_next

# Sign = 1 : arrival is in first half (longitude < 180)
# Sign = -1 : arrival is in second half (longitude > 180)
def get_sign(t2,t1,lat_prev,lon_prev,psi,U_0):
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
  zone = -1.
  if (lon_prev > np.pi):
    zone = 1.
  delta_t = np.arcsin(zone*value) - np.arcsin(np.sin(lat_prev)/beta) 

  delta_t = abs(delta_t)/alpha
  sign = 1.
  # Takes the modulo
  period = 2*np.pi/alpha
  diff = (t2-t1) % period
  #print "Dt, real",delta_t,diff
  # Check if we go to the other half, without going back on the current half
  if (diff > delta_t and diff < period*0.5 + delta_t):
    sign = -1.
  if (lon_prev > np.pi):
    sign = -sign
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
  # Change sense if we are in the other region (longitude > 180 degree)
  if (lon_prev > np.pi):
    alpha = -alpha
  sin_tmp = beta*np.sin(arcsin_tmp - alpha*(t2-t1))
  lat_prev = np.arcsin(sin_tmp)
  return lat_prev



###########################################

lon_prev = 1
U_0 = 1./174
t1 = 0
dt = 20.
nb_iter = 1

nb_lat = 80#179
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
#for i in range(1,2):#nb_lat):

name = "../plot/data_python/1.dat"
f = open(name,'w')

#for i in range(nb_lat):
  #for j in range(180,nb_lon):
i = 44
j = 200
k = (i-1)*nb_lon + j
lat_prevr = (90. - i*dlat)*np.pi/180
lon_prevr = j*dlon*np.pi/180
psi = get_psi(lat_prevr,lon_prevr,U_0)
if (abs(psi) < 1e-10):
  period = 0
else:
  alpha = np.sqrt(psi*U_0)
  period = int(2*np.pi/(alpha*dt))

for t in range(10):#period+1):
  t2 = t*dt
  #if (i == 90):
  #  lat_nextr = lat_prevr
  #  lon_nextr = lon_prevr
  #else:
  lat_nextr = get_lat_next(t2,t1,psi,lat_prevr,lon_prevr,U_0)
  sign = get_sign(t2,t1,lat_prevr,lon_prevr,psi,U_0)
  #if (j == 181):
    #print "latnext",lat_nextr 
    #print sign
  lon_nextr = get_lon_next(sign,lat_nextr,lat_prevr,lon_prevr,U_0)

  lat_next = lat_nextr*180/np.pi
  lon_next = lon_nextr*180/np.pi
  #print lat_next,lon_next
  i_lon = int(lon_next/dlon)
  i_lat = int((90.-lat_next)/dlat)
  ind = int((i_lat-1)*nb_lon + i_lon)
  lat_cur = 90. - i_lat*dlat
  lon_cur = i_lon*dlon
  s = "%f %f %f \n" % (lat_cur,lon_cur,q[k])
#  k = k + 1
  f.write(s)

f.close()
