import numpy as np
import math

# Input in radian
def get_psi(lat,lon,U_0):
  res = U_0*np.cos(lon*0.5)**2
  res = res*np.cos(lat)**2
  return res

# Input in radian
# Sign is the sign of the longitude 
def get_lon_prev(sign,lat_prev,lat_cur,lon_cur,U_0):
  reduced = lon_cur
  translate = False
  #print sign
  if (lon_cur> np.pi):
    translate = True
    reduced = 2*np.pi - lon_cur

  delta = np.cos(reduced*0.5)*np.cos(lat_cur)/np.cos(lat_prev)
  lon_prev = 2*np.arccos(delta)
  if (sign == -1):
      lon_prev = 2*np.pi - lon_prev
  return lon_prev

# Sign = 1 : arrival is in first half (longitude < 180)
# Sign = -1 : arrival is in second half (longitude > 180)
def get_sign(t2,t1,lat_cur,lon_cur,psi,U_0):
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
  zone = 1.
  if (lon_cur > np.pi):
    zone = -1.
  delta_t = np.arcsin(zone*value) - np.arcsin(np.sin(lat_cur)/beta) 

  delta_t = abs(delta_t)/alpha
  sign = 1.
  # Takes the modulo
  period = 2*np.pi/alpha
  diff = (t2-t1) % period
  #print "Dt, real",delta_t,diff
  # Check if we go to the other half, without going back on the cur_rent half
  if (diff > delta_t and diff < period*0.5 + delta_t):
    sign = -1.
  if (lon_cur > np.pi):
    sign = -sign
  return sign

def get_lat_prev(t2,t1,psi,lat_cur,lon_cur,U_0):
  a = 1.
  beta = np.sqrt(1 - psi/(a*U_0))
  alpha = np.sqrt(psi*U_0/(a**3))
  value = np.sin(lat_cur)/beta
  threshold = 1.0
  threshold_value = 0.99999999
  if (value > threshold):
    #print value,"> 1"
    value = threshold_value
  if (value < -threshold):
    #print value,"< -1"
    value = threshold_value
  arcsin_tmp = np.arcsin(value)
  # Change sense if we are in the other region (longitude > 180 degree)
  if (lon_cur > np.pi):
    alpha = -alpha
  sin_tmp = beta*np.sin(arcsin_tmp + alpha*(t2-t1))
  lat_cur = np.arcsin(sin_tmp)
  return lat_cur

def get_nb_mesh(i,nb_lat,nb_lat2):
  if (i > nb_lat2):
    i2 = nb_lat + 1 - i
  else:
    i2 = i
  nb_mesh = 3*(2*i2-1)
  return nb_mesh

###########################################

lon_prev = 1
U_0 = 1./174
t1 = 0
dt = 1.
nb_iter = 10

nb_lat = 180#89#179
nb_lat2 = 90
nb_lon = 360
nb_points = 6*nb_lat2*nb_lat2
q = np.zeros(nb_points)

dlat = 1.
dlon = 1.
k = 0
t2 = nb_iter*dt

# Init
name1 = "../plot/data/0.dat"
f = open(name1,'r')
for line in f.readlines():
  tmp = line.split()
  q[k] = tmp[2]
  k = k +1
f.close()
  
# Write init
#name = "../plot/data_python/0.dat"
#f = open(name,'w')
#k = 0
#for i in range(1,nb_lat+1):
#  nb_mesh = get_nb_mesh(i,nb_lat,nb_lat2)
#  dlon = 360./nb_mesh
#  lat_cur = (90. - i*dlat)
#  for j in range(1,nb_mesh+1):
#    lon_cur = j*dlon-180
#    #print lon_cur
#    s = "%f %f %f \n" % (lat_cur,lon_cur,q[k])
#    k = k + 1
#    f.write(s)
#
coef = np.pi/180

name = "../plot/data_python/"+str(nb_iter)+".dat"
f = open(name,'w')

#for i in range(1,nb_lat,5):
for i in range(1,nb_lat):
  lat_cur = (90. - i*dlat)
  lat_cur_r = lat_cur*coef
  nb_mesh = get_nb_mesh(i,nb_lat,nb_lat2)
  dlon = 360./nb_mesh
  for j in range(1,nb_mesh):
    lon_cur = j*dlon
    lon_cur_r = lon_cur*coef
    if (lat_cur_r == 0 and lon_cur_r == 0):
      lat_prev_r = lat_cur_r
      lon_prev_r = lon_cur_r
    else:
      psi = get_psi(lat_cur_r,lon_cur_r,U_0)
      
      lat_prev_r = get_lat_prev(t2,t1,psi,lat_cur_r,lon_cur_r,U_0)
      sign = get_sign(t2,t1,lat_cur_r,lon_cur_r,psi,U_0)
      lon_prev_r = get_lon_prev(sign,lat_prev_r,lat_cur_r,lon_cur_r,U_0)
    
    lat_prev = (lat_prev_r/coef)
    lon_prev = (lon_prev_r/coef)

    i_prev = int((90.-lat_prev)/dlat)
    nb_mesh_prev = get_nb_mesh(i_prev,nb_lat,nb_lat2)
    dlon_prev = 360./nb_mesh_prev
    j_prev = int(lon_prev/dlon_prev)
    if (i_prev > nb_lat2):
      i2 = nb_lat + 1 - i_prev
      k_prev = nb_points - 3*(i2-1)*(i2-1) - j_prev
    else:
      k_prev = 3*(i_prev-1)*(i_prev-1) + j_prev
    if (k_prev > nb_points-1):
      print "prev",lat_prev,lon_prev

    k= 3*(i-1)*(i-1) + j
    s = "%f %f %f \n" % (lat_cur,lon_cur-180.,q[k_prev])
    f.write(s)

f.close()
