# Wrong way : from reduced to reduced
import numpy as np
import math
import sys
import matplotlib.pyplot as p
from mpl_toolkits.basemap import Basemap
import subprocess as sub

# Analytical bidimensionnal (hourdin)

#######################################################################
# Functions definitions
#######################################################################

def init_q(nb_lat,nb_lon):
  # Init : the concentration goes from 0 to nb_lon/2 to 0
  # on [0,360]
  k = 0
  nb_points = nb_lat*nb_lon
  q = np.zeros(nb_points)
  
  for i in range(nb_lat):
    for j in range(1,nb_lon/2+2):
      q[k] = j
      k = k + 1
    for j in range(nb_lon/2,1,-1):
      q[k] = j;
      k = k + 1
  return q

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



###############################################################
# The script
###############################################################

lon_prev = 1
U_0 = 1./174
t1 = 0
dt = 20.
#nb_iter = int(2*np.pi/U_0)
if (len(sys.argv) > 1):
  nb_iter = int(sys.argv[1])
else:
  print "Missing the number of iterations !"
  print "Usage : "+str(sys.argv[0])+" NB_ITER"
  exit(2)

nb_lat = 180
nb_lon = 360

dlat = 1.
dlon = 1.
t2 = nb_iter*dt
nb_points = nb_lat*nb_lon
q = init_q(nb_lat,nb_lon)
lat = np.zeros(nb_points,float)
lon = np.zeros(nb_points,float)
q_new = np.zeros(nb_points,float)
coef = np.pi/180

k = 0
for i in range(nb_lat-1):
  lat_cur = (89. - i*dlat)
  lat_cur_r = lat_cur*coef
  for j in range(2,nb_lon+1,2):
    k = i*nb_lon + j-1
    # We need the longitude to be in [0,360]
    lon_cur = j*dlon
    lon_cur_r = lon_cur*coef
    if (lat_cur == 0 and lon_cur == 360):
        lat_prev_r = lat_cur_r
        lon_prev_r = lon_cur_r
    else:
      psi = get_psi(lat_cur_r,lon_cur_r,U_0)
      
      lat_prev_r = get_lat_prev(t2,t1,psi,lat_cur_r,lon_cur_r,U_0)
      sign = get_sign(t2,t1,lat_cur_r,lon_cur_r,psi,U_0)
      lon_prev_r = get_lon_prev(sign,lat_prev_r,lat_cur_r,lon_cur_r,U_0)
    
    lat_prev = (lat_prev_r*180/np.pi)
    lon_prev = (lon_prev_r*180/np.pi)
    
    i_prev = int((90.-lat_prev)/dlat)
    j_prev = int(lon_prev/dlon)
    k_prev = int((i_prev-1)*nb_lon + j_prev)

    if (k_prev > nb_points-1):
      print i_prev,j_prev,i,j,k_prev
    lat[k] = lat_cur
    lon[k] = lon_cur
    if (lon_cur > 180):
      lon[k] = lon_cur - 360

    q_new[k] = q[k_prev]
    k = k + 1

###################################################
# Write data for comparison
###################################################
name = "../data/data_python/lol_"+str(nb_iter)+".dat"
f = open(name,'w')
for k in range(nb_points):
  s = "%f %f %f \n" % (lat[k],lon[k],q_new[k])
  f.write(s)
f.close()

###################################################
# Plot
###################################################
map = Basemap(projection='mill',llcrnrlat=-90,urcrnrlat=90,\
                    llcrnrlon=-180,urcrnrlon=180)
map.drawcoastlines()
map.drawmapboundary()

x,y = map(lon,lat)
map.scatter(x,y, s=10, c=q_new, marker='o',edgecolors='none')
map.colorbar()
name_fig = 'bidim_num_'+str(nb_iter)+'.jpg'
p.savefig(name_fig)
sub.call(["eog", name_fig])
