import numpy as np
import math
import sys

###########################################
def init_q(nb_lat,nb_lon):
  nb_points = nb_lat*nb_lon
  q = np.zeros(nb_points)
  # Init : the concentration goes from 0 to nb_lon/2 to 0
  # on [0,360]
  k = 0
  for i in range(nb_lat):
    for j in range(1,21):
      q[k] = j
      k = k + 1
    k = k + nb_lon-20
#    print k
  return q

def init_u(nb_lat,nb_lon,dlat):
  u = np.zeros(nb_lat)
  coef = np.pi/180
  inv_norm = 1./180
  for i in range(nb_lat):
    lat_cur = (89.5 - i*dlat)
    lat_cur_r = lat_cur*coef
    u[i] = inv_norm*np.cos(lat_cur_r)
  return u

def get_nb_mesh(i,nb_lat,nb_lat2):
  if (i > nb_lat2):
    i2 = nb_lat + 1 - i
  else:
    i2 = i
  nb_mesh = 3*(2*i2-1)
  return nb_mesh


###########################################

U_0 = 1./174
t1 = 0
dt = 120.
if (len(sys.argv) > 1):
  nb_iter = int(sys.argv[1])
else:
  print "Missing the number of iterations !"
  print "Usage : "+str(sys.argv[0])+" NB_ITER"
  exit(2)


nb_lat = 180
nb_lon = 537#360

dlat = 120./179
dlon = 1.
k = 0
t2 = nb_iter*dt
coef = np.pi/180

q = init_q(nb_lat,nb_lon)
u = init_u(nb_lat,nb_lon,dlat)

name = "../plot/data_python/"+str(nb_iter)+".dat"
f = open(name,'w')
nb_lat2 = 90

for i in range(nb_lat):
  lat_cur = (89.5 - i*dlat)
  lat_cur_r = lat_cur*coef
  # Wind
  coslat = np.cos(lat_cur_r)
  u_cur = u[i]
  for j in range(nb_lon+1):
    # We need the longitude to be in [0,360]
    lon_cur = (j+0.5)*dlon
      
    lon_prev = lon_cur - u_cur*t2/coslat
    j_prev = int(lon_prev/dlon)
    k_prev = int(i*nb_lon + j_prev)

    s = "%f %f %f \n" % (lat_cur,lon_cur,q[k_prev])
    f.write(s)

f.close()
