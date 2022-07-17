import numpy as np
import math
import sys

###########################################
def init_q(nb_points):
  # Init
  q = np.zeros(nb_points)
  lat = np.zeros(nb_points)
  lon = np.zeros(nb_points)
  k = 0
  name1 = "../plot/data/0.dat"
  f = open(name1,'r')
  for line in f.readlines():
    tmp = line.split()
    lat[k] = tmp[0]
    lon[k] = tmp[1]
    q[k] = tmp[2]
    k = k +1
  f.close()
  
  if (k != nb_points):
    print "Not the same number of points in 0.dat !"
 
  return q,lat,lon

def init_u(nb_lat,nb_lon,dlat):
  u = np.zeros(nb_lat)
  coef = np.pi/180
  inv_norm = 1./180
  for i in range(nb_lat):
    lat_cur = (89. - i*dlat)
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

if (len(sys.argv) > 1):
  nb_iter = int(sys.argv[1])
else:
  print "Missing the number of iterations !"
  print "Usage : "+str(sys.argv[0])+" NB_ITER"
  exit(2)

lon_prev = 1
U_0 = 1./174
t1 = 0
dt = 120.

nb_lat = 180#89#179
nb_lat2 = 90
nb_lon = 360
nb_points = 6*nb_lat2*nb_lat2

dlat = 1.
dlon = 1.
t2 = nb_iter*dt

coef = np.pi/180
q,lat,lon = init_q(nb_points)
u = init_u(nb_lat,nb_lon,dlat)

name = "../plot/data_python/"+str(nb_iter)+".dat"
f = open(name,'w')
nb_total = 0

#for i in range(1,nb_lat,5):
for i in range(1,nb_lat):
  lat_cur = (90.5 - i*dlat)
  lat_cur_r = lat_cur*coef
  nb_mesh = get_nb_mesh(i,nb_lat,nb_lat2)
  dlon = 360./nb_mesh
  # Wind
  coslat = np.cos(lat_cur_r)
  u_cur = u[i-1]

  for j in range(1,nb_mesh+1):
    lon_cur = (j-0.5)*dlon
    lon_prev = lon_cur - u_cur*t2/coslat
    if (lon_prev < 0):
      lon_prev = 360. + lon_prev
    j_prev = int(lon_prev/dlon)
    k_prev = nb_total + j_prev
    if (k_prev > nb_points-1):
     print k_prev,i,j_prev
    #print lon_prev,lon[k_prev],lon[k_prev-1]
    s = "%f %f %f \n" % (lat_cur,lon_cur,q[k_prev])
    f.write(s)

  nb_total = nb_total + nb_mesh

f.close()
