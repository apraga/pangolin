import numpy as np
import math
import sys
import matplotlib.pyplot as p

# Zonal advection for lat/lon grid, then projected on
# reduced grid

###########################################
def init_q(nb_lat,nb_lon):
  nb_points = nb_lat*nb_lon
  q = np.zeros(nb_points)
  # Init : the concentration goes from 0 to nb_lon/2 to 0
  # on [0,360]
  k = 0
  nb_init = 20
  for i in range(nb_lat):
    for j in range(1,nb_init+1):
      q[k] = 80+j
      k = k + 1
    k = k + nb_lon-nb_init
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
nb_lon = 360#537#360

dlat = 360./nb_lon
dlon = 1.
k = 0
t2 = nb_iter*dt
coef = np.pi/180

q = init_q(nb_lat,nb_lon)
u = init_u(nb_lat,nb_lon,dlat)

nb_lat2 = 90
n = get_nb_mesh(90,nb_lat,nb_lat2)/3
nb_points = 6*nb_lat2*nb_lat2
#nb_points = nb_lat*nb_lon
lat = np.zeros(nb_points,float)
lon = np.zeros(nb_points,float)
q_new = np.zeros(nb_points,float)
k = 0

for i in range(nb_lat):
  lat_cur = (89.5 - i*dlat)
  lat_cur_r = lat_cur*coef
  # Wind
  coslat = np.cos(lat_cur_r)
  u_cur = u[i]
  nb_mesh = get_nb_mesh(i+1,nb_lat,nb_lat2)
  dlon2 = 360./nb_mesh
  for j in range(nb_mesh):
    lon_cur = (j+0.5)*dlon2
    lon_prev = lon_cur - u_cur*t2/coslat
    j_prev = int(lon_prev/dlon)
    k_prev = int(i*nb_lon + j_prev)

    lat[k] = lat_cur
    lon[k] = lon_cur
    q_new[k] = q[k_prev]
    k = k + 1
###################################################
# Write data for comparison
###################################################
name = "../plot/data_python/"+str(nb_iter)+".dat"
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
map.scatter(x,y, s=10, c=q, marker='o',edgecolors='none')

map.colorbar()

p.scatter(lon,lat, s=20, c=q_new, marker='o',edgecolors='none')

#p.xlim(0, 60)
#.ylim(0, 90)
p.colorbar()
#p.show()
p.savefig('fig.jpg')
