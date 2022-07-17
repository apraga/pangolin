import numpy as np

# Write vector field for gnuplot
nb_lat = 180
nb_lon = 360
f = open('../data/data_python/data_vectorfield.dat','w')
dlat = 1.
dlon = 1.
coef = np.pi/180

for i in range(0,nb_lat,2):
  lat_d = (90.-i*dlat)
  lat = lat_d*coef
  coslat = np.cos(lat)
  sinlat = np.sin(lat)
  for j in range(0,nb_lon,5):
    #lon_d = j*dlon-180.
    lon_d = j*dlon
    lon = lon_d*coef
    coslon = np.cos(lon*0.5)
    u = 2*coslat*sinlat*coslon*coslon*5
    v = -coslat*coslon*np.sin(lon*0.5)*5
    s = "%f %f %f %f \n" % (lon_d,lat_d,u,v)
    f.write(s)
f.close()


