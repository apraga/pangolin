import numpy as np

# Visu mesh
nb_lat = 90
f = open('mesh.txt','w')
for i in range(1,nb_lat+1,10):
  dlon = 120.0/(2*i-1)
  for j in range(0,3*(2*i-1),10):
    lat = (90 - i)*np.pi/180
    lon = (j*dlon)*np.pi/180
    x = np.cos(lat)*np.cos(lon)
    y = np.cos(lat)*np.sin(lon)
    z = np.sin(lat)
    s = "%f %f %f \n" % (x,y,z)
    f.write(s)
  #f.write("\n")
f.close()

