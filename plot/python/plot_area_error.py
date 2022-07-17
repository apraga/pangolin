from pylab import *
import numpy as np

# Plot the difference between exact area formulation and approximation
nb_lat2 = 80

approx = []
exact = []
dlat = 90./nb_lat2
coef = np.pi/180.
for i in range(1,nb_lat2+1):
  nb_cells = 2*i-1
  dlon = 120./nb_cells
  lat = 90.-(i-0.5)*dlat
  res = dlon*dlat*np.cos(lat*coef)
  approx.append(res)

  lat_n = 90.-(i-1)*dlat
  lat_s = 90.-i*dlat
  res = dlon*(np.sin(lat_n*coef) - np.sin(lat_s*coef))/coef
  exact.append(res)

lat = [90-i*dlat for i in range(1,nb_lat2+1)]
relative = [abs((x-x0)/x0) for x,x0 in zip(approx, exact)]
diff = [abs(x-x0) for x,x0 in zip(approx, exact)]
print "relative error", max(relative), min(relative)
print "absolute error", max(diff), min(diff)
plot(lat, diff)
show()

A_true = 360*(np.sin(90*coef) - np.sin((90-dlat)*coef))/coef 
print "Pole error : absolute", abs(exact[0] - A_true)
print "Pole error : relative", abs(exact[0] - A_true)/A_true
print A_true, exact[0]
