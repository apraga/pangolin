from pylab import *
import numpy as np

# Plot the difference between exact nb of cells and approximation
nb_lat2 = 80

approx = []
exact = []
dlat = 90./nb_lat2
coef = np.pi/180.
for i in range(1,nb_lat2+1):
  nb_cells = 2*i-1
  approx.append(nb_cells)

  lat = (i-0.5)*dlat
  res = 2*np.sin(0.5*dlat*coef)*np.sin(lat*coef)/(1 - np.cos(dlat*coef))
  exact.append(res)

lat = [90-i*dlat for i in range(1,nb_lat2+1)]
relative = [abs((x-x0)/x0) for x,x0 in zip(approx, exact)]
diff = [abs(x-x0) for x,x0 in zip(approx, exact)]
a=[x for x in relative if (x < 0.01) ]
pos = size(a)
print size(a), a, relative[pos], lat[pos]
print "relative error", max(relative), min(relative)
print "absolute error", max(diff), min(diff)
#plot(lat, diff)
plot(lat, relative)
pos = size(relative)-1
print relative[pos], approx[pos], exact[pos]
show()
