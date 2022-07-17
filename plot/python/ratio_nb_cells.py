from pylab import *
import numpy as np

def ratio(dlat):
  coef = np.pi/180.
  return (1 - np.cos(dlat*coef))/np.sin(dlat*coef)*2*np.pi/(dlat*coef)

tmp = 20
dx = 0.5
x = np.arange(0.01,10,0.01)

y = [ratio(cur) for cur in x]
z = [np.pi for cur in x]

err = [(ratio(cur)-np.pi)/np.pi for cur in x]

fig = figure()
ax = fig.add_subplot(1, 1, 1)
#ax.plot(x,y)
#ax.plot(x,z)
ax.plot(x,err)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('$\Delta\phi$', fontsize=20)
ax.set_ylabel('Relative error')
title("Number of cells at the pole $\\rightarrow \pi$", fontsize=20)
savefig("ratio_nb_cells.png")
show()
