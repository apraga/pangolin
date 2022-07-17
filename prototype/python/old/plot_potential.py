from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.arange(0, 2*np.pi, 0.1)

Y = np.arange(-np.pi*0.5, np.pi*0.5, 0.1)
X, Y = np.meshgrid(X, Y)
#Z = np.cos(X*0.5)**2*np.cos(Y)**2
#Z = 10/5*np.sin(X)**2*np.cos(Y)**2 - 2*np.pi/5*np.sin(Y)
Z =  - 2*np.pi/5*np.sin(Y)
#Z = np.sin(X)**2*np.cos(Y)**2 
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
            linewidth=0, antialiased=False)

#cset = ax.contour(X, Y, Z, zdir='z')#, offset=-10)

ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

#plt.savefig('potential.jpg')
