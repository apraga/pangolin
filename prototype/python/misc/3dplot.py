from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.arange(-180, 180, 10.)
Y = np.arange(-90.5, 90.5, 1.)

X, Y = np.meshgrid(X, Y)
coef = np.pi/180
#Z = np.sin(X*coef*0.5)**2*np.sin(Y*coef)**2
Z = 1-np.cos(X*coef)**2*np.cos(Y*coef)**2
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
        linewidth=0, antialiased=False)
ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
