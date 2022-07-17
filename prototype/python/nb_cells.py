###############################################################################
# Plot true and approximate number of cells for our reduced grid
###############################################################################
import numpy as np
import pylab as p

# Approximate cos
def approx(x,n_ref):
  #i = int(90.5+x);
  i = int(90.5-x);
  if (x > 90./3):
    res = 3*(2*i-1)
  else:
    res = n_ref
  return res


n = 90
#x = np.arange(-89.5, -0.5, 1.)
x = np.arange(89.5, 0.5, -1.)
coef = np.pi/180
n_ref = 360
dlambda_ref = 2*np.pi/n_ref
# Constant latitude step
dphi = dlambda_ref
kappa = 2*n_ref*(np.sin(dphi*0.5)/np.sin(dphi))
nb = np.array([kappa*np.cos(i*coef) for i in x])
nb2 = np.array([3*(2*i-1) for i in range(1,len(x)+1)])
nb3 = np.array([approx(i,n_ref) for i in x])
p.plot(x, nb,label='Exact', linewidth=2.0)
p.plot(x, nb3,label='Tronc',linestyle='-', linewidth=2.0)
p.plot(x, nb2,label='Approx',linestyle='--', linewidth=2.0)

xl = p.xlabel('Latitude $\Phi$')
yl = p.ylabel('n')
ttl = p.title('Number of cells with $\Delta \lambda=360$')
p.legend()
#p.ylim([0,1.2*n_ref])

p.savefig('nb_cells.png')
p.show()
