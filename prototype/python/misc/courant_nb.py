import numpy as np
import pylab as p

def get_nb(i):
  if (i > 90):
    return 2*(181-i) - 1
  return 2*i - 1

x = np.arange(20.0, 161.0, 1.)
nb = np.array([get_nb(i) for i in x])
coef = np.pi/180
U_0 = 1./180.
s = 120./(U_0*nb*np.cos((90.5-x)*coef))
p.plot(x, s)
xl = p.xlabel('Colatitude')
yl = p.ylabel('dt')
ttl = p.title('Pas de temps maximal')
#p.savefig('dt_max.eps')
p.show()
