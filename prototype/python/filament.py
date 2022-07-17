################################################################################
# Compute filament preservation for the reduced grid
# We assume the grid is stored as a single array ordered by latitudes
################################################################################
import lib_error
import lib_plot
import sys
import numpy as np
import pylab as plt

# Returns the spherical area = sum(A(k)) where q > tau on A(k)
def spherical_area(q,tau,nb_lat,nb_lat2):
  k = 0
  area = 0
  dlat = np.pi/nb_lat
  for i in range(1,nb_lat+1):
    nb_cells = lib_error.get_nb_cells(i,nb_lat2,nb_lat)
    lat = 0.5*np.pi - i*dlat
    dlon = 2*np.pi/nb_cells
    dS = dlon*(np.sin(lat+dlat)-np.sin(lat))
    for i in range(nb_cells):
      cur = q[k]
      if (cur >= tau):
        area = area + cur*dS
      k = k + 1
  return area

# Compute the diagnostic
def diagnostic(q,q0,tau,nb_lat,nb_lat2):
  l_f = 0.
  area0 = spherical_area(q0,tau,nb_lat,nb_lat2)
  precision = 1e-10
  # Compute the filament preservation diagnostic
  if (area0 > precision):
    area = spherical_area(q,tau,nb_lat,nb_lat2)
    l_f = 100.0*area/area0
  return l_f

# Return all diagnostics for all tau
def main(degree,folder,tau):
  prefix = 'filament_'+degree+'degree_'
  filename = folder+prefix+"0.dat"
  # q(t=0)
  lat0,lon0,q0 = lib_error.read_data(filename)
  # q(t)
  filename = folder+prefix+n_iter+".dat"
  lat,lon,q = lib_error.read_data(filename)
  
  lib_error.check_grids(lat0,lon0,lat,lon)
  
  n = lib_plot.file_len(filename)
  # Because nb_points = 6*nb_lat^2 on an hemisphere
  nb_lat2 = int(np.sqrt(n/6))
  nb_lat = 2*nb_lat2
  l_f = np.zeros(len(tau))
  i = 0
  for i in range(len(tau)):
    l_f[i] = diagnostic(q,q0,tau[i],nb_lat,nb_lat2)
  return l_f

################################################################################

nb_args = 1;
if (len(sys.argv) == nb_args+1):
  n_iter = sys.argv[1]
else:
  print "Usage : "+str(sys.argv[0])+" NB_ITER"
  exit(2)


folder = '../data/data_standard_tests/'
tau = np.arange(0.1,1.1,0.1)
l_f1 = main('1.0',folder,tau)
l_f2 = main('1.5',folder,tau)
limit = [100. for i in range(len(tau))]
plt.plot(tau,l_f1,label='$\Delta \lambda=1.0\degree$')
plt.plot(tau,l_f2,label='$\Delta \lambda=1.5\degree$')
plt.plot(tau,limit)
plt.legend(loc='upper right')
plt.title('Filaments diagnostics')
plt.savefig('filaments.png')


plt.show()
