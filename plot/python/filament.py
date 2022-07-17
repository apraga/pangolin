#-------------------------------------------------------------------------------
# Compute filament preservation for the reduced grid
# We assume the grid is stored as a single array ordered by latitudes
#-------------------------------------------------------------------------------
import error
import sys
import numpy as np
#import pylab as plt

# Returns the spherical area = sum(A(k)) where q > tau on A(k)
def spherical_area_old(q,tau,nb_lat,nb_lat2):
  k = 0
  area = 0
  dlat = np.pi/nb_lat
  for i in range(1,nb_lat+1):
    nb_cells = error.get_nb_cells(i,nb_lat2,nb_lat)
    lat = 0.5*np.pi - i*dlat
    dlon = 2*np.pi/nb_cells
    dS = dlon*(np.sin(lat+dlat)-np.sin(lat))
    for i in range(nb_cells):
      cur = q[k]
      if (cur >= tau):
        #area = area + cur*dS
        area = area + dS
      k = k + 1
  return area

# Sum of all areas where q >= tau
def spherical_area(q, lat, lon, tau):
  k_start = 0
  k = 0
  k_max = len(lat)-1
  area = 0
  coef = np.pi/180
  while (True):
    if (k > k_max):
      break
    # We cut the array into lines of constant latitude
    dlon = lon[k_start+1] - lon[k_start]
    nb_cells = error.get_nb_cells_dlon(dlon)

    k_next = k_start + nb_cells
    if (k_next < k_max):
      dlat = abs(lat[k_next] - lat[k_start])
    # We assume the radius is 1
    cur_lat = lat[k_start]*coef
    dS = 2*np.cos(cur_lat)*np.sin(0.5*dlat*coef)*dlon*coef
    for i in range(nb_cells):
      if (q[k] >= tau): area = area + dS
      k = k + 1
    k_start = k
    # Normalize
  return area



# Compute the diagnostic
def diagnostic_tau(q,q0, lat, lon, tau):
  l_f = 0.
  area0 = spherical_area(q0, lat, lon, tau)
  precision = 1e-10
  # Compute the filament preservation diagnostic
  if (area0 > precision):
    area = spherical_area(q, lat, lon, tau)
    l_f = 100.0*area/area0
  return l_f

# Return all diagnostics for all tau
# Needs data at t=0 and t=T/2
def diagnostics(file0, fileT2, degree, tau):
  # q(t=0)
  q0, lat0, lon0 = error.read_data(file0)
  # q(t)
  q, lat,lon = error.read_data(fileT2)
  
  error.check_grids(lat0,lon0,lat,lon)
  
  n = error.file_len(file0)
  # Because nb_points = 6*nb_lat^2 on an hemisphere
  #nb_lat2 = int(np.sqrt(n/6))
  #nb_lat = 2*nb_lat2
  l_f = np.zeros(len(tau))
  i = 0
  for i in range(len(tau)):
    #l_f[i] = diagnostic_tau(q,q0,tau[i],nb_lat,nb_lat2)
    l_f[i] = diagnostic_tau(q, q0, lat, lon, tau[i])
  return l_f

# Plot one diagonstic for a given resolution 
# resol is a string for the resolution
# style is a sequence for dash type
def plot_diagnostic(file0, fileT2, tau, labl, style, col):
  l_f1 = diagnostics(file0, fileT2, '1.5',tau)

  args = [('dashes',style), ('color',col) ]
  if (labl is not None): args.append(('label',labl))
  # Convert list to dictionary
  args = dict(args)
  plt.plot(tau,l_f1, '--', **args)

#-------------------------------------------------------------------------------

#folder = '/wkdir/pae2/praga/filaments_unlimited/'
#
#tau = np.arange(0.1,1.1,0.1)
#legend = ['1.5', '0.75', '0.375']
#style = [(8,6), (2,3), (12, 3)]
#labels = ['$\Delta \lambda='+i+'\degree$' for i in legend]
#col = ['b','g','r']
#k = 0
#
#for i in [40, 80, 160]:
#  file0 = folder+"ratio_"+str(i)+"lat_0.dat"
#  fileT2 = folder+"ratio_"+str(i)+"lat_T2.dat"
#
#  if (i == 80): labl = "unlimited"
#  else: labl = None
#  plot_diagnostic(file0, fileT2, tau, labl, (8,5), col[k])
#  k += 1
#
#folder = '/wkdir/pae2/praga/filaments_limited/'
#k = 0
#for i in [40, 80, 160]:
#  file0 = folder+"ratio_"+str(i)+"lat_0.dat"
#  fileT2 = folder+"ratio_"+str(i)+"lat_T2.dat"
#
#  plot_diagnostic(file0, fileT2, tau, labels[k], (None, None), col[k])
#  k += 1
#
#
#plt.xlabel('$\\tau$')
#plt.ylabel('$l_f$')
## Same limits as Lauritzen
#plt.ylim([0,140])
#limit = [100. for i in range(len(tau))]
#plt.plot(tau,limit,'black')
#plt.legend(loc='upper right')
#plt.title('Filaments preservation (un- and limited)')
#plt.savefig('filaments.pdf')
#
#plt.show()
