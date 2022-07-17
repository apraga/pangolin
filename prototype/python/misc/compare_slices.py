import numpy as np
import sys 
import matplotlib.pyplot as p
from mpl_toolkits.basemap import Basemap
import subprocess as sub

#############################################################

def file_len(fname):
  p = sub.Popen(['wc', '-l', fname], stdout=sub.PIPE, 
    stderr=sub.PIPE)
  result, err = p.communicate()
  if p.returncode != 0:
    raise IOError(err)
  return int(result.strip().split()[0])

def get_nb_mesh(i,nb_lat,nb_lat2):
  if (i > nb_lat2):
    i2 = nb_lat + 1 - i
  else:
    i2 = i
  nb_mesh = 3*(2*i2-1)
  return nb_mesh

def plot_latitude(i,lat,nb_lines):
  n = get_nb_mesh(i,nb_lat,nb_lat2)
  error = np.zeros(n)
  x = np.zeros(n)
  lat_cur = 90.5-i

  j = 0
  for k in range(nb_lines):
    if (lat[k] == lat_cur and j < n): 
      error[j] = diff[k]
      x[j] = lon[k]
      j = j + 1
  p.plot(x,error)

  
#############################################################


if (len(sys.argv) > 1):
  nb_iter = int(sys.argv[1])
else:
  print "Missing the number of iterations !"
  print "Usage : "+str(sys.argv[0])+" NB_ITER"
  exit(2)

name1 = "../data/data/"+str(nb_iter)+".dat"
name2 = "../data/data_python/"+str(nb_iter)+".dat"
nb_lines1 = file_len(name1)
nb_lines2 = file_len(name2)
if (nb_lines1 != nb_lines2):
  print "Number of lines different !"
nb_lines = min(nb_lines1,nb_lines2)

lat = np.zeros(nb_lines,float)
lon = np.zeros(nb_lines,float)
diff = np.zeros(nb_lines,float)
k = 0
max_q = -10000
f1 = open(name1,'r')
for line in f1.readlines():
  cur = line.split()
  lat[k] = float(cur[0])
  lon[k] = float(cur[1])
  if (lon[k] > 180):
    lon[k] = lon[k] - 360.
  diff[k] = float(cur[2])
  max_q = max(max_q,float(cur[2]))
  k = k + 1
f1.close()

k = 0
f2 = open(name2,'r')
for line in f2.readlines():
  cur = line.split()
  diff[k] = float(cur[2]) - diff[k]
  max_q = max(max_q,float(cur[2]))
  k = k + 1
f2.close()
diff = diff/max_q

###############################################################################
# Plotting
###############################################################################
nb_lat = 180
nb_lat2 = 90
nb_lon = 360
# Plot some latitude bands
i = 135
lat_cur = 90 - i
plot_latitude(i,lat,nb_lines)

output = 'eps'
name_fig = 'error_'+str(nb_iter)+'_lat_'+str(i)+'.'+output
p.title('Erreur pour latitude = '+str(lat_cur)+'degres')
p.savefig(name_fig)
if (output == 'jpg'):
  sub.call(["eog", name_fig])
