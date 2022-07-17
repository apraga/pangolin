#from pylab import *
import pylab as p
import numpy as np
import subprocess as sub

def file_len(fname):
  p = sub.Popen(['wc', '-l', fname], stdout=sub.PIPE, 
    stderr=sub.PIPE)
  result, err = p.communicate()
  if p.returncode != 0:
    raise IOError(err)
  return int(result.strip().split()[0])

def dirac(x):
  y = np.zeros(len(x))
  i = 0
  for k in x:
    if (k > 120 and k < 130):
      y[i] = 1.
    i = i + 1
  return y

name = "zonal_349.dat"
nb_lines = file_len(name)

lat = np.zeros(nb_lines,float)
lon = np.zeros(nb_lines,float)
q = np.zeros(nb_lines,float)
f = open(name,'r')
k = 0
for line in f.readlines():
  cur = line.split()
  lat[k] = float(cur[0])
  lon[k] = float(cur[1])
  q[k] = float(cur[2])
  k = k + 1
f.close()

# Plot
x = np.arange(min(lon),max(lon),0.1)
y = dirac(x)
p.plot(x,y,color='red')

p.scatter(lon, q, s=20, c='b', marker='o',
    faceted=False)

p.xlabel(r"Longitude", fontsize = 12)
p.ylabel(r"Ratio", fontsize = 12)
p.title(r"Profil de concentration (zonal, 349iter)", fontsize = 12)
p.show()
name_fig = 'profil_zonal.jpg'
p.savefig(name_fig)

