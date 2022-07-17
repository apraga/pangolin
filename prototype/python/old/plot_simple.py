import numpy as np
import sys 
#import pylab as p 
import matplotlib.pyplot as p

import subprocess

def file_len(fname):
  p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
    stderr=subprocess.PIPE)
  result, err = p.communicate()
  if p.returncode != 0:
    raise IOError(err)
  return int(result.strip().split()[0])

if (len(sys.argv) > 1):
  nb_iter = int(sys.argv[1])
else:
  print "Missing the number of iterations !"
  print "Usage : "+str(sys.argv[0])+" NB_ITER"
  exit(2)

name = "../plot/data_python/"+str(nb_iter)+".dat"
nb_lines = file_len(name)

lat = np.zeros(nb_lines,float)
lon = np.zeros(nb_lines,float)
q = np.zeros(nb_lines,float)
k = 0
f = open(name,'r')
for line in f.readlines():
  cur = line.split()
  lat[k] = float(cur[0])
  lon[k] = float(cur[1])
  q[k] = float(cur[2])
  k = k + 1
f.close()

###############################################################################
# Plotting
###############################################################################
p.subplot(111)
p.scatter(lon,lat, s=10, c=q, marker='o',edgecolors='none')

#p.xlim(0, 60)
#.ylim(0, 90)
p.colorbar()
#p.show()
p.savefig('compare.jpg')
#xlabel(r"Result", fontsize = 12)
#ylabel(r"Prediction", fontsize = 12)
