# Compute mean square error for each iteration and plot it
import sys
import os
import pylab as p
import os.path
import subprocess as sub
import numpy as np
from pylab import *

###############################################################################
# Get length of file
def file_len(fname):
  p = sub.Popen(['wc', '-l', fname], stdout=sub.PIPE, 
    stderr=sub.PIPE)
  result, err = p.communicate()
  if p.returncode != 0:
    raise IOError(err)
  return int(result.strip().split()[0])

# Compute the quadratic mean error
def get_error(file1,file2):
  nb_lines1 = file_len(file1)
  nb_lines2 = file_len(file2)
  if (nb_lines1 != nb_lines2): 
    print 'Not same number of lines for : '+file1+' and '+file2
    sys.exit()
  f1 = open(file1, 'r')
  f2 = open(file2, 'r')
  error = 0
  while True:
    line1 = f1.readline()
    line2 = f2.readline()
    if (not line2 or not line1): break

    cur1 = line1.split()
    lat1 = float(cur1[0])
    lon1 = float(cur1[1])
  
    cur2 = line2.split()
    lat2 = float(cur2[0])
    lon2 = float(cur2[1])

   # print cur1,cur2
    if (lat2 != lat1 or lon2 != lon1):
      print lat2,lat1,lon2,lon1,'Grid differs !'
      sys.exit()
    error = error + (float(cur2[2]) - float(cur1[2]))**2
  error = error/float(nb_lines1)
  return error
  f1.close()
  f2.close()


# Get number of iteration
def get_nb(name):
  size = len(name) - len('.dat')
  nb = int(name[0:size])
  return nb

def write_error(filename,version):
  rootdir='../data/data'+version+'/'
  pythondir='../data/data_python'+version+'/'
  
  # Search twin files
  twins = []
  for subdir, dirs, files in os.walk(rootdir):
    for file in files:
      file_python = pythondir+file
      if os.path.exists(file_python):
        twins.append(file)
  
  # Now get to work
  f = open(filename,'w')
  for file in twins:
    print file+'... ',
    file_python = pythondir+file
    x = get_nb(file)
    file= rootdir+file
    error = get_error(file,file_python)
    print ' '
    s = "%d %f \n" % (x,error)
    f.write(s)

  f.close()


###############################################################################
version = '_v1'
filename = 'error'+version+'.txt'
#write_error(filename,version)

# Plot
f = open(filename,'r')
k = 0
nb_lines = file_len(filename)
error = np.zeros(nb_lines)
x = np.zeros(nb_lines)
for line in f.readlines():
  cur = line.split()
  x[k] = cur[0]
  error[k] = cur[1]
  k = k + 1
f.close()
x.sort()
# Hack
error.sort()


# Second plot
version = '_v2'
filename = 'error'+version+'.txt'
f2 = open(filename,'r')
k = 0
nb_lines2 = file_len(filename)
error2 = np.zeros(nb_lines2)
x2 = np.zeros(nb_lines2)
for line in f2.readlines():
  cur = line.split()
  x2[k] = cur[0]
  error2[k] = cur[1]
  k = k + 1
f2.close()
x2.sort()
# Hack
error2.sort()


# Plot
#ax = p.subplot(111)
x2 = x2 * 43./80
p.plot(x,error,color='red',linewidth=2.5,label = 'C = 0.91')
p.plot(x2,error2,color='blue',linewidth=2.5,linestyle='--',label = 'C = 0.5')

p.xlabel("Iteration", fontsize = 12)
p.ylabel(r"Erreur", fontsize = 12)
p.title(r"Erreur quadratique moyenne", fontsize = 12)
p.legend(loc='lower right')
p.savefig('error_both.jpg')
p.show()

