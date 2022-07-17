#-------------------------------------------------------------------------------
# Compute error between analytical and numerical solution
# The normalized errors are : 2, oo, max
# refytical and numerical solutions are defined in two files 
# so we just read the files
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# Contains functions used for plotting errors
#-------------------------------------------------------------------------------

import sys
sys.path.append('../plot/python')
import subprocess as sub
#from lib_plot import file_len
import numpy as np

def check_grids(lat,lon,lat1,lon1):
  if (not np.array_equal(lat1,lat) or not np.array_equal(lon1,lon)):
    print "Files do not have the same grid"
    exit(2)

#-------------------------------------------------------------------------------
# Read a file on format (q, lat, lon) and return them
# It is faster to operate on an array
#-------------------------------------------------------------------------------
def read_data(filename):
  nb_lines = file_len(filename)
  f = open(filename, 'r')
  val1 = np.zeros(nb_lines,float)
  val2 = np.zeros(nb_lines,float)
  val3 = np.zeros(nb_lines,float)

  for k in range(nb_lines):
    line = f.readline()
    if (not line): break
    cur = line.split()
    val1[k] = float(cur[0])
    val2[k] = float(cur[1])
    val3[k] = float(cur[2])
  
  f.close()
  return val1,val2,val3

#-------------------------------------------------------------------------------
# Write a file on format (val1,val2,val3)
#-------------------------------------------------------------------------------
def write_data(filename,val1,val2,val3):
  nb_lines = len(val1)
  if (len(val2) != nb_lines or len(val3) != nb_lines):
    print "Cannot write data : arrays of different sizes" 
    exit(2)

  f = open(filename, 'w')
  for k in range(nb_lines):
    s = "%f %f %f \n" % (val1[k],val2[k],val3[k])
    f.write(s)
  f.close()


# Get the number of cells in a line at a given latitude from dlon
def get_nb_cells_dlon(dlon):
  nb_cells = int(360./dlon)
  remain = 360. - nb_cells*dlon
  error = 1.e-8
  if (dlon - remain < error):
    nb_cells = nb_cells + 1
  return nb_cells 

# Get the total number of cells in a line at latitude i
def get_nb_cells(i,nb_lat2,nb_lat):
  i2 = i 
  if (i > nb_lat2):
    i2 = nb_lat - i + 1
  return 3*(2*i2-1)

#-------------------------------------------------------------------------------
# Different types of errors
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Discret integral
# We assume all cells at the same latitude are consecutive and they have the 
# same delta lon
#-------------------------------------------------------------------------------
def spheric_integral(q,lat,lon):
  k_start = 0
  k = 0
  k_max = len(lat)-1
  result = 0
  coef = np.pi/180
  while (True):
    if (k > k_max):
      break
    # We cut the array into lines of constant latitude
    dlon = lon[k_start+1] - lon[k_start]
    nb_cells = get_nb_cells_dlon(dlon)

    k_next = k_start + nb_cells
    if (k_next < k_max):
      dlat = abs(lat[k_next] - lat[k_start])
    # We assume the radius is 1
    cur_lat = lat[k_start]*coef
    dS = 2*np.cos(cur_lat)*np.sin(0.5*dlat*coef)*dlon*coef
    for i in range(nb_cells):
      result = result + q[k]*dS
      k = k + 1
    k_start = k
    # Normalize
  return result/(4*np.pi)

#-------------------------------------------------------------------------------
# Error with 2-norm
#-------------------------------------------------------------------------------
def error_2(q_ref,q_num,lat,lon):
  num = (q_num - q_ref)**2
  denom = q_ref**2
  int1 = spheric_integral(num,lat,lon)
  int2 = spheric_integral(denom,lat,lon)
  return np.sqrt(int1/int2)

#-------------------------------------------------------------------------------
# Error with max norm
#-------------------------------------------------------------------------------
def error_oo(q_ref,q_num,lat,lon):
  diff = abs(q_num - q_ref)
  q_abs = abs(q_ref)
  return max(diff)/max(q_abs)

#-------------------------------------------------------------------------------
# Quadratic mean error (weighted)
#-------------------------------------------------------------------------------
def error_quad(q_ref,q_num,lat,lon):
  res = spheric_integral((q_num - q_ref)**2, lat, lon)
  return res


#-------------------------------------------------------------------------------
# Wrapper
#-------------------------------------------------------------------------------
def get_error(ref,num,err):
  q_ref, lat,lon = read_data(ref)
  q_num, lat1,lon1 = read_data(num)
  
  #check_grids(lat,lon,lat1,lon1)
  
  if (err == '2'):
    res = error_2(q_ref,q_num,lat,lon)
  elif (err == 'oo'):
    res = error_oo(q_ref,q_num,lat,lon)
  elif (err == 'quad'):
    res = error_quad(q_ref,q_num,lat,lon)
  else:
    print 'Error type is not valid', err
    res = -1;
  return res


#------------------------------------------------------------------------------
# Get number of lines in file
#------------------------------------------------------------------------------
def file_len(fname):
  p = sub.Popen(['wc', '-l', fname], stdout=sub.PIPE, 
    stderr=sub.PIPE)
  result, err = p.communicate()
  if p.returncode != 0:
    raise IOError(err)
  return int(result.strip().split()[0])
