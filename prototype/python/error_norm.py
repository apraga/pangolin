###############################################################################
# Compute error between analytical and numerical solution
# The normalized errors are : 2, oo, max
# Analytical and numerical solutions are defined in two files 
# so we just read the files
###############################################################################

import sys
import os
import os.path
import lib_error


###############################################################################
# Main
###############################################################################
nb_args = 3
#folder = '../data/data_standard_tests/'
folder = '../data/data/'
#folder = '../data/data_standard_tests_int/'
if (len(sys.argv) == nb_args+1):
  error = sys.argv[1]
  anal = sys.argv[2]
  num = sys.argv[3]
else:
  print "Usage : "+str(sys.argv[0])+" ERROR ANALYTICAL NUMERICAL"
  print "ERROR = 2, oo, quad"
  print "ANALYTICAL, NUMERICAL : files contained in "+folder
  exit(2)

# Plot
anal = folder+anal
num = folder+num
res = lib_error.get_error(anal,num,error);

print res
