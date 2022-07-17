import lib_error
import numpy as np
import pylab as p
import sys

folder = '../data/data_standard_tests/'
filename = 'plotting_error.txt'
#readdata = True
readdata = False
if (readdata):
  print "Reading file instead of computing again"
  try :
    x,y1,y2 = lib_error.read_data(filename)
  except IOError as e:
    print "File does not exist."
 
else:
  x = np.array([120,240,360,480,720])
  y1 = np.zeros(len(x))
  y2 = np.zeros(len(x))
  
  init = "gaussian"
  k = 0
  for i in x:
    anal = folder+init+"_"+str(i)+"_0.dat"
    num = folder+init+"_"+str(i)+"_T.dat"
  
    y1[k] = lib_error.get_error(anal,num,'2')
    y2[k] = lib_error.get_error(anal,num,'oo')
    k = k +1

  lib_error.write_data(filename,x,y1,y2)

# Delta lambda
x2 = 360./x
# Reorder
x3 = [i for i in range(len(x))]
a = np.log(0.5)/2.25
y3 = [0.7*np.exp(a*i) for i in x3]

# Real values for x
p.xticks(x3,x2)

p.gca().set_yscale('log')
p.plot(x3,y1,label='2-norm')
p.plot(x3,y2,label='oo-norm')
p.plot(x3,y3,label='first-order')
p.legend(bbox_to_anchor=(1., 1.))
p.ylabel('Error')
p.xlabel('Resolution (degree)')
p.title('Convergence plot')
p.savefig('convergence_plot.png')
p.show()
