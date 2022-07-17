import numpy

def split_single(step,list_middle,f):
  i = 1
  middle = 0.5/step
  while (middle < 1):
    s = "%f %f %d %d \n" % (middle,1/step,1,i)
    #f.write(s)
    middle += 1/step
    i += 1


def split(k,f):
  step_up = float(2*k-1)
  step_down = float(2*k+1)
  i = 0
  j = 0
  cur = 0
  prev = 0
  while (cur < 1):
    prev = cur
    next_up = (i+1)/step_up
    next_down = (j+1)/step_down
    i_prev = i
    j_prev = j
    if (next_up < next_down):
      i += 1
      cur = next_up
    else:
      j += 1
      cur = next_down
    middle = 0.5*(prev + cur) 
    step = cur -prev
    s = "%f %f %d %d \n" % (middle,step,i_prev+1,j_prev+1)
    #f.write(s)

###############
nb_lat = 90
list_middle = numpy.zeros((nb_lat,2*nb_lat - 1))
list_step = numpy.zeros((nb_lat,2*nb_lat - 1))
list_index_up = numpy.zeros((nb_lat,2*nb_lat - 1))
list_index_down = numpy.zeros((nb_lat,2*nb_lat - 1))
f = open('grid.txt','w')
split_single(3.0,list_middle,f)

for k in range(2,91):
 split(k,f)


