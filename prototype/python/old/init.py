import numpy

def concentration(nb_lat,nb_lon):
  C = numpy.zeros((nb_lat,nb_lon))
  for i in range(nb_lat):
    for j in range(5):
      C[i,j] = 1
  return C

# For all lat, compute the step and position
def grid(nb_lat):
  step = numpy.zeros((nb_lat))
  pos = numpy.zeros((nb_lat,2*nb_lat-1))
  nb = 2*nb_lat - 1
  for i in range(1,nb_lat+1):
    tmp = 1./(2*i - 1)
    step[i-1] = tmp
    for j in range(2*i-1):
      pos[i-1,j] = (j + 0.5)*tmp
  return (step,pos)

## Returns the middle of 2 nodes, the step, and the number of the mesh (previous
## and next)
#def nodes(nb_lat):
#  nb = 4*nb_lat - 1
#  prev = numpy.zeros(2*nb_lat-1)
#  next = numpy.zeros(2*nb_lat+1)
#  node_pos = numpy.zeros(4*nb_lat)
#
#  node_prev = numpy.zeros((nb_lat,4*nb_lat-1))
#  node_next = numpy.zeros((nb_lat,4*nb_lat-1))
#  node_step = numpy.zeros((nb_lat,4*nb_lat-1))
#  node_middle = numpy.zeros((nb_lat,4*nb_lat-1))
#  for i in range(1,nb_lat+1):
#    # nodes position according to latitude
#    for j in range(2*i-1):
#      prev[j] = j/(2*i-1)
#    for j in range(2*i+1):
#      next[j] = j/(2*i+1)
#    
#    node_pos[0] = 0
#
#    # Compute a node, each at a time
#    for j in range(4*i-1):
#      min1 = 1
#      min2 = 1
#      for k in range(1,2*i-1):
#        min1 = min(min1,prev[k])
#      for k in range(1,2*i+1):
#        min2 = min(min2,next[k])
#      space = min(min1,min2)
#      node_pos[j] = node_pos[j-1] + space
#
#      for k in range(2*i-1):
#        if (prev[k] > 0.99999999999):
#          prev[k] -= space
#        if (prev[k] < 1e-20):
#          prev[k] = 1
#
#      for k in range(2*i+1):
#        if (next[k] > 0.99999999999):
#          next[k] -= space
#        if (next[k] < 1e-20):
#          next[k] = 1
#
#    node_pos[4*i-1] = 1
#    for j in range(1,4*i):
#      node_step[i-1,j] = node_pos[j] - node_pos[j-1]
#      node_middle[i-1,j]= 0.5*(node_pos[j] + node_pos[j-1])
#
#    # indices
#    for j in range(4*i-1):
#      node_prev[i-1,j] = numpy.floor(node_middle[i-1,j]*(2*i-1))+1
#      node_next[i-1,j] = numpy.floor(node_middle[i-1,j]*(2*i+1))+1
#
#  return (node_middle,node_step,node_prev,node_next)
