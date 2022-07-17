#-------------------------------------------------------------------------------
# Class for the reduced grid
#-------------------------------------------------------------------------------

import math

class Grid():

  # Arguments : Grid( [dlat [nb_lat] ] )
  def __init__(self, *args, **kwargs):
    self.dlat = 1. if (len(args) < 1) else args[0]
    self.nb_lat = int(180/self.dlat) if (len(args) < 2) else args[1]
    # Number of points on the north hemispher
    # Warning : can be different from nb_lat/2
    self.nb_lat2 = int(90./self.dlat)
    if (self.nb_lat2 > self.nb_lat):
      raise SystemExit("Not enough points on the North Hemisphere. We need at"\
      " least "+str(self.nb_lat2)+" points.")

  #-----------------------------------------------------------------------------
  # Size and number of cells
  #-----------------------------------------------------------------------------
  # Total number of points on the whole sphere
  def get_size(self,i=-1):
    return 3*self.get_size_zone(i)

  def get_size_zone(self,i=-1):
    if (i < 0):
      i = self.nb_lat

    if (i*self.dlat > 90.):
      total = 2*self.nb_lat2**2 - (180 - i)**2 
    else:
      total = i**2
    return total 

  # Get the number of cells at latitude i in all zones
  def get_nb_cells(self,i):
    return 3*self.get_nb_cells_zone(i)

  def get_nb_cells_zone(self,i):
    if (i > self.nb_lat2):
      i = 180 - i + 1
    return (2*i-1)

  #-----------------------------------------------------------------------------
  # Neighbours
  # Convention : i,j are positive
  #-----------------------------------------------------------------------------
  # Compute global position from local
  def local_to_global(self,i,j):
    return self.get_size(i-1) + j

  # Find adjacent neighbours in the (lat,lon) format
  def find_neighbours_adj(self,i,j,neighbours):
    nb_cells = self.get_nb_cells(i)
    j_prev = j-1 if (j > 1)  else nb_cells
    j_next = j+1 if (j < nb_cells) else 1
    # Does not matter here if j <= 0 or j >= nb_cells 
    neighbours.extend([[i,j_prev],[i,j_next]])

  # Find neighbours on the previous latitude in the (lat,lon) format
  def find_neighbours_prev(self,i,j,neighbours,nb_cells_zone):
    if (i > 1):
      ratio = self.get_nb_cells_zone(i-1)/float(nb_cells_zone)
      j1 = int(math.floor(ratio*(j-1)+1))
      j2 = int(math.ceil(ratio*j))
      neighbours.extend([[i-1,j1]])
      if (j1 < j2):
        neighbours.extend([[i-1,j2]])

  # Find neighbours on the next latitude in the (lat,lon) format
  def find_neighbours_next(self,i,j,neighbours,nb_cells_zone):
    if (i < self.nb_lat):
      ratio = self.get_nb_cells_zone(i+1)/float(nb_cells_zone)
      j1 = int(math.floor(ratio*(j-1)+1))
      j2 = int(math.ceil(ratio*j))
      neighbours.extend([[i+1,j1]])
      if (j1+1 < j2):
        neighbours.extend([[i+1,j1+1]])
      neighbours.extend([[i+1,j2]])

  # Find the neigbours on a zone !
  # Returns a list with pairs in the (lat,lon) format 
  # i,j > 0
  def find_neighbours_zone(self,i,j):
    if (j == 0):
      SystemExit('Need j > 0 for neighbours search')
    nb_cells_zone = self.get_nb_cells_zone(i)
    neighbours = []
    self.find_neighbours_adj(i,j,neighbours)
    self.find_neighbours_prev(i,j,neighbours,nb_cells_zone)
    self.find_neighbours_next(i,j,neighbours,nb_cells_zone)
    return neighbours

  # Finds the indices of the neighbours array 
  # !! If a neighbours is outside a zone, we put -1 as its index
  def position_to_indice(self,neighbours):
    out = []
    for x in neighbours:
      # k = (i-1)^2 + j
      line = self.get_size_zone(x[0]-1)
      if (x[1] < self.get_nb_cells_zone(x[0])+1):
        out.append(line+x[1]-1)
      else:
        out.append(-1)
    return out

## Compute local position from global with k in [0,nt-1]
  ## i in [1,n] and j in [0,nb_cells -1]
  ## We want k = (i-1)^2 + j 
  #def global_to_local(self,k):
  #  # If i = 1
  #  if (k < self.get_nb_cells(1)):
  #    return 1,k

  #  i = int(np.sqrt(k))
  #  j = k - (i-1)**2
  #  # Performs a check
  #  while (j > self.get_nb_cells(i)):
  #    i = i + 1
  #    j = k - (i-1)**2
  #  return i,j



