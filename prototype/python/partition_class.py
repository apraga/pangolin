#-------------------------------------------------------------------------------
# Class for the partition
# Strategy : cuts in same width bands
#-------------------------------------------------------------------------------

import numpy as np 
import os
import random
import matplotlib.pyplot as plt

class Partition:

  def __init__(self, *args, **kwargs):
    if (len(args) < 2): 
      raise SystemExit('Partition needs a grid and number of partitions !')
    else:
      self.grid = args[0]
      # We partition only a zone (there are 3)
      self.size = self.grid.get_size_zone()
      self.part = np.zeros(self.size)
      self.nb_part = args[1]
      self.strategy = 0
      self.mp_dir = '../metapost/' 
      self.init_colors()

  # Display a part of the array
  def display(self,i1,i2):
    for i in range(i1,i2+1):
      print self.part[i],

  #-----------------------------------------------------------------------------
  # Partitioning one zone
  #-----------------------------------------------------------------------------
  def split_to_squares(self,square_side):
    part_cur = 0
  
    # First band is a single partition (easy)
    for i in range(self.grid.get_size_zone(square_side)):
      self.part[i] = part_cur
  
    nb_regions = 3
    i = square_side + 1
    part_cur = 1
    cur = square_side**2
    # Others bands are split in regions
    for k in range(1,self.nb_part2):
      mid = (nb_regions-1)/2
  
      # We examine each line and assing the cells to a region
      for line in range(square_side):
        nb_cells = self.grid.get_nb_cells_zone(i)
        cur = self.partition_lat(cur,part_cur,mid,nb_cells,nb_regions,\
            square_side)
        i = i + 1

      part_cur = part_cur + nb_regions 
      nb_regions = nb_regions + 2
    self.nb_regions_beforelast = nb_regions - 2
    return part_cur,cur

  def place_remaining_points(self,part_cur,cur,square_side):
    
    i_cur = self.nb_part2*square_side + 1
    nb_lat_left = self.grid.nb_lat - i_cur + 1

    nb_points_left = self.grid.get_size_zone() - square_side**2*self.nb_part2**2
    nb_square_left = int(float(nb_points_left/(square_side**2)))
    rect_height = nb_lat_left
    rect_width = square_side**2 / nb_lat_left
    nb_regions = nb_points_left / rect_width
    print rect_width,rect_height
    mid = 0

    self.strategy = 3
    part_save = part_cur
    for i in range(i_cur,self.grid.nb_lat+1):
      nb_cells = self.grid.get_nb_cells_zone(i)
      cur = self.partition_lat(cur,part_cur,mid,nb_cells,nb_regions,\
          rect_width)
      part_cur = part_save

  def partitioning(self):
    # We partition into equal squares first
    # Number of squares
    self.nb_part2 = int(np.sqrt(self.nb_part))

    # Length of a side
    square_side = int((self.grid.nb_lat-1)/self.nb_part2)

    part_cur,cur = self.split_to_squares(square_side)
    # Then place the remaining points in a last band
    self.place_remaining_points(part_cur,cur,square_side)
    
  #-------------------------------------------------------------------------------
  # Assing each cell in a line to a partition
  #-------------------------------------------------------------------------------
  def partition_lat(self,cur,part_cur,mid,nb_cells,nb_regions,split_size):
    # Rest and quotient
    r = nb_cells % nb_regions
    q = nb_cells / nb_regions
    diff = r
  
    left = 0
    size = 0
    band_cur = 0
    # band_cur = -1 serves as initialisation
    size,diff = self.next_size(size,q,diff,-1,mid,nb_regions,nb_cells,split_size,\
        part_cur)
    right = left + size - 1
  
    for k in range(nb_cells):
      self.part[cur] = part_cur + band_cur
      if (k == right):
        left = left + size
        size,diff = self.next_size(size,q,diff,band_cur,mid,nb_regions,\
            nb_cells,split_size,part_cur)
        right = left + size - 1
        band_cur = band_cur + 1
      cur = cur+1
    return cur

  #-------------------------------------------------------------------------------
  # Set the size of the next partition at a given latitude
  #-------------------------------------------------------------------------------
  def next_size(self,size,q,diff,band_cur,mid,nb_regions,nb_cells,split_size,\
      part_cur):
    if (self.strategy == 1):
      # The error in the division is set gradually
      size = q
      if (diff > 0):
        size = size + 1 
        diff = diff - 1
    elif (self.strategy == 2):
      # No division but each partition is a square (except at the center, where
      # it is triangular)
      size = split_size
      if (band_cur == mid-1):
        size = nb_cells - (nb_regions-1)*split_size
    elif (self.strategy == 3):
      actual_part = part_cur + band_cur
      # If not juste before the last
      if (actual_part < self.nb_part - 2):
        size = split_size
      else:
        # Fill until the end for last partition
        size = nb_cells - (nb_regions-1)*split_size
    else:
      # The error in the division is set at the center
      size = q
      if (band_cur == mid-1 and diff > 0):
        size = size + diff 
    return size,diff
  
  #-------------------------------------------------------------------------------
  # Output
  #-------------------------------------------------------------------------------
  # Call metapost and ghostview
  def display(self):
    cur_dir = os.getcwd()
    os.chdir(self.mp_dir)
    mp_file = 'partition.mp'
    os.system('mpost '+mp_file)
    mp_file_out = 'partition-1.mps'
    os.system('gv '+mp_file_out)



  # Write (cell number,partition number, colors) 
  # Display is done by metapost
  def metapost_output(self):
    self.list_partitions = self.mp_dir+'list_partitions.txt'
    f = open(self.list_partitions,'w')
    k = 0
    for cur in self.part:
      x = int(cur) % 16
      s = "%d,%d,%f,%f,%f \n" % (k,x,self.colors[x][0],self.colors[x][1],\
          self.colors[x][2])
      f.write(s)
      k = k + 1
    f.close()


  #-------------------------------------------------------------------------------
  # Init color map for display
  #-------------------------------------------------------------------------------
  def init_colors(self):
        # Red
    self.colors = [[1.00, 0.00, 0.00], \
      # Green        
      [0.00, 1.00, 0.00],         \
      # Yellow       
      [1.00, 1.00, 0.00],         \
      # Blue         
      [0.00, 0.00, 1.00],         \
      # Magenta      
      [1.00, 0.00, 1.00],         \
      # Cyan         
      [0.00, 1.00, 1.00],         \
      # Orange       
      [1.00, 0.50, 0.20],         \
      # Olive        
      [0.30, 0.55, 0.00],         \
      # Dark pink    
      [0.72, 0.47, 0.47],         \
      # Sea blue     
      [0.33, 0.33, 0.81],         \
      # Pink         
      [1.00, 0.63, 0.63],         \
      # Violet       
      [0.62, 0.44, 0.65],         \
      # Pale green   
      [0.60, 0.80, 0.70],         \
      # Brown        
      [0.47, 0.20, 0.00],         \
      #Turquoise    
      [0.00, 0.68, 0.68],         \
      #! Purple     
      [0.81, 0.00, 0.40]] 
    self.nb_colors = len(self.colors)
 
  #-------------------------------------------------------------------------------
  # Communication costs
  #-------------------------------------------------------------------------------
  # Communication cost for all partitions in one pass
  def comm_cost(self):
    k = 0
    cost = np.zeros(self.nb_part)
    nb_neighbours = np.zeros(self.nb_part)
    volume = np.zeros(self.nb_part)
    for i in range(1,self.grid.nb_lat+1):
      nb_cells = self.grid.get_nb_cells_zone(i)
      for j in range(1,nb_cells+1):
        cur = self.part[k]
        tmp = self.grid.find_neighbours_zone(i,j)
        neighbours = self.grid.position_to_indice(tmp)
        # If the neighbours is not in the current partition, then we add it to
        # the cost
        for x in neighbours:
          if(self.part[x] == cur):
            volume[cur] = volume[cur] + 1
          else:
            nb_neighbours[cur]= nb_neighbours[cur]+1
        k = k + 1

    # Beware, the volume must be divided by 2, as the communications are one way
    # only (so counted twice)
    volume = volume * 0.5
    for k in range(self.nb_part):
      cost[k] = float(nb_neighbours[k])/volume[k]
    print "volume=",volume
    print "Nb neighbours=",nb_neighbours
    #print cost
    return cost

  # Simple plot for cost
  def plot_cost(self,cost):
    fig = plt.figure()
    plt.plot(cost)
    plt.title("Communication cost/volume")
    plt.xlabel("Partition number")
    plt.ylabel("Ratio")
    #plt.show()
    plt.savefig("partition_cost.png")

