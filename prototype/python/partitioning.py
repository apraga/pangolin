#-------------------------------------------------------------------------------
# Compute a partitioning and write it to a file for postprocessing with
# metapost
#-------------------------------------------------------------------------------

import grid_class
import partition_class


#-------------------------------------------------------------------------------
# Partition the reduced grid
#-------------------------------------------------------------------------------
g = grid_class.Grid(1.,90)
p = partition_class.Partition(g,133)
p.strategy = 2
p.partitioning()
cost = p.comm_cost()
p.plot_cost(cost)

p.metapost_output()
p.display()
