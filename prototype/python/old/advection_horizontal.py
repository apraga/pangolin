#import sys
#sys.path.append("~/code/prototype/python")
import init # Local 

nb_lat = 90
nb_lon = 90
delta_t = 20

# Init concentration
C = init.concentration(nb_lat,nb_lon)
(step,pos) = init.grid(nb_lat)
#(node_middle,node_step,node_prev,node_next) = init.nodes(nb_lat)
#print node_middle


