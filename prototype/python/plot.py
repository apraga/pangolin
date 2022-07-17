import sys 
import lib_plot # personal library

###############################################################################
# Wrapper for plotting function
###############################################################################
nb_args = 3;
if (len(sys.argv) == nb_args+1):
  data = sys.argv[1]
  plot_type = sys.argv[2]
  nb_iter = int(sys.argv[3])
else:
  print "Usage : "+str(sys.argv[0])+" DATA PLOT NB_ITER"
  print "data = num, anal"
  print "plot = mill, ortho,zoom_ortho, npolar, spolar"
  exit(2)

if (data == 'num'):
  folder = "../matlab/data/"
else:
  folder = "../matlab/data/data_python/"

name = folder+str(nb_iter)+".dat"
title = "Cosine bell distribution" 
name_fig = 'bidim_'+data+'_'+str(nb_iter)+'.png'

lib_plot.plot_data(name,plot_type,title,name_fig)
