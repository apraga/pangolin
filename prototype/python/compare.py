####################################################
# Compare analytical and numerical results
####################################################
import numpy as np
import sys 
import subprocess as sub
import lib_plot # personal library


nb_args = 2;
if (len(sys.argv) == nb_args+1):
  plot_type = sys.argv[1]
  nb_iter = int(sys.argv[2])
else:
  print "Usage : "+str(sys.argv[0])+" PLOT NB_ITER"
  print "plot = mill, ortho, npolar, spolar"
  exit(2)

#version = "_C=0.5"
#name1 = "../data/data"+version+"/"+str(nb_iter)+".dat"
#name2 = "../data/data_python"+version+"/"+str(nb_iter)+".dat"
name1 = "../data/data/"+str(nb_iter)+".dat"
name2 = "../data/data_python/"+str(nb_iter)+".dat"

nb_lines1 = lib_plot.file_len(name1)
nb_lines2 = lib_plot.file_len(name2)
if (nb_lines1 != nb_lines2):
  print "Number of lines different !"
nb_lines = min(nb_lines1,nb_lines2)

lat = np.zeros(nb_lines,float)
lon = np.zeros(nb_lines,float)
diff =  np.zeros(nb_lines,float)

k = 0
max_q = -10000
f1 = open(name1,'r')
for line in f1.readlines():
  cur = line.split()
  lat[k] = float(cur[0])
  lon[k] = float(cur[1])
  if (lon[k] > 180):
    lon[k] = lon[k] - 360.
  diff[k] = float(cur[2])
  max_q = max(max_q,float(cur[2]))
  k = k + 1
f1.close()

k = 0
f2 = open(name2,'r')
for line in f2.readlines():
  cur = line.split()
  diff[k] = float(cur[2]) - diff[k]
  max_q = max(max_q,float(cur[2]))
  k = k + 1
f2.close()

diff = diff/max_q

###############################################################################
# Plotting
###############################################################################
name_fig = 'bidim_error_'+str(nb_iter)+'.jpg'
title = "Relative error" 

lib_plot.plot(lat,lon,diff,plot_type,title,name_fig)

## On a map
##map = Basemap(projection='mill',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=0,urcrnrlon=190,resolution='c')
##map = Basemap(llcrnrlon=-180,llcrnrlat=-80,urcrnrlon=180,urcrnrlat=80,projection='mill')
#if (plot_type == 'ortho'):
#  m = Basemap(projection='ortho',lat_0 = 50, lon_0 = -100,\
#      resolution = 'l', area_thresh = 1000.)
#elif (plot_type == 'npolar'):
#  m = Basemap(projection='npstere',boundinglat=70,lon_0=270,resolution='l')
#elif (plot_type == 'spolar'):
#  m = Basemap(projection='spstere',boundinglat=-70,lon_0=270,resolution='l')
#else:
#  m = Basemap(projection='mill',llcrnrlat=-90,urcrnrlat=90,\
#             llcrnrlon=-180,urcrnrlon=180)
#
#m.drawcoastlines()
#m.drawmapboundary()
#
#x,y = m(lon,lat)
#m.scatter(x,y, s=10, c=diff, marker='o',edgecolors='none')
#
#m.colorbar()
#p.title('Erreur analytique - numerique (projection de Miller)')
#name_fig = 'bidim_error_'+str(nb_iter)+'.jpg'
##name_fig = 'bidim_error_'+str(nb_iter)+'.eps'
#p.savefig(name_fig)
#sub.call(["eog", name_fig])
