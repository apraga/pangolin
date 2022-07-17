#######################################################################
# Contains everything used for plotting numerical and analytical data
#######################################################################

import numpy as np
import matplotlib.pyplot as p
import subprocess as sub
from mpl_toolkits.basemap import Basemap

###############################################################################
# Plotting the data in a file 
###############################################################################
def plot_data(filename,plot_type,title,name_fig):
  nb_lines = file_len(filename)
  
  lat = np.zeros(nb_lines,float)
  lon = np.zeros(nb_lines,float)
  q = np.zeros(nb_lines,float)
  k = 0
  f = open(filename,'r')
  for line in f.readlines():
    cur = line.split()
    lat[k] = float(cur[1])
    lon[k] = float(cur[2])
    q[k] = float(cur[0])
    if (lon[k] > 180):
      lon[k] = lon[k] - 360.
    k = k + 1
  f.close()
  plot(lat,lon,q,plot_type,title,name_fig)

###############################################################################
# Plotting function : need latitude, longitude, data
###############################################################################
def plot(lat,lon,q,plot_type,title,name_fig):
  res = 'l';
  # Orthographic projection
  if (plot_type == 'ortho'):
    lat_0=0.; 
    lon_0=0.
    m = Basemap(projection='ortho',lat_0 = lat_0, lon_0 = lon_0,\
        resolution = 'l', area_thresh = 1000.)
  # Orthographic projection with zoom
  elif (plot_type == 'zoom_ortho'):
    lat_0=0.; lon_0=0.
    m1 = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution=None)
    
    width = m1.urcrnrx - m1.llcrnrx 
    height = m1.urcrnry - m1.llcrnry
    
    coef = 0.8
    width = width*coef
    height = height*coef
    m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution=res,\
        llcrnrx=-0.5*width,llcrnry=-0.5*height,urcrnrx=0.5*width,urcrnry=0.5*height)

  # Polar projection (north)
  elif (plot_type == 'npolar'):
    m = Basemap(projection='npstere',boundinglat=70,lon_0=270,resolution=res)
  # Polar projection (south)
  elif (plot_type == 'spolar'):
    m = Basemap(projection='spstere',boundinglat=-70,lon_0=270,resolution=res)
  # Miller projection 
  else:
    m = Basemap(projection='mill',llcrnrlat=-90,urcrnrlat=90,\
               llcrnrlon=-180,urcrnrlon=180,resolution=res)
  
  m.drawcoastlines()
  m.drawmapboundary()
  #m.drawparallels(np.arange(-90, 90, 30),linewidth=2.0,dashes=[50,1])
  #m.drawmeridians(np.arange(-180, 180, 20),linewidth=2.0,dashes=[50,1])
  
  x,y = m(lon,lat)
  m.scatter(x,y, s=7, c=q, marker='s',edgecolors='none')
  
  m.colorbar()
  p.title(title)
  p.savefig(name_fig)
  # As we cannot see it...
  sub.call(["eog", name_fig])

###############################################################################
# Get number of line in file
###############################################################################
def file_len(fname):
  p = sub.Popen(['wc', '-l', fname], stdout=sub.PIPE, 
    stderr=sub.PIPE)
  result, err = p.communicate()
  if p.returncode != 0:
    raise IOError(err)
  return int(result.strip().split()[0])
