#----------------------------------------------------------------------
# Contains everything used for plotting numerical and analytical data
# Shoud be compatible for python 2 and 3
#----------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as p
#from scipy.interpolate import griddata
from matplotlib.mlab import griddata
import subprocess as sub
from mpl_toolkits.basemap import Basemap, cm

prec = 1.e-13
#------------------------------------------------------------------------------
# Plotting the data in a file 
# Find the format automatically (either concentration or winds)
#------------------------------------------------------------------------------
def plot_data(filename, projection, interpol, iserr, is_winds, title, name_fig):

  m  = define_projection(projection)

  if (is_winds):
    print("data format is : winds")
    u, v, lat_u, lon_u, lat_v, lon_v = read_winds(filename)
    draw_vectors(u, v, lat_u, lon_u, lat_v, lon_v, m)
  else:
    print("data format is : concentration")
    q, lat, lon = read_concentration(filename)

    if interpol is None:
      print("Scatter data")
      draw_scatter(q, lat, lon, m, projection, iserr)
    else:
      print("Interpolation")
      draw_contour(q, lat, lon, m, projection, iserr)

  finalize(m, title, name_fig, is_winds)

#------------------------------------------------------------------------------
# Gives zonal file and retrieve meridional file
# Format : u/v, cell center (lat, lon)
#------------------------------------------------------------------------------
def read_winds(u_fname):
  data_u = [line.strip().split() for line in open(u_fname)]
  v_fname = u_fname.replace("u_", "v_")
  data_v = [line.strip().split() for line in open(v_fname)]

  # Select data
  u = [float(a) for (a,b,c) in data_u]
  lat_u = [float(b) for (a,b,c) in data_u]
  lon_u = [float(c) for (a,b,c) in data_u]

  v = [float(a) for (a,b,c) in data_v]
  lat_v = [float(b) for (a,b,c) in data_v]
  lon_v = [float(c) for (a,b,c) in data_v]

  if (len(lat_v) != len(lat_u) or len(lon_u) != len(lon_v)):
    raise Exception("We can only plot winds with both composantes on the same\
    point")
  #Convert to -180, 180 for plot
  lon_u = [ cur - 360 if (cur > 180) else cur for cur in lon_u]
  lon_v = [ cur - 360 if (cur > 180) else cur for cur in lon_v]

  return u, v,lat_u, lon_u, lat_v, lon_v


#------------------------------------------------------------------------------
# Format : cell value, cell center (2)
#------------------------------------------------------------------------------
def read_concentration(filename):
  data = [line.strip().split() for line in open(filename)]

  # Select data
  q = [float(a) for (a,b,c) in data]
  lat = [float(b) for (a,b,c) in data]
  lon = [float(c) for (a,b,c) in data]

  #Convert to -180, 180 for plot
  lon = [ cur - 360 if (cur > 180) else cur for cur in lon]

  return q, lat, lon


def define_projection(plot_type):
  res = 'l';
  if (plot_type is None):
    print("No projection")
    return None

  # Orthographic projection
  if (plot_type == 'ortho'):
    m = Basemap(projection='ortho',lat_0 =90, lon_0 =0,\
        resolution = 'l', area_thresh = 1000.)
  # Orthographic projection with zoom
  elif (plot_type == 'zoom_ortho'):
    lat_0=90.0
    lon_0=0.
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

  return m


def draw_scatter(q, lat, lon, m, plot_type, iserr):
  if (plot_type == 'ortho'):
    size=10
  else:
    size=5

  # For error, no fixed range
  #opts = {'vmin':0.,'vmax':1.} if (iserr is None) else {}
  # Custom limits : colormap and axe
  opts = {} if (iserr is None) else {'vmin':-0.5,'vmax':0.1}
  #p.ylim((-90,90))

  if m is not None:
    x,y = m(lon,lat)
    m.scatter(x,y, s=size, c=q, marker='s',edgecolors='none', **opts)
  else:
    p.scatter(lon, lat, s=size, c=q, marker='s',edgecolors='none',**opts)

def draw_contour(q, lat, lon, m, plot_type, iserr):
  npts = 500
  loni = np.linspace(0.,360.,npts)
  lati = np.linspace(-90., 90.,npts)
  qi = griddata((lon, lat), q, (loni[None,:], lati[:,None]), method='cubic')

  # For error plat, fixed range
  #opts = {} if (iserr is None) else {'levels':np.linspace(-1.,1, 200)}
  # Setting the levels is enough for colormap min and max
  opts = {'levels':np.linspace(0.,1, 50)}


  if plot_type is not None:
    loni, lati = np.meshgrid(loni, lati)
    x, y = m(loni, lati)
  
    # Filled, use contour for non-filled
    cs = m.contourf(x, y, qi, 15,cmap=p.cm.jet)
  else:
    cs = p.contourf(loni, lati, qi,cmap=p.cm.jet, **opts)

# For projection, needs interpolation
def draw_vectors(u, v, lat_u, lon_u, lat_v, lon_v, m):
  if (m is None):
    p.quiver(lon_u, lat_u, u, v, color='r')#, scale=1./0.3)#width=0.002)#headlength=7)
    return

  lat_u = np.asarray(lat_u)
  lon_u = np.asarray(lon_u)
  print("Interpolation for vectors")
  npts = 500
  loni = np.linspace(-180.,180.,npts)
  lati = np.linspace(-90., 90.,npts)
  #v = [0. for x in v]
  #u = [0. for x in u]
  u_i = griddata(lon_u, lat_u, u, loni, lati)
  v_i = griddata(lon_v, lat_v, v, loni, lati)

  uproj,vproj,xx,yy = \
      m.transform_vector(u_i, v_i, loni, lati,101,101,returnxy=True)#,masked=True)

  m.quiver(xx, yy, uproj, vproj)#

def finalize(m, title, name_fig, is_winds):
  if m is None and not is_winds:
    p.colorbar()
  elif not is_winds: 
    m.colorbar()
  p.title(title)
  # Saving it and display it with eog if issue

  p.savefig(name_fig)
  p.show()
  #sub.call(["eog", name_fig])

#------------------------------------------------------------------------------
# Get number of lines in file
#------------------------------------------------------------------------------
def file_len(fname):
  p = sub.Popen(['wc', '-l', fname], stdout=sub.PIPE, 
    stderr=sub.PIPE)
  result, err = p.communicate()
  if p.returncode != 0:
    raise IOError(err)
  return int(result.strip().split()[0])
