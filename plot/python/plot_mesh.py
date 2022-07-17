from mpl_toolkits.basemap import Basemap
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def draw_polygon(i, j, dlat, dlon, m, lat_0, lon_0, color=None):
  poly = add_to_polygon(90.-i*dlat, (j-1)*dlon, 90.-(i-1)*dlat,
      (j-1)*dlon, m, lat_0, lon_0)
  poly = poly + add_to_polygon(90.-(i-1)*dlat, (j-1)*dlon, 90.-(i-1)*dlat,
      j*dlon, m, lat_0, lon_0)
  poly = poly + add_to_polygon(90.-(i-1)*dlat, j*dlon, 90.-i*dlat, j*dlon, m,
      lat_0, lon_0)
  poly = poly + add_to_polygon(90.-i*dlat, j*dlon, 90.-i*dlat, (j-1)*dlon, m,
      lat_0, lon_0)

  if (len(poly) > 0):
    if (color):
      poly = Polygon(poly, facecolor=color, alpha=0.4)
    else:
      poly = Polygon(poly, facecolor='white')
    plt.gca().add_patch(poly)


def add_to_polygon(lat0, lat1, lon0, lon1, m, lat_0, lon_0):
  resolution = 10
  lats = np.linspace(lat0, lon0, resolution )
  lons = np.linspace(lat1, lon1, resolution )
  lats_f = []
  lons_f = []
  coef = math.pi/180
  for i in range(len(lats)):
    cos_c = np.sin(lat_0*coef)*np.sin(lats[i]*coef) + np.cos(lat_0*coef)*np.cos(lats[i]*coef)*np.cos((lons[i] - lon_0)*coef)
    cond = True
    if (proj == 'ortho'): cond = cos_c >= 0
    #if (cos_c >= 0):
    if cond:
      lats_f.append(lats[i])
      lons_f.append(lons[i]-offset)
  x,y = m(lons_f, lats_f)
  res = zip(x, y)
  return res


nb_lat2 = 10
nb_lat = 2*nb_lat2
dlat = 90./nb_lat2
nb_cells = nb_lat**2
poly = []
k = 0
offset = 0

lat_0, lon_0 = 30, 0
proj = 'robin'
has_color = True

# Need to be in [-180,180]
if (proj != 'ortho'): offset = 180
m = Basemap(projection=proj,lon_0=lon_0, lat_0=lat_0)
#m.drawcoastlines()
m.drawmapboundary()


for i in range(1, nb_lat+1):
  nb_lon2 = 2*i-1
  if (i > nb_lat2) :
    nb_lon2 = 2*(nb_lat - i + 1)-1
  dlon = 120./nb_lon2
  nb_lon = 3*nb_lon2

  for j in range(1, nb_lon+1):
    if (has_color):
      if (i < nb_lat2+1):
        if (j < nb_lon2+1): 
          color = 'red'
        elif (j < 2*nb_lon2+1):
          color = 'blue'
        else:
          color = 'green'
      else:
        if (j < nb_lon2+1): 
          color = 'gray'
        elif (j < 2*nb_lon2+1):
          color = 'yellow'
        else:
          color = 'white'

      draw_polygon(i, j, dlat, dlon, m, lat_0, lon_0, color)
    else:
      draw_polygon(i, j, dlat, dlon, m, lat_0, lon_0)

#plt.show()
name = "grid"
if (has_color):
  name = name+"_color"

name = name+".pdf"
plt.savefig(name, format='pdf')
