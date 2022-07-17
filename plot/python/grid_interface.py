#-------------------------------------------------------------------------------
# Number of cells at line i on all sectors
# n_lat2 : number of latitude lines on an hemisphere
#-------------------------------------------------------------------------------
def nb_cells(i, n_lat2):
  n = 2*i-1
  if (i > n_lat2) :
    n = 2*(2*n_lat2 - i + 1)-1
  return n
 

#-------------------------------------------------------------------------------
# Print each cell corners in a file
#-------------------------------------------------------------------------------
def print_grid():

  for i in range(1, n_lat+1):
    n_cells = nb_cells(i, n_lat2)
    dlon = 360./n_cells
    lat = 90.-(i-1)*dlat
    for j in range(1, n_cells+1):
      print lat, (j-1)*dlon
      print lat, j*dlon
      print lat - dlat, j*dlon
      print lat - dlat, (j-1)*dlon


#-------------------------------------------------------------------------------
# Find the cell grid corresponding to the current point (must be in degree)
# lat0 must be between 90 and -90
# lon0 must be between 0 and 360
# Cell number starts from 1
#-------------------------------------------------------------------------------
def find_cell(lat0, lon0):
  i = int((90.-lat0)/dlat) + 1
  n_cells = nb_cells(i, n_lat2)
  dlon = 360./n_cells
  j = int(lon0/dlon) + 1
  print "Found cell", i, j
  # Check
  if (lat0 > 90.-(i-1)*dlat or lat0 < 90.-i*dlat):
    print "do not fit for lat. May be on the border"
  if (lon0 < (j-1)*dlon or lon0 > j*dlon):
    print "do not fit for lon. May be on the border"

#-------------------------------------------------------------------------------
# Define n_eq number of cells at the equator. This will be approximated
#-------------------------------------------------------------------------------
# Modify this
n_eq = 3

#-------------------------------------------------------------------------------
n_lat2 = ((n_eq / 3) + 1)/2
n_lat = 2*n_lat2
n_eq2 = nb_cells(n_lat2, n_lat2)
dlat = 90./n_lat2
print "Number of latitudes", n_lat
print "True number of cells at the equator", n_eq2
print "Dlat", dlat


#print_grid()
#find_cell(-89.5, 89.5)
