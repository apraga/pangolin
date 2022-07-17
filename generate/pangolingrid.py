#-------------------------------------------------------------------------------
# Class for creating a Pangolin grid and correcting winds
# Also contains a class for generating ratio and winds for test cases
# Should be compatible for python 2 and 3
#-------------------------------------------------------------------------------

from ctypes import *
import numpy as np
import random
try:
  import h5py
  has_hdf5 = True
except ImportError:
  has_hdf5 = False

# Winds are defined in m/s

if (has_hdf5):
    hdf5_type = h5py.h5t.NATIVE_DOUBLE
fformat = "%23.16f"
wformat = fformat+" "+fformat+" "+fformat+"\n"
coef = np.pi/180.
precision = 1e-13;
# Earth radius in meters
R_earth = 6.3172e6;
R = R_earth

#-------------------------------------------------------------------------------
# Pangolin grid
# Note : use new-class for calling parent constructor later
#-------------------------------------------------------------------------------
class PangolinGrid(object):
  """ An area-preserving grid.
  
  The grid is complety defined by the number of latitudes. We also store half
  the number of latitudes.
  
  This class contains function for generating a grid from analytical cases,
  interpolating data to and from this grid. We can also correct the winds for a
  null divergence. Most of the functions needs an input file to work on.case
  
  We assume the cell height is constant, so the cell area is not completely
  preserved
  """

  def __init__(self, wind_format):
    """ Default constructor. Needs output cell format : center or borders for
    winds
    """
    self.nb_lat = 0
    if (wind_format != "center" and wind_format != "borders"):
      raise Exception("cell format must be center or borders")
    self.wind_format = wind_format
  
  def init(self, nb_lat, io_type="hdf5"):
    """ nb_lat must be even. Generate a grid with 6*nb_lat^2 cells
    """
    if (not isinstance(nb_lat, int) or nb_lat < 1):
      raise Exception("nb_lat must be a strictly positive integer")
    if (nb_lat % 2 != 0):
      raise Exception("Needs an even number of latitudes")
    if (io_type != "hdf5" and io_type != "ascii"):
      raise Exception("IO type is either hdf5 or ascii")
    if (io_type == "hdf5" and not has_hdf5):
      raise Exception("Module h5py not found and required")

    self.nb_lat = nb_lat
    self.nb_lat2 = nb_lat/2
    print("nb lat2 %d"% self.nb_lat2)
    self.dlat = 90./self.nb_lat2
    # Number of cells
    self.size_q = 6*self.nb_lat2**2
    # Number of zonal winds
    self.size_u = 6*self.nb_lat2**2 + self.nb_lat
    # Number of merid winds. This is found by summing the number of nodes
    # (3.(4n-1)) up to n-1, multiply by 2 to have both hemisphere and add the
    # special case at the Equator (3*(2n-1))
    self.size_v = 12*self.nb_lat2**2 - 12*self.nb_lat2 + 3
    # HDF5 is the default, but ascii is also possible
    self.io_type = io_type
    self.ext = ".h5"
    if (io_type != "hdf5"):
      self.ext = ".dat"

  def is_valid_line(self, line):
    return (not "#" in line and line.strip())

  def valid_lines(self, f):
    """ Ignore comments and empty lines
    """
    for line in f:
      if (self.is_valid_line(line)):
        yield line

  def next_valid_line(self, fi):
    line = fi.readline()
    while (not self.is_valid_line(line)):
      line = fi.readline()
    return line

  def nb_lines(self, fname):
    """ Number of valid lines in a file (omits comments and blank lines)
    """
    nb = 0
    with open(fname) as f:
      for line in self.valid_lines(f):
        nb += 1
    return nb
     

  def init_from_file(self, fname):
    """ Create an instance from a file containing the cell ratio.
    We just count the number of cells (ignoring the comments and empty lines)
    Also check the file format.
    """
    k = 0
    print("Init from file")
    with open(fname) as f:
      for line in self.valid_lines(f):
        col = line.split()
        if (len(col) != 3):
          raise Exception("Wrong line in file : should be (q, lat, lon) \n"+\
              "Line : "+line)
        k += 1
   
    nb_lat2 = int(np.sqrt(k/6))
    n = 6*nb_lat2**2
    if ( n != k):
      raise Exception("Wrong number of cells in : needs 6n^2 cells\n"+\
          "Found "+str(n))
    self.init(2*nb_lat2)

  def iterator(self):
    """ Iterate over all the cells and return a couple (i, j)
    """
    for i in range(1, self.nb_lat+1):
      for j in range(1, self.nb_cells(i)+1):
        yield [i, j]
 

  def nb_cells(self, i):
    """ Total number of cells at a latitude line (numbered from the north pole)
    """
    return 3*self.nb_cells_sector(i)

  def nb_cells_sector(self, i):
    """ Number of cells in one sector (there are 3)
    """
    if (i < 1 or i > self.nb_lat):
      raise Exception("i outside range")

    if (i > self.nb_lat2):
      i = self.nb_lat - i + 1
    return 2*i-1

  def nb_nodes(self, i):
    """ Number of nodes at latitude 90 -i , i being > 0
    """
    return 3*self.nb_nodes_sector(i)

  def nb_nodes_sector(self, i):
    if (i < 1 or i > self.nb_lat-1):
      raise Exception("i outside range")

    if (i == self.nb_lat2):
      return self.nb_cells_sector(i)

    if (i > self.nb_lat2):
      i = self.nb_lat - i 
    return 4*i-1

  def north_neighbours(self, i, j):
    """ North neighbours of cells at (i, j). Returns an array of 2 integers for
    the max and min
    """
    if (i < self.nb_lat2 + 1):
      return self.north_neighbours_northern(i, j)
    elif (i == self.nb_lat2 + 1):
      return [j]
    else:
      i = self.nb_lat - i + 1
      return self.south_neighbours_northern(i, j)

  def south_neighbours(self, i, j) :
    """ South neighbours of cells at (i, j). Returns all neighbours
    """
    if (i < self.nb_lat2):
      return self.south_neighbours_northern(i, j)
    elif (i == self.nb_lat2):
      return [j];
    else:
      i = self.nb_lat - i + 1
      return self.north_neighbours_northern(i, j)

  def switch_first_sector(self, i, j):
    nb_cells_sector = self.nb_cells_sector(i)
    j2 = j % nb_cells_sector
    if (j2 == 0): j2 = nb_cells_sector
    sector = (j-1)/nb_cells_sector

    mid = (nb_cells_sector + 1)/2
    return [j2, sector, mid]

  def north_neighbours_northern(self, i, j):
    """ Only works for the northern hemisphere
    """
    nb_cells = self.nb_cells(i)
    nb_cells_next = self.nb_cells(i-1)
    if (j < 1 or j > max(nb_cells, nb_cells_next)):
      raise Exception("j outside range")

    [j2, sector, mid] = self.switch_first_sector(i, j)

    nb_cells = self.nb_cells_sector(i)
    nb_cells_next = self.nb_cells_sector(i-1)
    if (j2 == 1):
      j_neighb = [1]
    elif (j2 == nb_cells):
      j_neighb = [nb_cells_next]
    elif (j2 == mid):
      j_neighb = [j2-1]
    else:
      offset = 0
      if (j2 > mid):
        offset = -1
      j_neighb = [j2-1 + offset, j2 + offset]

    # Going back to complete grid
    j_neighb = [x + sector*nb_cells_next for x in j_neighb]
    return j_neighb


  def south_neighbours_northern(self, i, j):
    """ Only works for the northern hemisphere
    """
    nb_cells = self.nb_cells(i)
    nb_cells_next = self.nb_cells(i+1)
    if (j < 1 or j > max(nb_cells, nb_cells_next)):
      raise Exception("j outside range")

    [j2, sector, mid] = self.switch_first_sector(i, j)

    if (j2 == nb_cells):
      j_neighb = [j2]
    elif (j2 == mid):
      j_neighb = [j2, j2+1, j2 + 2]
    else:
      offset = 0
      if (j2 > mid):
        offset = 1
      j_neighb = [j2 + offset, j2 + 1 + offset]

    # Going back to complete grid
    nb_cells_sector = self.nb_cells_sector(i+1)
    j_neighb = [x + sector*nb_cells_sector for x in j_neighb]
    #if (i == 2): print "i j jneighb", i, j, j_neighb
    return j_neighb

  def generator_center(self):
    """ Returns a list of grid cell centers
    """
    for i in range(1, self.nb_lat+1):
      nb_lon = self.nb_cells(i)
      lat = 90 - (i-0.5)*self.dlat
      dlon = 360./nb_lon
      # Constant for all j
      for j in range(nb_lon):
        lon = (j+0.5)*dlon
        yield lat, lon

  def generator_zonal(self):
    """ Returns a list of grid cell center of zonal borders
    """
    for i in range(1, self.nb_lat+1):
      lat = 90 - (i-0.5)*self.dlat
      nb_lon = self.nb_cells(i)
      dlon = 360./nb_lon

      # Write only the left border and the right at the end of the grid
      # which should be the same as the first by periodicity
      for j in range(nb_lon+1):
        lon = j*dlon
        yield lat, lon

  def generator_merid(self):
    """ Returns a list of grid cell center of meridional borders
    """
    # We could use the cell neighbours to find the interface length
    for i in range(1, self.nb_lat):
      lat = 90. - i*self.dlat
      nb_lon = self.nb_cells(i)
      dlon = 360./nb_lon
      dlon_next = self.get_dlon(i+1)
      for j in range(1, nb_lon+1):
        j_neighb = self.south_neighbours(i, j)
        for k in j_neighb:
          lon_max = min(j*dlon, k*dlon_next)
          lon_min = max((j-1)*dlon, (k-1)*dlon_next)
          lon = lon_min + 0.5*(lon_max - lon_min)
          yield lat, lon

  def generator(self, dtype):
    """ Generate cells coordinates for ratio or winds
    """
    if (dtype == "ratio"):
      coords = self.generator_center()
    else:
      if (self.wind_format == "center"):
        coords = self.generator_center()
      else:
        if (dtype == "zonal"):
          coords = self.generator_zonal()
        else:
          coords = self.generator_merid()
    return coords

  def generate_ratio(self, qfile, case):
    """ Write cell concentration at the cell centers and cell center coordinates
    """
    qfile += self.ext
    if (self.io_type == "hdf5"):
      f = h5py.File(qfile,'w')
      self.generate_data_hdf5(f, case, self.size_q, "ratio")
    else:
      self.generate_ratio_ascii(qfile, case)

  def generate_data_hdf5(self, f, case, size, data_type, t1=None, t2=None):
    """ Generic data write : works for ratio and winds
        Write an array for each coordinates and data
        File must be opened 
    """
    #Faster to use arrays, but it suppose we have enough memory
    center_lat = np.zeros(size)
    center_lon = np.zeros(size)
    d = np.zeros(size)
    coords = self.generator(data_type)
    k = 0
    for (lat, lon) in coords:
      center_lat[k] = lat
      center_lon[k] = lon
      if (data_type == "ratio"):
        d[k] = case.get_ratio(lat, lon)
      elif (data_type == "zonal"):
        d[k] = case.get_u(lat, lon, t1, t2)
      elif (data_type == "merid"):
        d[k] = case.get_v(lat, lon, t1, t2)
      k += 1

    # Write array to hdf5 file
    name = "ratio"
    if (data_type == "zonal"): name = "u"
    if (data_type == "merid"): name = "v"
    dset_d = f.create_dataset(name, (size,), dtype=hdf5_type)
    dset_d[:] = d[:]
    dset_lat = f.create_dataset("lat", (size,), dtype=hdf5_type) 
    dset_lat[:] = center_lat[:]
    dset_lon = f.create_dataset("lon", (size,), dtype=hdf5_type)
    dset_lon[:] = center_lon[:]
    f.close()

  def generate_ratio_ascii(self, qfile, case):
    f = open(qfile, 'w')
    print("Generating %s"% qfile)
    coords = self.generator("ratio")
    for (lat, lon) in coords:
      val = case.get_ratio(lat, lon)
      self.write_data(val, lat, lon, f)
    f.close()

  def write_data(self, val, lat, lon, f):
      f.write(wformat % (val, lat, lon))

  def generate_u(self, qfile, case, t1=None, t2=None):
    """ Write zonal winds at the east and west cell borders
    """
    qfile += self.ext
    if (self.io_type == "hdf5"):
      f = h5py.File(qfile,'w')
      self.generate_data_hdf5(f, case, self.size_u, "zonal", t1, t2)
    else:
      self.generate_u_ascii(qfile, case, t1, t2)
 
  def generate_u_ascii(self, qfile, case, t1=None, t2=None):
    print("Generating %s"% qfile)
    f = open(qfile, 'w')
    coords = self.generator("zonal")
    for (lat, lon) in coords:
      val = case.get_u(lat, lon, t1, t2)
      self.write_data(val, lat, lon, f)

    f.close()

  def generate_v(self, qfile, case, t1=None, t2=None):
    """ Write meridional winds at the east and west cell borders
    """
    qfile += self.ext
    if (self.io_type == "hdf5"):
      f = h5py.File(qfile,'w')
      self.generate_data_hdf5(f, case, self.size_v, "merid", t1, t2)
    else:
      self.generate_v_ascii(qfile, case, t1, t2)
 
  def generate_v_ascii(self, qfile, case, t1=None, t2=None):
    print("Generating %s"% qfile)
    f = open(qfile, 'w')
    coords = self.generator("merid")
    for (lat, lon) in coords:
      val = case.get_v(lat, lon, t1, t2)
      self.write_data(val, lat, lon, f)

    f.close()

  def check_v_size(self, v_in):
    """ Check the number of lines is as expected
    """
    if (self.io_type == "hdf5"):
      f = h5py.File(v_in,'r')
      nb = f["v"].size
      f.close()
    else:
      nb = self.nb_lines(v_in)

    if (nb != self.size_v):
      raise Exception("Wrong number of lines in "+u_in+"\n"
          +"Expected "+str(size)+ " got "+str(nb))
     
  def check_u_size(self, u_in):
    """ Check the number of lines is as expected
    """
    if (self.io_type == "hdf5"):
      f = h5py.File(u_in,'r')
      nb = f["u"].size
      f.close()
    else:
      nb = self.nb_lines(u_in)

    if (nb != self.size_u):
      raise Exception("Wrong number of lines in "+u_in+"\n"
          +"Expected "+str(size)+ " got "+str(nb))
     
  def correct_winds(self, u_in, u_out, v_in, v_out):
    """Correct meridional and zonal winds for null-divergence
    This is the interface, to  assure the function are called properly
    We need 4 files (2 input winds, 2 output winds)
    Winds sing is switched from positive=north to positive=south
    """
    v_in += self.ext
    v_out += self.ext
    u_in += self.ext
    u_out += self.ext
    if (self.wind_format != "borders"):
      raise Exception("winds format should be borders for correction")
    self.check_v_size(v_in)
    self.correct_v(v_in, v_out)

    self.check_u_size(u_in)
    self.correct_u(u_in, u_out, v_in, v_out)

  def sum_merid_fluxes(self, i, v_val):
    """ Sum of southern fluxes for a latitude line
    """
    flux = 0
    k_n = 0
    for j in range(1, self.nb_cells(i)+1):
      flux2, k_n = self.merid_fluxes(i, j, i+1, v_val, k_n)
      flux += flux2
    return flux


  def correct_v(self, finput, foutput):
    """ Correct meridional winds (read from a file) for a null divergence
    "" Also convert from m/s to degree/s on unit sphere
    """
    if (self.io_type == "hdf5"):
      self.correct_v_hdf5(finput)
    else:
      self.correct_v_ascii(finput, foutput)

  # Correct winds inplace
  def correct_v_hdf5(self, finput):
    print("Correcting "+finput)
    f = h5py.File(finput,'r+')
    lat = f["lat"]
    val = f["v"]

    k = 0
    # Except for last latitude line
    for i in range(1,self.nb_lat):
      n = self.nb_nodes(i)
      tmp = val[k:k+n]
      tmp *= 180./(R_earth*np.pi)
      flux = self.sum_merid_fluxes(i, tmp)
    
      # We sum all the fluxes and divides by the line surface
      dl = 360.*np.cos(lat[k]*np.pi/180.) 
      div = flux/dl
      # Then we correct the values
      # Also, we switch the wind convention > 0 means now towards south
      val[k:k+n] = - (tmp - div)
      k += n
     
    f.close()


  def correct_v_ascii(self, finput, foutput):
    print("Correcting "+finput)
    # Local arrays
    val = np.zeros(1)
    lat = np.zeros(1)
    lon = np.zeros(1)
    fo = open(foutput, 'w')
    with open(finput,'r') as fi:
      # Except for last latitude line
      for i in range(1,self.nb_lat):
        n = self.nb_nodes(i)
        for j in range(n):
          val.resize(n)
          lat.resize(n)
          lon.resize(n)
          line = self.next_valid_line(fi)
          [val[j], lat[j], lon[j]] = line.split()

        val *= 180./(R_earth*np.pi)

        flux = self.sum_merid_fluxes(i, val)
      
        # We sum all the fluxes and divides by the line surface
        dl = 360.*np.cos(lat[j]*np.pi/180.) 
        div = flux/dl
        # Then we correct the values
        # Also, we switch the wind convention > 0 means now towards south
        for j in range(n):
          fo.write(wformat % (-(val[j] - div), lat[j], lon[j]))
     
    fo.close()

  def read_array_u(self, i, fi_u):
    """ Read a latitude line and store zonal winds values, with the positions
    """
    nb_cells = self.nb_cells(i)
    val = np.zeros(nb_cells+1)
    lon = np.zeros(nb_cells+1)
    lat = np.zeros(nb_cells+1)
    for j in range(nb_cells+1):
      line = self.next_valid_line(fi_u)
      [val[j], lat[j], lon[j]] = line.split()
    return [val, lat, lon]

  def read_array_v(self, i, fi_v):
    """ Read a latitude line and store meridional winds values between i and i+1
    """
    nb_nodes = self.nb_nodes(i)
    v_val = np.zeros(nb_nodes)
    for j in range(nb_nodes):
      line = self.next_valid_line(fi_v)
    
      # We only want the winds value
      [v_val[j], _, _] = line.split()
    return v_val

  def get_dlon(self, i):
    """ Cell dlon at given latitude line
    """
    dlon = 360./self.nb_cells(i)
    return dlon

  def interface_length(self, i, j, i_n, j_n):
    """ Meridional interface length between (i,j) and (i_n, j_n) in degree
    """
    dlon = self.get_dlon(i)
    dlon_next = self.get_dlon(i_n)
    lon_max = min(j*dlon, j_n*dlon_next)
    lon_min = max((j-1)*dlon, (j_n-1)*dlon_next)

    if (i_n > i):
      lat = 90. - i*self.dlat
    else:
      lat = 90. - i_n*self.dlat

        # Testing
    #res = R*(lon_max - lon_min)*np.cos(lat*np.pi/180.)
    res = (lon_max - lon_min)*np.cos(lat*np.pi/180.)
    return res

  def sum_merid_fluxes_cell(self, i, j, v_val_prev, v_val_next, k_p, k_n):
    """ Sum of all meridional fluxes for a cell. Add inwards fluxes, substract
    outwards.
    """
    flux = 0
    # We may have None arrays, so check it
    if (v_val_prev is not None):
      flux2, k_p = self.merid_fluxes(i, j, i-1, v_val_prev, k_p)
      flux += flux2
    if (v_val_next is not None):
      flux2, k_n = self.merid_fluxes(i, j, i+1, v_val_next, k_n)
      flux -= flux2
    return flux, k_p, k_n

  def merid_fluxes(self, i, j, i_neighb, v_val, k):
    """ Computes either northern meridional fluxes or south,
    according to i_neighb, for a cell.
    k is the index in the whole winds array, we need to update it
    """
    flux = 0
    if (i_neighb > i):
      j_neighb_n = self.south_neighbours(i, j)
    else:
      j_neighb_n = self.north_neighbours(i, j)

    for l in j_neighb_n:
      interf = self.interface_length(i, j, i_neighb, l)
      #print k, v_val.size, j_neighb_n, v_val[k]
      flux += v_val[k]*interf
      k += 1
    return flux, k

  def correct_u_line_ascii(self, fi_u, i, v_val_prev, v_val_next, fo):
    """ From the previous and next meridional winds values, we correct u so that
    the sum of fluxes in a cell must be null
    We may not have a previous or next winds values, so we can have None instead
    """
    # Read array of values first
    [u_val, lat, lon] = self.read_array_u(i, fi_u)

    u_val *= 180./(R_earth*np.pi)
    # First value is unchanged
    fo.write(wformat % (u_val[0], lat[0], lon[0]))
    # We need to get the index back
    k_p = 0
    k_n = 0
    n = self.nb_cells(i)
    for j in range(1, n):
      flux, k_p, k_n = self.sum_merid_fluxes_cell(i, j, v_val_prev, 
          v_val_next, k_p, k_n)

      #u_val[j] = u_val[j-1] + flux/(R*self.dlat)
        # Testing
      u_val[j] = u_val[j-1] + flux/(self.dlat)
      fo.write(wformat % (u_val[j], lat[j], lon[j]))

    # Periodicity
    u_val[n] = u_val[0]
    fo.write(wformat % (u_val[n], lat[n], lon[n]))

  def correct_u(self, u_in, u_out, v_in, v_out):
    """ Correct zonal winds (read from a file) for a null divergence
    Must be called after correct_v
    """
    if (self.io_type == "hdf5"):
      self.correct_u_hdf5(u_in, v_in)
    else:
      self.correct_u_ascii(u_in, u_out, v_out)

  def correct_u_hdf5(self, u_in, v_in):
    # Data is written inplace
    print("Correcting "+u_in)
    v_val_next = np.zeros(1)
    v_val_prev = np.zeros(1)

    fi_u = h5py.File(u_in,'r+')
    fi_v = h5py.File(v_in,'r')

    u_val = fi_u["u"]
    u_val[:] *= 180./(R_earth*np.pi)
    v_val = fi_v["v"]

    # Initialization : first line
    n = self.nb_nodes(1)
    n2 = self.nb_cells(1)
    v_val_next = v_val[0:n]
    u_val = self.correct_u_line_hdf5(u_val, 0, n2, 1, None, v_val_next)
    
    # Main loop from 2 to n-1
    k = n
    k2 = n2+1
    for i in range(2, self.nb_lat):
      v_val_prev.resize(len(v_val_next))
      v_val_prev = np.copy(v_val_next)

      n = self.nb_nodes(i)
      n2 = self.nb_cells(i)
      v_val_next = v_val[k:k+n]
      u_val = self.correct_u_line_hdf5(u_val, k2, k2+n2, i, v_val_prev, v_val_next)
      #if (i == 2): print "final", u_val[0:3]
      k += n
      # Add one for eastmost boundary
      k2 += n2 + 1
    
    #Finalize : last line
    v_val_prev.resize(len(v_val_next))
    v_val_prev = np.copy(v_val_next)
    n2 = self.nb_cells(self.nb_lat)
    u_val = self.correct_u_line_hdf5(u_val, k2, k2+n2, self.nb_lat, v_val_prev, None)

    fi_u.close()
    fi_v.close()

  def correct_u_line_hdf5(self, u_val, start, end, i, v_val_prev, v_val_next):
    """ From the previous and next meridional winds values, we correct u so that
    the sum of fluxes in a cell must be null
    We may not have a previous or next winds values, so we can have None instead
    Here, we modify the array directly 
    Start and end define the subarray modified here
    """

    # We need to get the index back
    k_p = 0
    k_n = 0
    # First value is unchanged, last value is changed later
    n = self.nb_cells(i)
    for j in range(1, n):
      flux, k_p, k_n = self.sum_merid_fluxes_cell(i, j, v_val_prev, 
          v_val_next, k_p, k_n)

      u_val[start+j] = u_val[start+j-1] + flux/self.dlat

    # Periodicity
    u_val[end] = u_val[start]
    return u_val


  def correct_u_ascii(self, u_in, u_out, v_in):
    print("Correcting "+u_in)
    v_val_next = np.zeros(1)
    v_val_prev = np.zeros(1)

    fo = open(u_out, 'w')
    fi_u = open(u_in, 'r')
    fi_v = open(v_in, 'r')

    # Initialization : first line
    v_val_next = self.read_array_v(1, fi_v)
    self.correct_u_line_ascii(fi_u, 1, None, v_val_next, fo)
    
    # Main loop
    for i in range(2, self.nb_lat):
      v_val_prev.resize(len(v_val_next))
      v_val_prev = np.copy(v_val_next)

      v_val_next = self.read_array_v(i, fi_v)
      self.correct_u_line_ascii(fi_u, i, v_val_prev, v_val_next, fo)
    
    #Finalize : last line
    v_val_prev.resize(len(v_val_next))
    v_val_prev = np.copy(v_val_next)

    self.correct_u_line_ascii(fi_u, self.nb_lat, v_val_prev, None, fo)

    fo.close()
    fi_u.close()
    fi_v.close()


#-------------------------------------------------------------------------------
# Analytical cases : generic class
# Note : use new-class for calling parent constructor later
#-------------------------------------------------------------------------------

class AnalyticalCase(object):
  """ Define an analytical case
  * zonal : zonal advection, no mass correction
  * merid : meridional advection, no mass correction
  * solid_rot : solid rotation around an axis with mass correction
  * hourdin : analytical 2d rotation
  """
  def __init__(self, case, beta=None):
    self.name = case
    print(case)
    if (case == "zonal"):
      self.case = _ZonalCase()
    elif (case == "meridional"):
      self.case = _MeridCase()
    elif (case == "solid_rot"):
      self.case = _SolidRotCase()
    elif (case == "hourdin"):
      self.case = _HourdinCase(beta)
    elif (case == "gaussian_hills"):
      self.case = _GaussianHills()
    elif (case == "cosine_bells"):
      self.case = _CosineBells()
    elif (case == "cosine_bells_corr"):
      self.case = _CosineBellsCorrelated()
    elif (case == "custom"):
      self.case = _Custom()
    elif (case == "random"):
      self.case = _Random()
    else:
      raise Exception("Wrong case")
    self.display()

  def display(self):
    """ Print test case specifications
    """
    self.case.display();

  def set_t1(self, t):
    if (self.name == "hourdin" or self.name == "custom"):
      self.case.set_t1(t)
    else:
      raise Exception("Only for hourdin case")

  def set_t2(self, t):
    if (self.name == "hourdin"):
      self.case.set_t2(t)
    else:
      raise Exception("Only for hourdin case")


  def get_ratio(self, lat, lon):
    """ Return cell ratio at a given latitude and longitude.
    "" Input in degrees.
    """
    return self.case.get_ratio(lat, lon);

  def get_u(self, lat, lon, t1=None, t2=None):
    """ Return zonal winds at a given latitude and longitude
    """
    if (self.case == "random"):
      print "No winds for random"
    else:
      return self.case.get_u(lat, lon,t1, t2);

  def get_v(self, lat, lon, t1=None, t2=None):
    """ Return meridional winds at a given latitude and longitude
    """
    if (self.case == "random"):
      print "No winds for random"
    else:
      return self.case.get_v(lat, lon, t1, t2);


#-------------------------------------------------------------------------------
# Zonal advection of a meridional concentration band
#-------------------------------------------------------------------------------
class _ZonalCase(AnalyticalCase):
  # Made for a number of cells = 179 at the equator (so 2*90 latitudes)
  def __init__(self):
    self.width = 10.
    self.norm = 60./179.;

  def set_band(self, width):
    self.width = 10.

  def display(self):
    print("Zonal advection.")
    print("Band concentration of width=", self.width)

  def get_ratio(self, lat, lon):
    val = 0.
    if (lon > 0. and lon < self.width):
      val = 1.
    return val

  # Zonal winds depend from the latitude
  def get_u(self, lat, lon, t1=None, t2=None):
    return self.norm*np.cos(lat*coef);

  def get_v(self, lat, lon, t1=None, t2=None):
    return 0


#-------------------------------------------------------------------------------
# Meridional advection of a zonal concentration band
#-------------------------------------------------------------------------------
class _MeridCase(AnalyticalCase):
  def __init__(self):
    self.lat0 = 50.
    self.lat1 = 40.
    self.norm = 1./160.

  def set_band(self, lat0, lat1):
    if (lat0 < lat1):
      raise Exception("lat1 must be less than lat0")
    self.lat0 = lat0
    self.lat1 = lat1

  def get_ratio(self, lat, lon):
    val = 0.
    if (lat < self.lat0 and lat > self.lat1):
      val = 1.
    return val

  def display(self):
    print("Meridional advection.")
    print("Band concentration in [", self.lat0, ",", self.lat1, "]")

  def get_u(self, lat, lon, t1=None, t2=None):
    return 0

  # Constant winds
  def get_v(self, lat, lon, t1=None, t2=None):
    return self.norm

#-------------------------------------------------------------------------------
# Solid rotation around a custom axis
#-------------------------------------------------------------------------------
class _SolidRotCase(AnalyticalCase):
  # Must have lat0 > lat1
  def __init__(self):
    self.alpha = 0.
    self.beta = 0.5*np.pi
    self.norm = 0.0001*180./np.pi    
    self.lat0 = 10.
    self.lat1 = -10.
    # Used to speed up computation
    self.cos_a = np.cos(self.alpha)
    self.sin_a = np.sin(self.alpha)
    self.cos_b = np.cos(self.beta)
    self.sin_b = np.sin(self.beta)
    self.lib = cdll.LoadLibrary("./libanalytical.so")


  def set_angles(self, alpha, beta):
    self.alpha = alpha
    self.beta = beta

  def set_band(self, lat0, lat1):
    if (lat0 < lat1):
      raise Exception("lat1 must be less than lat0")
    self.lat0 = lat0
    self.lat1 = lat1

  def get_ratio(self, lat, lon):
    val = 0.
    if (lat < self.lat0 and lat > self.lat1):
      val = 1.
    return val

  def display(self):
    print("Solid rotation around the axis %f, %f" 
        % (self.alpha, self.beta))

    def get_u(self, lat, lon, t1=None, t2=None):
      return self.find_wind(lat, lon, True)

  def get_v(self, lat, lon, t1=None, t2=None):
    return self.find_wind(lat, lon, False)

  # Rotation around the z-axis, then around y'1 (euler notation)
  # Call the fortran function (~3 times faster)
  def find_wind(self, lat, lon, is_zonal):
    method = self.lib.find_wind_rotation_
    method.restype = c_double
    norm = c_double(self.norm)
    lat = c_double(lat*np.pi/180.)
    lon = c_double(lon*np.pi/180.)
    alpha = c_double(self.alpha)
    beta = c_double(self.beta)
    iszonal = c_bool(is_zonal)
    res = method(byref(norm), byref(lat), byref(lon), byref(alpha), byref(beta),
        byref(iszonal)) 
    return res


#-------------------------------------------------------------------------------
# Analytical 2d advection
#-------------------------------------------------------------------------------
class _HourdinCase(AnalyticalCase):
  # Must have lat0 > lat1
  def __init__(self, beta=None):
    self.U_0 = 229.7 #m./s and allows for a period of 2 days at the center
    self.t1 = 0
    self.t2 = 0
    # If we want to do a solid rotation of the test case, this is the rotation
    # axis, defined by Euler angles. Only works for -pi/2
    self.beta = 0 if (beta is None) else -np.pi*0.5
    #self.beta = -np.pi*0.5
    self.lib = cdll.LoadLibrary("./libanalytical.so")

  def display(self):
    print("Hourdin advection (snail)")

  # Analytical winds. Input in degrees, output in m/s
  def get_u(self, lat, lon, t1=None, t2=None):
    # Go back to find initial point (before rotation)
    lat, lon = lat*coef, lon*coef
    lat_p, lon_p =  self.new_spherical_coord(lat, lon)

    u_p = 2*self.U_0*np.cos(lat_p)*np.sin(lat_p)*np.cos(0.5*lon_p)**2
    if (self.beta != -np.pi*0.5): return u_p

    u = self.u_theta(lon_p)

    v_p = -self.U_0*np.cos(lat_p)*np.cos(0.5*lon_p)*np.sin(0.5*lon_p)
    v  = self.u_phi(lat_p, lon_p)

    x, y, z = [u_p*u[i] + v_p*v[i] for i in range(len(u))]

    # Rotate it
    x, y, z = self.new_coord(x, y, z, -self.beta)

    # Projection on u_phi
    [x_u, y_u, z_u] = self.u_theta(lon)
    u = x*x_u + y*y_u + z*z_u
    return u

  # Analytical winds. Input in degrees, output in m/s
  def get_v(self, lat, lon, t1=None, t2=None):
    # Eventual solid rotation
    lat, lon = lat*coef, lon*coef
    lat_p, lon_p =  self.new_spherical_coord(lat, lon)

    v_p = -self.U_0*np.cos(lat_p)*np.cos(0.5*lon_p)*np.sin(0.5*lon_p)
    # No rotation
    if (self.beta != -np.pi*0.5): return v_p

    v  = self.u_phi(lat_p, lon_p)

    # Get 2d initial vector
    u_p = 2*self.U_0*np.cos(lat_p)*np.sin(lat_p)*np.cos(0.5*lon_p)**2
    u = self.u_theta(lon_p)


    x, y, z = [u_p*u[i] + v_p*v[i] for i in range(len(u))]

    # Rotate it
    x, y, z = self.new_coord(x, y, z, -self.beta)

    # Projection on u_phi
    [x_v, y_v, z_v] = self.u_phi(lat, lon)
    v = x*x_v + y*y_v + z*z_v
    return v

  def set_t1(self, t):
    """ Update starting time in analytical version 
    """
    self.t1 = t

  def set_t2(self, t):
    """ Update final time in analytical version 
    """
    self.t2 = t

  def get_ratio(self, lat, lon):
    # Eventual solid rotation
    lat, lon =  self.new_spherical_coord(lat*coef, lon*coef)
    lat_prev, lon_prev = self.find_prev_lat_lon(lat, lon)

    q_prev = self.get_gaussian_ratio(lat_prev, lon_prev)
    return q_prev

  # Return gaussian concentration at given cell. Input in degree
  def get_gaussian_ratio(self, lat, lon):

    # Switch from [0,2pi] to [-pi,pi]
    #lon2 = lon2 / coef
    #if (lon2 > 180): lon2 = lon2 - 360
    if (lon > 180): lon = lon - 360
    
    return np.exp(-lon**2/(2*70**2))

  # Input and output in [0, 2pi]
  # Interface to fortran
  def new_spherical_coord(self, lat, lon, beta0=None):
    method = self.lib.new_spherical_coordinates_
    method.restype = c_double
    lat, lon = c_double(lat), c_double(lon)
    alpha = c_double(0)
    beta = c_double(self.beta) if (beta0 is None) else c_double(beta0)
    lat_rot, lon_rot = c_double(0), c_double(0)

    res = method(byref(lat_rot), byref(lon_rot), byref(lat), byref(lon), 
        byref(alpha), byref(beta))
    return lat_rot.value, lon_rot.value

  # Returns a triplet for unitary vector
  def u_theta(self, lon):
    return self.cartesian_coord(0., lon + np.pi*0.5)

  # Returns a triplet for unitary vector. Phi is the latitude.
  def u_phi(self, lat, lon):
    return self.cartesian_coord(lat + np.pi*0.5, lon)


  # Interface to fortran
  def cartesian_coord(self, lat, lon):
    method = self.lib.cartesian_coordinates_
    method.restype = c_int
    x, y, z = c_double(0), c_double(0), c_double(0)
    lat, lon = c_double(lat), c_double(lon)

    res = method(byref(x), byref(y), byref(z), byref(lat), byref(lon))
    return x.value, y.value, z.value

  # Interface to fortran
  def new_coord(self, x, y, z, beta):
    method = self.lib.new_coordinates_
    method.restype = c_int
    x1, y1, z1 = c_double(0), c_double(0), c_double(0)
    x, y, z = c_double(x), c_double(y), c_double(z)
    beta = c_double(beta)

    res = method(byref(x1), byref(y1), byref(z1), byref(x), byref(y), byref(z),
        byref(c_double(0)), byref(beta))
    return x1.value, y1.value, z1.value

  # Find the psi for the current position. Input in radian
  def find_psi(self, lat, lon):
    res = R_earth*self.U_0*np.cos(lon*0.5)**2
    res = res*np.cos(lat)**2
    return res

  # Find where the tracer in the current cell came from. For that, we go back on
  # a trajectory (psy constant). Input and output in radian
  def find_prev_lat_lon(self, lat_r, lon_r):
    psi = self.find_psi(lat_r,lon_r)
    lat_prev_r, value, alpha, beta = self.find_prev_lat(psi, lat_r, lon_r)
    half = self.find_half(lat_r, lon_r, value, alpha, beta, psi)
    lon_prev_r = self.find_prev_lon(half, lat_prev_r, lat_r, lon_r)

    lat_prev = lat_prev_r*180./np.pi
    lon_prev = lon_prev_r*180./np.pi
    return lat_prev, lon_prev

  # Arcsin is defined in [-pi/2,pi/2] so we need to know where the starting 
  # point is :
  # half = 1 : first half (longitude < 180)
  # half = -1 : second half (longitude > 180)
  def find_half(self, lat, lon, value, alpha, beta, psi):
    zone = 1.
    if (lon > np.pi):
      zone = -1.
    frac = np.sin(lat)/beta
    if (frac > 1-precision):
      frac = 1.
    if (frac < -1+precision):
      frac = -1.

    delta_t = np.arcsin(zone*value) - np.arcsin(frac)
  
    delta_t = abs(delta_t)/alpha
    half = 1.
    # Takes the modulo
    period = 2*np.pi/alpha
    diff = (self.t2-self.t1) % period
    # Check if we go to the other half, without going back on the cur_rent half
    if (diff > delta_t and diff < period*0.5 + delta_t):
      half = -1.
    if (lon > np.pi):
      half = -half
    return half

  # Previous latitude, at time t1 knowing the position at time t2. Input in
  # radians.
  # We solve arcsin(sin psi1/beta) = arcsin( sin psi2/beta) + alpha*sign*(t2-1)
  # Contrary to what Hourdin says in its paper, it's not the sign of lon but
  # rather the sign of cos(0.5*lon)*cos(lat). 
  def find_prev_lat(self, psi, lat, lon):
    a = R_earth
    beta = np.sqrt(1 - psi/(a*self.U_0))
    alpha = np.sqrt(psi*self.U_0/(a**3))
    value = np.sin(lat)/beta 

    #Restrict to [-1,1]
    value = min(value, 1.)
    value = max(value, -1.)
    arcsin_tmp = np.arcsin(value)

    # Need longitude sign in equation
    sign = -1 if (np.cos(0.5*lon)*np.cos(lat) < 0) else 1
    #Should be in [-1,1]
    sin_tmp = beta*np.sin(arcsin_tmp + alpha*sign*(self.t2-self.t1))
    lat = np.arcsin(sin_tmp)

    return lat, value, alpha, beta

  # Previous longitude, at time t1 knowing the position at time t2. Input in
  # radians.
  def find_prev_lon(self, half, lat_prev, lat, lon):
    reduced = lon
    #print half
    #if (lon > np.pi):
    #  reduced = 2*np.pi - lon
  
    delta = np.cos(reduced*0.5)*np.cos(lat)/np.cos(lat_prev)
    if (delta > 1-precision):
      delta = 1.
    if (delta < -1+precision):
      delta = -1.
    lon_prev = 2*np.arccos(delta)
    #if (half == -1):
    #    lon_prev = 2*np.pi - lon_prev
    return lon_prev

#-------------------------------------------------------------------------------
# Test suite of Lauritzen et al. Contains common data
#-------------------------------------------------------------------------------
class _TestSuite(AnalyticalCase):
  def __init__(self):
    self.lon1 = 5*np.pi/6 #+ 2*np.pi/3
    self.lat1 = 0 
    self.lon2 = 7*np.pi/6 #+ 2*np.pi/3
    self.lat2 = 0
    self.T = 12*24*3600 # seconds 

  # Analytical winds. Input in degrees.
  def get_u(self, lat, lon, t1, t2):
    lat = lat*coef
    lon2 = lon*coef - 2*np.pi*t1/self.T
    coef1 = 10*R/self.T;
    coef2 = 2*np.pi*R/self.T;
    # Winds are taken at half the timestep
    cost = np.cos(np.pi*0.5*(t1+t2)/self.T)
    res = coef1*np.sin(lon2)**2 * np.sin(2*lat)*cost + coef2*np.cos(lat)

    #if (lat > 89*coef):
      #print "u", res, lat, lon, t1, t2
    return res

  # Analytical winds. Input in degrees.
  def get_v(self, lat, lon, t1, t2):
    lat = lat*coef
    lon2 = lon*coef - 2*np.pi*t1/self.T
    # Change the sign so positive winds mean advection to the south
    coef1 = 10*R/self.T;
    cost = np.cos(np.pi*0.5*(t1+t2)/self.T)
    res = coef1* np.sin(2*lon2)*np.cos(lat)*cost
    return res

#-------------------------------------------------------------------------------
# Gaussian hills case, subclass of _TestSuiteCase
#-------------------------------------------------------------------------------
class _GaussianHills(_TestSuite):

  def __init__(self):
    # Call parent's constructor
    super(_GaussianHills,self).__init__()

    self.set_hills_centers()
    self.h_max = 0.95
    self.b = 5

  def set_hills_centers(self):
    self.x1, self.y1, self.z1 = self.to_cartesian(self.lat1, self.lon1)
    self.x2, self.y2, self.z2 = self.to_cartesian(self.lat2, self.lon2)

  # Input in radians. Convert to an unit sphere
  def to_cartesian(self, lat, lon):
    return np.cos(lon)*np.cos(lat), np.sin(lon)*np.cos(lat), \
        np.sin(lat)


  # Ratio for one of the two gaussian hills (given by the id)
  # Input in radians
  def gaussianhill(self, lat, lon, trolol):
    x, y , z = self.to_cartesian(lat, lon)
    if (trolol == 1):
      val = (x - self.x1)**2 + (y - self.y1)**2 + (z - self.z1)**2
    else:
      val = (x - self.x2)**2 + (y - self.y2)**2 + (z - self.z2)**2

    return self.h_max*np.exp(-self.b*val)

  def display(self): print("Gaussian hills")

  def get_ratio(self, lat, lon):
    lat = lat*coef
    lon = lon*coef

    res = self.gaussianhill(lat, lon, 1) + self.gaussianhill(lat, lon, 2)
    return res

#-------------------------------------------------------------------------------
# Cosine bells hills case, subclass of _TestSuiteCase
#-------------------------------------------------------------------------------
class _CosineBells(_TestSuite):

  def __init__(self):
    # Call parent's constructor
    super(_CosineBells,self).__init__()
    self.r = R_earth/2
    self.b = 0.1
    self.c = 0.9
    self.h_max = 1

  # Great-circle distance. Input in radian
  def great_circle_dist(self, lat, lon, case):
    lat_tmp, lon_tmp = self.lat1, self.lon1
    if (case == 2):
      lat_tmp, lon_tmp = self.lat2, self.lon2

    res = np.sin(lat_tmp)*np.sin(lat) + \
          np.cos(lat_tmp)*np.cos(lat)*np.cos(lon-lon_tmp)
    return R_earth*np.arccos(res)

  # Create 1 cosine bell. Input in radian
  def bell(self, r_i):
    res = 0.5*self.h_max*(1+np.cos(np.pi*r_i/self.r))
    return res

  def display(self): print("Cosine Bells")

  # Create 2 cosine bells. Input in radian
  def get_ratio(self, lat, lon):
    lat = lat*coef
    lon = lon*coef

    r1 = self.great_circle_dist(lat, lon, 1)
    r2 = self.great_circle_dist(lat, lon, 2)

    if (r1 < self.r):
      res = self.b + self.c*self.bell(r1)
    elif (r2 < self.r):
      res = self.b + self.c*self.bell(r2)
    else:
      res = self.b
    return res

#-------------------------------------------------------------------------------
# Create a distribution correlated to CosineBells
#-------------------------------------------------------------------------------
class _CosineBellsCorrelated(_TestSuite):

  def __init__(self):
    # Call parent's constructor
    super(_CosineBellsCorrelated,self).__init__()
    self.bells = _CosineBells()
    self.a1 = -0.8
    self.b1 = 0.9

  def display(self): print("Cosine Bells correlated")

  # Create 2 cosine bells. Input in radian
  def get_ratio(self, lat, lon):
    # No need to cast to radians
    res = self.bells.get_ratio(lat, lon)
    res = self.a1*res**2 + self.b1

    return res

#-------------------------------------------------------------------------------
# For I/O tests : random data
#-------------------------------------------------------------------------------
class _Random(_TestSuite):

  def __init__(self):
    # Call parent's constructor
    super(_Random,self).__init__()

  def display(self): print("Random")

  def get_ratio(self, lat, lon):
    return random.random()

#-------------------------------------------------------------------------------
# For testing : v = 0, u = constant, and cos(lon) as concentration
#-------------------------------------------------------------------------------
class _Custom(_TestSuite):

  def __init__(self):
    # Call parent's constructor
    super(_Custom,self).__init__()
    # Winds are set such as 1 min empties half a cell witth 80 lat2
    self.norm = 2*np.pi*R_earth/(3*159) # dx
    self.norm = 0.5*self.norm/60.

  def display(self): print("Custom")

  def get_u(self, lat, lon, t1, t2):
    return self.norm

  def get_v(self, lat, lon, t1, t2):
    return 0
 
  # One cosine bells
  def get_ratio(self, lat, lon):
    return abs(np.cos(lon*coef))

