import subprocess as sub
import numpy as np
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------
# Get number of line in file
#-------------------------------------------------------------------------------
def file_len(fname):
  p = sub.Popen(['wc', '-l', fname], stdout=sub.PIPE, 
    stderr=sub.PIPE)
  result, err = p.communicate()
  if p.returncode != 0:
    raise IOError(err)
  return int(result.strip().split()[0])

#-------------------------------------------------------------------------------
# Read data. First argument is the filename, others are the number of columns, 
# and the type of each (1 for int, float otherwise)
# We ignore the first int
#-------------------------------------------------------------------------------
def read_data(*args):
  nb_lines = file_len(args[0])

  f = open(args[0], 'r')
  val = []
  size =len(args)
  # Ignore first arg and second
  for i in range(2, size):
    if (args[i] == 1):
      val.append(np.zeros(nb_lines,int))
    else:
      val.append(np.zeros(nb_lines,float))

  for k in range(nb_lines):
    line = f.readline()
    if (not line): break
    cur = line.split()

    for i in range(2, size):
      val[i-2][k] = cur[i-1]
      #if (args[i] == 1):
        #val[i-2][k] = int(cur[i-1])
      #else:
        #val[i-2][k] = float(cur[i-1])

  f.close()
  return val

#-------------------------------------------------------------------------------
# Plot max, min, men
#-------------------------------------------------------------------------------
def plot_data(ax1, nb_procs, maxt, mint, meant):
  ax1.plot(nb_procs, maxt, label="Max")
  ax1.plot(nb_procs, mint, label="Min")
  ax1.plot(nb_procs, meant, label="Mean")

#-------------------------------------------------------------------------------
# Compute and set into the array max, min mean time
#-------------------------------------------------------------------------------
def max_min_mean(time, k, max_time, min_time, mean_time):
  max_time[k-1] = max(time)
  min_time[k-1] = min(time)
  mean_time[k-1] = sum(time)/k

#-------------------------------------------------------------------------------
# Process different times. Returns array : total, simulation, advection
# computing, communication time, chemistry time
#-------------------------------------------------------------------------------
def process_times(add_chemistry, nb_procs_max):
  nb = 4
  if (add_chemistry):
    nb = 5
  max_time = np.zeros((nb, nb_procs_max))
  min_time = np.zeros((nb, nb_procs_max))
  mean_time = np.zeros((nb, nb_procs_max))
  
  nb_cells = np.zeros(nb_procs_max)
  
  nb_procs = range(1, nb_procs_max+1)
  # Time for chemistry in mocage
  chemistry = 0.0125
  # For all cells
  chem_all = np.zeros(nb_procs_max)
  
  for k in nb_procs:
    # Read total time and simulation time
    filename = "../output/totaltime_"+str(k)
    # Read rank, total time, simulation time
    time_t, time = read_data(filename, 1, 0, 0)
    
    # Read advection computing time
    filename = "../output/computingtime_"+str(k)
    # Read rank, advection time, number of cells
    time_a, tmp = read_data(filename, 1, 0, 1)
  
    if (add_chemistry):
      # Each cell has a chemistry step
      chem_all = chemistry*tmp
      # We add it to the total time and simulation time
      time_t = time_t + chem_all
      time = time + chem_all
 
    # Total time
    max_min_mean(time_t, k, max_time[0], min_time[0], mean_time[0])
    # Simulation time
    max_min_mean(time, k, max_time[1], min_time[1], mean_time[1])
    # Advection computing time
    max_min_mean(time_a, k, max_time[2], min_time[2], mean_time[2])

    # Indice for the max, gives the number of cells
    ipos = np.where(time_a == max_time[2][k-1])
    nb_cells[k-1] = tmp[ipos[0][0]]
  
    # Deduce communication time
    time_c = time - time_a
  
    # Don't forget to substract chemistry here too
    if (add_chemistry):
      time_c = time_c - chem_all
  
    max_min_mean(time_c, k, max_time[3], min_time[3], mean_time[3])

    if (add_chemistry):
      max_min_mean(chem_all, k, max_time[4], min_time[4], mean_time[4])
    #if  (min_time_c[k-1] < 0):
      #print time, time_a, time - time_a

  return max_time, min_time, mean_time, nb_cells, nb_procs

#-------------------------------------------------------------------------------
# Plot different time. Input = total, simulation, advection computing, comm
# times
#-------------------------------------------------------------------------------
def multi_plot_graph(max_time, min_time, mean_time, nb_cells, nb_procs):

  figsize = (8*1.2,6)
  fig = plt.figure(1, figsize=figsize)

#  ax1 = fig.add_subplot(221)
#  plot_data(ax1, nb_procs, max_time[0], min_time[0], mean_time[0])
#  ax1.set_xlabel("Number of process")
#  ax1.set_ylabel("Total time")
#  ax1.legend()#loc=4)
#
#  ax2 = fig.add_subplot(222)
#  plot_data(ax2, nb_procs, max_time[1], min_time[1], mean_time[1])
#  ax2.set_xlabel("Number of process")
#  ax2.set_ylabel("Advection time")
#  ax2.legend()#loc=2)

  #plt.savefig("simutime_64.png")

  #ax3 = fig.add_subplot(223)
  ax3 = fig.add_subplot(211)
  plot_data(ax3, nb_procs, max_time[2], min_time[2], mean_time[2])
  ax3.set_xlabel("Number of process")
  ax3.set_ylabel("Computing time")
  ax3.legend()#loc=4)

  #ax4 = fig.add_subplot(224)
  ax4 = fig.add_subplot(212)
  plot_data(ax4, nb_procs, max_time[3], min_time[3], mean_time[3])
  ax4.set_xlabel("Number of process")
  ax4.set_ylabel("Communication time")
  ax4.legend()#loc=4)

  fig.suptitle("Advection time, 1000 iterations")

  plt.savefig("process_times_64.png")

  # Time versus number of cells 
  fig = plt.figure(2)
  m, b = np.polyfit(nb_cells, max_time[2], 1)
  plt.plot(nb_cells, m*nb_cells+b, color='red')
  plt.scatter(nb_cells, max_time[2])
  print m, b
  plt.xlabel("Nb cells")
  plt.ylabel("Advection computing time")
  plt.savefig("times_vs_nb_cells.png")
  plt.show()

#-------------------------------------------------------------------------------
# Plot bar chart (stacked)
#-------------------------------------------------------------------------------
def plot_barchart(max_time, min_time, mean_time, nb_procs, add_chemistry, ylog):
  nb_procs_max = len(max_time[0])

  # Comm time at the bottom
  plt.bar(nb_procs, height=max_time[3], color='red', label="Communication",
      log=ylog)
  bottoms = max_time[3]
  # Chemistry
  if (add_chemistry):
    plt.bar(nb_procs, bottom = bottoms, height=max_time[4], color='green',
        label="Chemistry", log=ylog)
    bottoms = bottoms + max_time[4]

  # Advection computing time
  plt.bar(nb_procs, bottom = bottoms, height=max_time[2], color='blue',
  label="Computing advec", log=ylog)
  plt.xlabel("Nb partitions")
  plt.ylabel("Max time (s)")
  plt.legend()
  if (add_chemistry):
    plt.savefig("bar_chemistry.png")
  else:
    plt.savefig("bar.png")
  plt.show()
################################################################################

add_chemistry = True
#add_chemistry = None
nb_procs_max = 64


max_time, min_time, mean_time, nb_cells, nb_procs = process_times(add_chemistry, nb_procs_max)
ylog = False
multi_plot_graph(max_time, min_time, mean_time, nb_cells, nb_procs)
#plot_barchart(max_time, min_time, mean_time, nb_procs, add_chemistry, ylog)

#mu, sigma = 200, 25
## create a new data-set
#x = mu + sigma*np.random.randn(1000,3)
#

