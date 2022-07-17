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
# Read a file on format (int, nb_float*float) and return them
# It is faster to operate on an array
#-------------------------------------------------------------------------------
def read_data(filename, nb_floats):
  nb_lines = file_len(filename)
  if (nb_floats > 2):
    print "No more than 2 float columnms"
    exit

  f = open(filename, 'r')
  val1 = np.zeros(nb_lines,float)
  if (nb_floats > 1):
    val2 = np.zeros(nb_lines,float)

  for k in range(nb_lines):
    line = f.readline()
    if (not line): break
    cur = line.split()
    val1[k] = float(cur[1])
    if (nb_floats > 1):
      val2[k] = float(cur[2])
 
  f.close()
  if (nb_floats > 1):
    return val1, val2
  else:
    return val1

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
  mean_time[k-1] = sum(time)/len(time)

###############################################################################

nb_procs_max = 3#64
max_time = np.zeros(nb_procs_max)
min_time = np.zeros(nb_procs_max)
mean_time = np.zeros(nb_procs_max)

max_time_t = np.zeros(nb_procs_max)
min_time_t = np.zeros(nb_procs_max)
mean_time_t = np.zeros(nb_procs_max)

max_time_a = np.zeros(nb_procs_max)
min_time_a = np.zeros(nb_procs_max)
mean_time_a = np.zeros(nb_procs_max)

max_time_c = np.zeros(nb_procs_max)
min_time_c = np.zeros(nb_procs_max)
mean_time_c = np.zeros(nb_procs_max)

nb_procs = (180, 360, 540)#range(1, nb_procs_max+1)
print nb_procs

k = 0
for l in nb_procs:
  # Read total time and simulation time
  filename = "../output_grid/totaltime_64_"+str(l)
  time_t, time = read_data(filename, 2)
  max_min_mean(time_t, k, max_time_t, min_time_t, mean_time_t)

  max_min_mean(time, k, max_time, min_time, mean_time)

  # Read advection computing time
  filename = "../output_grid/computingtime_64_"+str(l)
  time_a = read_data(filename, 1)
  max_min_mean(time_a, k, max_time_a, min_time_a, mean_time_a)

  # Communication time
  time_c = time - time_a
  max_min_mean(time_c, k, max_time_c, min_time_c, mean_time_c)
  k = k + 1


fig = plt.figure()

ax1 = fig.add_subplot(221)
plot_data(ax1, nb_procs, max_time_t, min_time_t, mean_time_t)
ax1.set_xlabel("Number of process")
ax1.set_ylabel("Total time")
ax1.legend()#loc=4)

ax2 = fig.add_subplot(222)
plot_data(ax2, nb_procs, max_time, min_time, mean_time)
ax2.set_xlabel("Number of process")
ax2.set_ylabel("Advection time")
ax2.legend()#loc=2)

#plt.savefig("simutime_64.png")

ax3 = fig.add_subplot(223)
plot_data(ax3, nb_procs, max_time_a, min_time_a, mean_time_a)
ax3.set_xlabel("Number of process")
ax3.set_ylabel("Advection computing time")
ax3.legend()#loc=4)

ax4 = fig.add_subplot(224)
plot_data(ax4, nb_procs, max_time_c, min_time_c, mean_time_c)
ax4.set_xlabel("Number of process")
ax4.set_ylabel("Advection communication time")
ax4.legend()#loc=4)

fig.suptitle("Process times, 1000 iterations")

#print time_c[17], time_a[17], time[17]
plt.savefig("process_times_64.png", dpi=(400))
plt.show()
