import pylab as plt
import numpy as np

def plot_vol_max_last(n_min, has_legend=None, has_xlabel=None):
  n = []
  vol_max = []
  vol_max_last = []
  f = open("../output/area.dat")
  for line in f.readlines():
    cur = line.split()
    if (int(cur[0]) > n_min):
      n.append(cur[0])
      vol_max.append(cur[1])
      vol_max_last.append(cur[2])
  f.close()
  amax, = plt.plot(n, vol_max)
  amax_last, = plt.plot(n, vol_max_last)
  if (has_legend):
    plt.legend([amax, amax_last] ,["Maximal volume", "Max. volume on last band"])
  if (has_xlabel):
    plt.xlabel("nb partitions")
  plt.ylabel("nb cells")

plt.subplot(211)
plt.title("Comparing last band load (limit = 0.99)")
plot_vol_max_last(10, has_legend=True)
plt.subplot(212)
plot_vol_max_last(100, has_xlabel=True)
plt.savefig("vol_max_last.png")
plt.show()
