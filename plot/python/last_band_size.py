import pylab as plt
import numpy as np

# Number of partitions on the band (except from last)
def get_nb_part_on_lat(lat):
  return 2*lat - 1

# Total Number of partitions 
def get_nb_part_total(lat):
  return lat**2

nb_lat = 90
n_max = 200
n_start = 10
n = range(n_start, n_max + 1)

v_square = np.zeros(len(n))
v_rest_tri = np.zeros(len(n))
v_rest_rect = np.zeros(len(n))
i = 0
for nb in n:
  # Number of bands with squared partitions
  nb_bands = int(np.sqrt(nb))
  if (nb_bands**2 >= nb):
    # Square side
    side = nb_lat / nb_bands
    # Width for the last band
    width_rect = 0
    width_tri = 0
    height = 0
  else:
    side = (nb_lat - 1) / nb_bands
    # Height for last band
    height = nb_lat - nb_bands*side

    # Number of partitions on last band
    nb_part_prev = get_nb_part_on_lat(nb_bands)
    nb_tot_prev = get_nb_part_total(nb_bands)
    nb_part_last  = nb - nb_tot_prev
    #print nb, nb_part_last

    # Width for the central part (a triangle cut)
    width_tri = get_nb_part_total(nb_bands + height) - get_nb_part_total(nb_bands) 
    if (nb_part_prev >= nb_part_last):
      q = nb_part_prev / nb_part_last 
    # nb prev + 1 = nb last
    else:
      q = nb_part_prev / (nb_part_last -1)
      # smaller half
      width_tri = int(width_tri*0.5)

    # Width for the rectangular partitions
    width_rect = q*side
    #width = min(q*side, size_mid)
  v_rest_rect[i] = width_rect*height
  v_rest_tri[i] = width_tri*height
  v_square[i] = side**2

    
  i = i + 1

plt.plot(n, v_rest_rect)
plt.plot(n, v_rest_tri)
plt.plot(n, v_square)
plt.show()
#i = 195
#print "ratio 195", v_rest[i]/v_square[i]
for i in range(len(n)):
  if (n[i] == 101):
    print "rect, tri", v_rest_rect[i], v_rest_tri[i] 
  #r = v_rest[i]/v_square[i]
  #if (r > 0.9):
    #print "sup", n[i], r
