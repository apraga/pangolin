from pylab import *
import error
import numpy as np

# Lists are given as input, folder too
def compute_cv_rate(y2, yoo, folder):
  for i in [20, 40, 80, 160, 320]:
    ref = folder+"ratio_"+str(i)+"lat_0.dat"
    num = folder+"ratio_"+str(i)+"lat_T.dat"
    y2.append(error.get_error(ref,num,'2'))
    yoo.append(error.get_error(ref,num,'oo'))

recompute = False
if (recompute):
  y2 = []
  yoo = []
  compute_cv_rate(y2, yoo, "/wkdir/pae2/praga/cv_rate_limited/")

  y2_un = []
  yoo_un = []
  compute_cv_rate(y2_un, yoo_un, "/wkdir/pae2/praga/cv_rate_unlimited/")

  print y2
  print yoo
  print y2_un
  print yoo_un
else:
  y2 = [0.544993590304237, 0.35845506972377916, 0.20504514119456335, \
      0.10960394644940914, 0.057931199147530976]
  yoo = [0.61266265674649567, 0.4449382984258699, 0.28655756752829564, \
      0.16874888109052535, 0.093743047032628812]

  y2_un = [0.51406095934337126, 0.34039150800595641, 0.19774043956352386, \
      0.10817778583935432, 0.057774599215461457]
  yoo_un = [0.57817441373082856, 0.42200716226272711, 0.27339299327253436, \
      0.16394648910520548, 0.092522276803568115]

x0 = [3,1.5,0.75,0.375,0.1875]
# Log scale for x
x = np.log(x0)

first = []
sec = []
for cur in x:
  first.append(np.exp(cur+ 1))
  sec.append(np.exp(2*cur- 3))

ax = subplot(1,1,1)
ax.set_yscale('log')

plot(x,y2, 'b-', label="$l_2 $ limited", marker='v')
plot(x,yoo, 'r-', label="$l_{oo} $ limited", marker='^')
plot(x,y2_un, 'b--', label="$l_2 $ unilimited", marker='v')
plot(x,yoo_un, 'r--',label="$l_{oo}$ unilimited", marker='^')

plot(x,first, 'g')#, label="1st order")
plot(x,sec, 'black')#, label="2nd order")

# Revert ordering and add some offset
offset_l = np.log(x0[0] + 0.4)
offset_r = np.log(x0[len(x)-1] - 0.04)
xlim([offset_l, offset_r])
ylim([10e-6,10])
xticks(x,x0)

xlabel("$\Delta \lambda$")
ylabel("error")
title("Convergence rate")
legend(loc="lower left")
savefig("cv_rate.pdf")
show()
