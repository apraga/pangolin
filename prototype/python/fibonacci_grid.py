###############################################################################
# Fibonacci area-preserving grid 
###############################################################################
import numpy as np
import pylab as p

# Approximate cos
def spiral(i):
  r0 = 1.
  invphi = 2./(1+np.sqrt(5))
  r = r0*np.sqrt(i-0.5)
  theta = 2*np.pi*i*invphi
  x = r*np.cos(theta)
  y = r*np.sin(theta)
  return x,y

def fibo(ref):
  prev_prev = 0
  prev = 1
  ref[0] = prev_prev
  ref[1] = prev
  for i in range(2,len(ref)):
    ref[i] = prev + prev_prev
    prev_prev = prev
    prev = ref[i]

n = 5000
x = np.zeros(n+1,float)
y = np.zeros(n+1,float)
color = np.zeros(n+1,float)
count = 0
n2 = 20
ref = np.zeros(n2+1,int)
fibo(ref)
start = 8
end = 15

for i in range(1,n+1):
  x[i],y[i] = spiral(i)
  for j in range(start,end+1):
    if (i % ref[j] == 0):
      color[i] = (j-start+2)*10
      break
  else:
    color[i] = 10

p.scatter(x, y,c=color)

#xl = p.xlabel('Latitude $\Phi$')
#yl = p.ylabel('n')
#ttl = p.title('Number of cells with $\Delta \lambda=360$')
#p.legend(bbox_to_anchor=(1., 0.3))
#p.ylim([0,1.2*n_ref])

#p.savefig('nb_cells.jpg')
p.show()
