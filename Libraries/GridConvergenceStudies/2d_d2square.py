from pylab import *
import sys
sys.path.append("../")
import time
from libraries import *
close("all")
#####################################################################
## This is a script to determine the order of accuracy of a
## 2D square mesh. This will test the differencing
## algorithm for a 2D case. The field that we compare derivatives
## of is 
## u = sin(x) + sin(y)
#####################################################################
## Create grid

## Create grid
t1 = time.time()
nruns = 5
x_error = zeros(nruns)
y_error = zeros(nruns)
dxm = zeros(nruns)
dym = zeros(nruns)
nxa = zeros(nruns)
nya = zeros(nruns)
nxa[:] = [11,31,71,101,301]
nya[:] = [7,27,57,87,251]
for L in range(0,int(nruns)):
  print(L)
  #xperd = 0
  #yperd = 1
  nr = int(nxa[L])
  ntheta = int(nya[L])
  x = linspace(0,1,nxa[L])
  y = linspace(0,1,nya[L])
  z = linspace(0,1,1)
  y,x,z = meshgrid(y,x,z)

  ## Dummy field
  class grid():
      A = gridMetrics(x,y,z)
      xperd = 0
      yperd = 0
  u = 10*sin(x) + 4*sin(y)
  up_x = -10*sin(x)
  up_y = -4*sin(y)
  gradU = d_dX(d_dX(u,grid)[:,:,:,0],grid)
  gradV = d_dX(d_dX(u,grid)[:,:,:,1],grid)
  x_error[L] = mean(abs(( (gradU[:,:,:,0] - up_x)) ) )
  y_error[L] = mean(abs(( (gradV[:,:,:,1] - up_y)) ) )
  dxm[L] = mean(abs(x[1::,:,:] - x[0:-1,:,:]))
  dym[L] = mean(abs(y[:,1::,:] - y[:,0:-1,:]))

t2 = time.time()
print('time = ' + str(t2 - t1))
plot(dxm,x_error,label = 'x Error')
plot(dxm,dxm**2,'--',label = '$\mathcal{O}(2)$')
plot(dym,y_error,label = 'y Error')
plot(dym,dym**2,'--',label = '$\mathcal{O}(2)$')                                 
xscale('log')
yscale('log')
legend(loc=2)
xlabel(r'$dx$')
ylabel(r'mean(|E|)')
show()