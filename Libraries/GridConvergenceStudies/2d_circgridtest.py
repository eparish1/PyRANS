from pylab import *
import sys
sys.path.append("../../Libraries/")
from libraries import *
import time
close("all")

#####################################################################
## This is a script to determine the order of accuracy of a
## 2D mesh around a circle. This will test the grid metrics 
## algorithm for a 2D case. The field that we compare derivatives
## of is 
## u = sin(x) + sin(y)
#####################################################################
## Create grid
t1 = time.time()
nruns = 8
x_error = zeros(nruns)
y_error = zeros(nruns)
dxm = zeros(nruns)
dym = zeros(nruns)
nxa = zeros(nruns)
nya = zeros(nruns)
nxa[:] = [11,31,71,101,301,501,1001,2001]
nya[:] = [7,27,57,87,251,451,851,1551]
for L in range(0,int(nruns)):
  print(L)
  #xperd = 0
  #yperd = 1
  nr = int(nxa[L])
  ntheta = int(nya[L])
  Rg = 40.**(1./nr)
  r = zeros(nr)
  r[0] = 1.
  for i in range(0,nr-1):
      r[i+1] = r[i] + r[i]*(Rg - 1)
  #r = linspace(0.1,3,nr)
  nz = 1
  theta = linspace(0,2*pi,ntheta)
  x = zeros((nr,ntheta,nz))
  y = zeros((nr,ntheta,nz))
  z = zeros((nr,ntheta,nz))
  for i in range(0,nr):
      for j in range(0,ntheta):
         for k in range(0,nz):
              x[i,j,k] = r[i]*cos(theta[j])
              y[i,j,k] = r[i]*sin(theta[j])
              z[i,j,k] = k


  ## Dummy field
  class grid():
      A = gridMetrics(x,y,z)
      xperd = 0
      yperd = 0
  u = 10*sin(x) + 4*sin(y)
  up_x = 10*cos(x)
  up_y = 4*cos(y)
  gradU = d_dX(u,grid)
  x_error[L] = mean(abs(( (gradU[:,:,:,0] - up_x)) ) )
  y_error[L] = mean(abs(( (gradU[:,:,:,1] - up_y)) ) )
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
xlabel(r'$dx,dy$')
ylabel(r'mean(|E|)')
show()