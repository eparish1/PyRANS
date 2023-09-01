from pylab import *
import sys
sys.path.append("../")
import time
from libraries import *
close("all")
###############
t1 = time.time()
from mpi4py import *
n = 13
print(2**n)
x = linspace(0,1,2**11)
y = linspace(0,1,2**9)
z = linspace(0,1,1)
y,x,z = meshgrid(y,x,z)
nx,ny,nz = shape(x)
class grid():
    A = gridMetrics(x,y,z)
    xperd = 0
    yperd = 0
u = 10*sin(x) + 4*sin(y)
up_x = 10*cos(x)
up_y = 4*cos(y)
gradU = d_dX(u,grid)
print(norm(gradU[:,:,:,0] - 10*cos(x)))
t2 = time.time()
print(t2 - t1)
#print('From Proc ' + str(MPI_RANK), 'x - ', str(xs), ' ' ,str(xe))
