from pylab import *
import sys
sys.path.append("../")
import time
from libraries import *
close("all")
###############
from mpi4py import *
n = 10


##############################################
################# CODE TO BREAK UP A MESH AND 
################# DISTRIBUTE TO PROCS 
comm = MPI.COMM_WORLD
MPI_RANK = comm.Get_rank()
size = comm.Get_size()
nmx = 1
nmy = 1
size = nmx*nmy
if (MPI_RANK == 0):
  t1 = time.time()
  #################
  overlap = 1
  x = linspace(0,1,1000)
  y = linspace(0,1,1000)
  z = linspace(0,1,1)
  y,x,z = meshgrid(y,x,z)
  nx,ny,nz = shape(x)
  uG = zeros((nx,ny,nz))
  grid_indexG = zeros((size,3,2,2))
  grid_sizeG = zeros((size,2,2))
  grid_index = zeros((3,2,2))
  xrange = zeros(nmx)
  xrangeF = zeros(nmx)
  xrange[0] = int(nx/nmx) + overlap
  xrangeF[0] = int(nx/nmx) 
  for i in range(1,nmx-1):
    xrange[i] = int(nx/nmx) + overlap*2
    xrangeF[i] = int(nx/nmx)
  xrange[-1] = nx - (sum(xrange[0:-1]) - (nmx-2)*overlap*2 -overlap) + overlap
  xrangeF[-1] = nx - (sum(xrangeF[0:-1])) 
  yrange = zeros(nmy)
  yrangeF = zeros(nmy)
  yrange[0] = int(ny/nmy) + overlap
  yrangeF[0] = int(ny/nmy) 
  for i in range(1,nmy-1):
    yrange[i] = int(ny/nmy) + overlap*2
    yrangeF[i] = int(ny/nmy) 
  yrange[-1] = ny - (sum(yrange[0:-1]) - (nmy-2)*overlap*2 -overlap) + overlap
  yrangeF[-1] = ny - ( sum(yrangeF[0:-1]) )

  xranges = zeros(nmx)
  xranges[:] = xrange[:]
  for j in range(1,nmy):
    xrange = append(xrange,xranges[:])


  grid_indexG[0,0,0,0] = 0
  grid_indexG[0,1,0,0] = 0
  grid_indexG[0,0,1,0] = xrange[0]
  grid_indexG[0,1,1,0] = yrange[0]
  y_ind = 0.
  y_ind2 = 0.
  for j in range(0,nmy):
    x_ind = 0.
    x_ind2 = 0.
    for i in range(0,nmx):
        grid_indexG[j*nmx + i,0,0,0] = x_ind
        grid_indexG[j*nmx + i,0,1,0] = x_ind + xrange[i]
        grid_indexG[j*nmx + i,1,0,0] = y_ind
        grid_indexG[j*nmx + i,1,1,0] = y_ind + yrange[j]
        grid_indexG[j*nmx + i,0,0,1] = x_ind2
        grid_indexG[j*nmx + i,0,1,1] = x_ind2 + xrangeF[i]
        grid_indexG[j*nmx + i,1,0,1] = y_ind2
        grid_indexG[j*nmx + i,1,1,1] = y_ind2 + yrangeF[j]
        x_ind += xrange[i] - 2.*overlap
        x_ind2 += xrangeF[i]
    y_ind += yrange[j] - 2.*overlap  
    y_ind2 += yrangeF[j]

#############################
  for i in range(1,size):
        xs = grid_indexG[i,0,0,0]
        xe = grid_indexG[i,0,1,0]
        grid_sizeG[i,0,0] = xe - xs
        grid_sizeG[i,0,1] = grid_indexG[i,0,1,1] - grid_indexG[i,0,0,1]
        ys = grid_indexG[i,1,0,0]
        ye = grid_indexG[i,1,1,0]
        grid_sizeG[i,1,0] = ye - ys
        grid_sizeG[i,1,1] = grid_indexG[i,1,1,1] - grid_indexG[i,1,0,1]
        comm.send((x[xs:xe,ys:ye,:],y[xs:xe,ys:ye,:],z[xs:xe,ys:ye,:],grid_indexG[i,:,:,:],grid_sizeG[i,:,:]),dest=i, tag=13)
  grid_index[:,:,:] = grid_indexG[0,:,:,:]
  xs = grid_index[0,0,0]
  xe = grid_index[0,1,0]
  ys = grid_index[1,0,0]
  ye = grid_index[1,1,0]
  grid_size = zeros((2,2))
  grid_size[0,0] = xe - xs
  grid_size[0,1] = grid_indexG[0,0,1,1] - grid_indexG[0,0,0,1]
  grid_size[1,0] = ye - ys
  grid_size[1,1] = grid_indexG[0,1,1,1] - grid_indexG[0,1,0,1]
  gradU = zeros((nx,ny,nz,3))
  xloc,yloc,zloc = (x[xs:xe,ys:ye,:],y[xs:xe,ys:ye,:],z[xs:xe,ys:ye,:])
  MPI_grid = (xs,xe,ys,ye,grid_sizeG,grid_indexG)

else:
    xloc,yloc,zloc,grid_index,grid_size = comm.recv(source=0, tag=13)
class grid():
   A = gridMetrics(xloc,yloc,zloc)
   xperd = 0
   yperd = 0


start_x = grid_index[0,0,1] - grid_index[0,0,0]  ## get starting index of vals we want
end_x = start_x + grid_size[0,1] ## add local grid size
start_y = grid_index[1,0,1] - grid_index[1,0,0]  ## get starting index of vals we want
end_y = start_y + grid_size[1,1] ## add local grid size
MPI_grid = (start_x,end_x,start_y,end_y,grid_size)
#print('Proc' + str(MPI_RANK) + ', y range ' + str(start_y) + ' ' + str(end_y))
u = 10*sin(xloc) + 4*sin(yloc)
up_x = 10*cos(xloc)
up_y = 4*cos(yloc)
gradUL = d_dX(u,grid)


#################### CODE TO GATHER VALUES ####################################
def MPI_ALLgather(u,MPI_grid):
  if (MPI_RANK != 0):
    start_x,end_x,start_y,end_y,grid_size = MPI_grid
    dummy = zeros( ( grid_size[0,1], grid_size[1,1], 1  ) )
    dummy[:,:,:] = u[start_x:end_x,start_y:end_y,:]
    comm.Send(dummy, dest=0, tag=13)
  if (MPI_RANK == 0):
    start_x,end_x,start_y,end_y,grid_size = MPI_grid
    uG[grid_indexG[0,0,0,1]:grid_indexG[0,0,1,1],grid_indexG[0,1,0,1]:grid_indexG[0,1,1,1],:] = u[start_x:end_x,start_y:end_y,:]
    for i in range(1,size):
        dummy = zeros( ( grid_sizeG[i,0,1], grid_sizeG[i,1,1],1  ) )
        comm.Recv(dummy,source=i, tag=13)
        uG[grid_indexG[i,0,0,1]:grid_indexG[i,0,1,1],grid_indexG[i,1,0,1]:grid_indexG[i,1,1,1]] = dummy


def MPI_RANK1_send(uG):
  for i in range(1,size):
      xs = grid_indexG[i,0,0,0]
      xe = grid_indexG[i,0,1,0]
      grid_sizeG[i,0,0] = xe - xs
      grid_sizeG[i,0,1] = grid_indexG[i,0,1,1] - grid_indexG[i,0,0,1]
      ys = grid_indexG[i,1,0,0]
      ye = grid_indexG[i,1,1,0]
      grid_sizeG[i,1,0] = ye - ys
      grid_sizeG[i,1,1] = grid_indexG[i,1,1,1] - grid_indexG[i,1,0,1]
      comm.send(uG[xs:xe,ys:ye,:],dest=i, tag=13)  
  xs = grid_indexG[0,0,0,0]
  xe = grid_indexG[0,0,1,0]
  ys = grid_indexG[0,1,0,0]
  ye = grid_indexG[0,1,1,0]
  u = uG[xs:xe,ys:ye,:]
  return u

def MPI_RANKN_recv():
    u = comm.recv(source=0,tag=13) 
    return u



t2 = time.time()


## Gather onto 1 proc
MPI_ALLgather(gradUL[:,:,:,0],MPI_grid)
if (MPI_RANK == 0):
  print('Error 1= ' + str(norm(uG - 10*cos(x))))
### Distribute
if (MPI_RANK ==0):
    u = MPI_RANK1_send(uG)
else:
    u = MPI_RANKN_recv()
t3 = time.time()
MPI_ALLgather(u,MPI_grid)
gradUL = d_dX(u,grid)
MPI_ALLgather(gradUL[:,:,:,0],MPI_grid)
#                  
if (MPI_RANK == 0):
  print('Error 2= ' + str(norm(uG + 10*sin(x))))
  print(t2 - t1)
  print(t3 - t2)