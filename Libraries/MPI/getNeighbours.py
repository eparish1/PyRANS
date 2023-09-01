from pylab import *
from mpi4py import MPI
comm = MPI.COMM_WORLD
MPI_RANK = comm.Get_rank()
size = comm.Get_size()

### Code to get the neighbouring processors for a square/rectangular decomposition
### 10/20

def getNeighbours(comm,MPI_RANK,size,nmx,nmy):

  if (size != nmx*nmy):
      if (MPI_RANK == 0):
        print('Error, number of decompositions does not equal the number of processors')
      quit()

  if (MPI_RANK == 0):
    neighbour_x_L = zeros(nmx*nmy)
    neighbour_x_R = zeros(nmx*nmy)
    neighbour_y_u = zeros(nmx*nmy)
    neighbour_y_d = zeros(nmx*nmy)

    counter = 0
    for j in range(0,nmy):
        for i in range(0,nmx-1):
            neighbour_x_R[counter] = counter + 1 
            counter = counter + 1
        neighbour_x_R[counter ] = -4
        counter = counter + 1
    counter = 0
    for j in range(0,nmy):
        neighbour_x_L[counter ] = -3
        counter = counter + 1
        for i in range(1,nmx):
            neighbour_x_L[counter] = counter -1 
            counter = counter + 1

    counter = 0
    for i in range(0,nmx):
        counter = i
        for j in range(0,nmy-1):
            neighbour_y_u[counter] = counter + nmx
            counter = counter + nmx
        neighbour_y_u[counter] = -2
        #counter = counter -nmx + 1
    
    counter = 0
    for i in range(0,nmx):
        counter = i
        neighbour_y_d[counter] = -1
        counter = counter + nmx
        for j in range(1,nmy):
            neighbour_y_d[counter] = counter - nmx
            counter = counter + nmx
        
  if (MPI_RANK == 0):
      neighbour = zeros(4,dtype = int)
      for i in range(1,size):
          comm.send((neighbour_x_L[i],neighbour_x_R[i],neighbour_y_d[i],neighbour_y_u[i]),dest=i,tag=13)
      neighbour[0] = neighbour_x_L[0]
      neighbour[1] = neighbour_x_R[0]
      neighbour[2] = neighbour_y_d[0]
      neighbour[3] = neighbour_y_u[0]
      del neighbour_x_L,neighbour_x_R,neighbour_y_d,neighbour_y_u,counter
  else:
      neighbour = zeros(4,dtype = int)
      neighbour[0],neighbour[1],neighbour[2],neighbour[3] = comm.recv(source=0, tag=13)
      #counter = counter -nmx + 1
  return neighbour 


def refreshNeighbourVals(u,neighbour):
    u[:,:] = MPI_RANK
    overlap = 1
    if (neighbour[0] >= 0):
      comm.send(u[2,:],dest=neighbour[0],tag=14)
      #print(neighbour[0])
      
    if (neighbour[1] >= 0):
      comm.send(u[-3,:],dest=neighbour[1],tag=13)
    #  
    if (neighbour[2] >= 0):
      comm.send(u[:,2],dest=neighbour[2],tag=15)
    #  
    if (neighbour[3] >= 0):
      comm.send(u[:,-3],dest=neighbour[3],tag=16)
    
    #________________________
    if (neighbour[0] >= 0):
      u[0,:] = comm.recv(source=neighbour[0],tag=13)
      
    if (neighbour[1] >= 0):
      u[-1,:] = comm.recv(source=neighbour[1],tag=14)
    #  
    if (neighbour[2] >= 0):
      u[:,0] = comm.recv(source=neighbour[2],tag=16)
    #  
    if (neighbour[3] >= 0):
      u[:,-1] = comm.recv(source=neighbour[3],tag=15)
      
      
def MPI_recv_mesh(comm):
  xloc,yloc,grid_index,grid_size = comm.recv(source=0, tag=13)
  return xloc,yloc,grid_index,grid_size
        
def MPI_send_mesh(nmx,nmy):
  from mesh import x,y
  overlap = 1
  nx,ny = shape(x)
  uG = zeros((nx,ny))
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
  xloc,yloc= (x[xs:xe,ys:ye],y[xs:xe,ys:ye],z[xs:xe,ys:ye])
  return xloc,yloc,grid_index,grid_size,grid_indexG,grid_sizeG