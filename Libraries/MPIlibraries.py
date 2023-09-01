from pylab import *
from mpi4py import MPI
##############################################
################# CODE TO BREAK UP A MESH AND 
################# DISTRIBUTE TO PROCS 
def initiateMPI():
  comm = MPI.COMM_WORLD
  MPI_RANK = comm.Get_rank()
  size = comm.Get_size()
  return comm,MPI_RANK,size


def MPI_send_mesh(nmx,nmy,nmz):
  comm,MPI_RANK,size = initiateMPI()
  #################
  from mesh import x,y,z
  overlap = 2
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
  return xloc,yloc,zloc,grid_index,grid_size,grid_indexG,grid_sizeG


def MPI_recv_mesh(comm):
  xloc,yloc,zloc,grid_index,grid_size = comm.recv(source=0, tag=13)
  return xloc,yloc,zloc,grid_index,grid_size


### Code to distribute solution vector from proc one to all other procs
def MPI_sendQ(rhoG,rhoUG,rhoVG,rhoWG,rhoEG,grid_indexG,grid_sizeG):
  comm,MPI_RANK,size = initiateMPI()
  for i in range(1,size):
      xs = grid_indexG[i,0,0,0]
      xe = grid_indexG[i,0,1,0]
      grid_sizeG[i,0,0] = xe - xs
      grid_sizeG[i,0,1] = grid_indexG[i,0,1,1] - grid_indexG[i,0,0,1]
      ys = grid_indexG[i,1,0,0]
      ye = grid_indexG[i,1,1,0]
      grid_sizeG[i,1,0] = ye - ys
      grid_sizeG[i,1,1] = grid_indexG[i,1,1,1] - grid_indexG[i,1,0,1]
      comm.send((rhoG[xs:xe,ys:ye,:],rhoUG[xs:xe,ys:ye,:],rhoVG[xs:xe,ys:ye,:],rhoWG[xs:xe,ys:ye,:],rhoEG[xs:xe,ys:ye,:]),dest=i, tag=13)  
  xs = grid_indexG[0,0,0,0]
  xe = grid_indexG[0,0,1,0]
  ys = grid_indexG[0,1,0,0]
  ye = grid_indexG[0,1,1,0]
  rho = rhoG[xs:xe,ys:ye,:]
  rhoU = rhoUG[xs:xe,ys:ye,:]
  rhoV = rhoVG[xs:xe,ys:ye,:]
  rhoW = rhoWG[xs:xe,ys:ye,:]
  rhoE = rhoEG[xs:xe,ys:ye,:]
  return rho,rhoU,rhoV,rhoW,rhoE
### Code for other procs to recieve solution vector
def MPI_recvQ():
    comm,MPI_RANK,size = initiateMPI()
    rho,rhoU,rhoV,rhoW,rhoE = comm.recv(source=0,tag=13) 
    return rho,rhoU,rhoV,rhoW,rhoE


### Code for sending variables from all procs to one proc
def MPI_reduceQ(rho,rhoU,rhoV,rhoW,rhoE,grid_index,grid_size):
    comm,MPI_RANK,size = initiateMPI()
    start_x = grid_index[0,0,1] - grid_index[0,0,0]  ## get starting index of vals we want
    end_x = start_x + grid_size[0,1] ## add local grid size
    start_y = grid_index[1,0,1] - grid_index[1,0,0]  ## get starting index of vals we want
    end_y = start_y + grid_size[1,1] ## add local grid size
    comm.send((rho[start_x:end_x,start_y:end_y,:],rhoU[start_x:end_x,start_y:end_y,:],rhoV[start_x:end_x,start_y:end_y,:],rhoW[start_x:end_x,start_y:end_y,:],rhoE[start_x:end_x,start_y:end_y,:]), dest=0, tag=13)

def MPI_gatherQ(rho,rhoU,rhoV,rhoW,rhoE,rhoG,rhoUG,rhoVG,rhoWG,rhoEG,grid_index,grid_size,grid_indexG,grid_sizeG):
    comm,MPI_RANK,size = initiateMPI()
    start_x = grid_index[0,0,1] - grid_index[0,0,0]  ## get starting index of vals we want
    end_x = start_x + grid_size[0,1] ## add local grid size
    start_y = grid_index[1,0,1] - grid_index[1,0,0]  ## get starting index of vals we want
    end_y = start_y + grid_size[1,1] ## add local grid size
    rhoG[grid_indexG[0,0,0,1]:grid_indexG[0,0,1,1],grid_indexG[0,1,0,1]:grid_indexG[0,1,1,1],:] = rho[start_x:end_x,start_y:end_y,:]
    rhoUG[grid_indexG[0,0,0,1]:grid_indexG[0,0,1,1],grid_indexG[0,1,0,1]:grid_indexG[0,1,1,1],:] = rhoU[start_x:end_x,start_y:end_y,:]
    rhoVG[grid_indexG[0,0,0,1]:grid_indexG[0,0,1,1],grid_indexG[0,1,0,1]:grid_indexG[0,1,1,1],:] = rhoV[start_x:end_x,start_y:end_y,:]
    rhoWG[grid_indexG[0,0,0,1]:grid_indexG[0,0,1,1],grid_indexG[0,1,0,1]:grid_indexG[0,1,1,1],:] = rhoW[start_x:end_x,start_y:end_y,:]
    rhoEG[grid_indexG[0,0,0,1]:grid_indexG[0,0,1,1],grid_indexG[0,1,0,1]:grid_indexG[0,1,1,1],:] = rhoE[start_x:end_x,start_y:end_y,:]
    for i in range(1,size):
        dummy = zeros( ( grid_sizeG[i,0,1], grid_sizeG[i,1,1],1  ) )
        rhoG[grid_indexG[i,0,0,1]:grid_indexG[i,0,1,1],grid_indexG[i,1,0,1]:grid_indexG[i,1,1,1]],\
        rhoUG[grid_indexG[i,0,0,1]:grid_indexG[i,0,1,1],grid_indexG[i,1,0,1]:grid_indexG[i,1,1,1]],\
        rhoVG[grid_indexG[i,0,0,1]:grid_indexG[i,0,1,1],grid_indexG[i,1,0,1]:grid_indexG[i,1,1,1]],\
        rhoWG[grid_indexG[i,0,0,1]:grid_indexG[i,0,1,1],grid_indexG[i,1,0,1]:grid_indexG[i,1,1,1]],\
        rhoEG[grid_indexG[i,0,0,1]:grid_indexG[i,0,1,1],grid_indexG[i,1,0,1]:grid_indexG[i,1,1,1]] = comm.recv(source=i, tag=13)
    return rhoG,rhoUG,rhoVG,rhoWG,rhoEG


