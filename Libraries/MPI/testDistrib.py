from pylab import *
import sys
sys.path.append("../")
import time
from libraries import *
from MPIlibraries import *
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
nmx = 2
nmy = 2
nmz = 1
size = nmx*nmy
if (MPI_RANK == 0):
  from mesh import x,y,z
  nxg,nyg,nzg = shape(x)
  rhoG = cos(x)
  rhoUG = cos(x)*sin(y)
  rhoVG = cos(x)*sin(y)*2
  rhoWG = cos(x)*sin(y)*3
  rhoEG = cos(x)*sin(y)*5




if (MPI_RANK ==0):
  xloc,yloc,zloc,grid_index,grid_size,grid_indexG,grid_sizeG = MPI_send_mesh(nmx,nmy,nmz)
else:
  xloc,yloc,zloc,grid_index,grid_size = MPI_recv_mesh(comm)

nxl,nyl,nzl = shape(xloc)
rho = zeros((nxl,nyl,nzl))
rhoU = zeros((nxl,nyl,nzl))
rhoV = zeros((nxl,nyl,nzl))
rhoW = zeros((nxl,nyl,nzl))
rhoE = zeros((nxl,nyl,nzl))
if (MPI_RANK == 0):
    rho[:,:,:],rhoU[:,:,:],rhoV[:,:,:],rhoW[:,:,:],rhoE[:,:,:] = MPI_sendQ(rhoG,rhoUG,rhoVG,rhoWG,rhoEG,grid_indexG,grid_sizeG)

else:
  rho[:,:,:],rhoU[:,:,:],rhoV[:,:,:],rhoW[:,:,:],rhoE[:,:,:] = MPI_recvQ()

rho[:,:,:] = rho[:,:,:]**2*cos(xloc)
if (MPI_RANK == 0):
  rhoG2,rhoUG2,rhoVG2,rhoWG2,rhoEG2 = MPI_gatherQ(rho,rhoU,rhoV,rhoW,rhoE,nxg,nyg,nzg,grid_index,grid_size,grid_indexG,grid_sizeG)
else:
  MPI_reduceQ(rho,rhoU,rhoV,rhoW,rhoE,grid_index,grid_size)


if (MPI_RANK == 0):
  figure(1)
  print(norm(rhoG[:,:,0]**2*cos(x[:,:,0]) - rhoG2[:,:,0]))