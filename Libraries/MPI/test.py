from pylab import *
from mpi4py import MPI
from getNeighbours import getNeighbours, refreshNeighbourVals
comm = MPI.COMM_WORLD
MPI_RANK = comm.Get_rank()
size = comm.Get_size()
nmx = 2
nmy = 2
u = zeros((10,10))
neighbour = getNeighbours(comm,MPI_RANK,size,nmx,nmy)
refreshNeighbourVals(u,neighbour)

if (MPI_RANK == 1):
    pcolormesh(u.transpose())
    colorbar()
    show()