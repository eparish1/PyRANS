from pylab import *
from evtk.hl import gridToVTK 
import time
sys.path.append("../../Libraries/")
from libraries import * ## add our libraries for the code (differencing schemes, etc)
## Create grid
xperd = 0
yperd = 1
nr = 177
ntheta = 132
dtheta = 2.*pi/(ntheta)
Rg = 40.**(1./nr)
r = zeros(nr)
r[0] = 1.
for i in range(0,nr-1):
  r[i+1] = r[i] + r[i]*(Rg - 1)
#r = linspace(0.1,3,nr)
nz = 1
theta = linspace(dtheta/2.,2*pi-dtheta/2.,ntheta)
x = zeros((nr,ntheta,nz))
y = zeros((nr,ntheta,nz))
z = zeros((nr,ntheta,nz))
for i in range(0,nr):
    for j in range(0,ntheta):
        for k in range(0,nz):
            x[i,j,k] = r[i]*cos(theta[j])
            y[i,j,k] = r[i]*sin(theta[j])
            z[i,j,k] = k


#x = linspace(0,1,100)
#y = linspace(0,1,100)
#z = linspace(0,1,1)
#y,x,z = meshgrid(y,x,z)
## For the circle mesh we also need the normals at the surface of the cylinder
#  for velocity extrapolation
nx,ny,nz = shape(x)
normals = zeros((nx,ny,nz,3))
normMag = zeros((nx,ny,nz))

# compute magnitude
#normMag = ( (y[0:-1,1::,0] - y[0:-1,0:-1,0])**2 +  (x[0:-1,1::,0] - x[0:-1,0:-1,0])**2)**0.5
# now get individual normalized x and y mags
# x normal = y_ip1 - y_i
normals[:,1:-1,0,0] =0.5*( y[:,2::,0] - y[:,0:-2,0]) #/ normMag
normals[:,0,0,0] = 0.5*( y[:,1,0] - y[:,-2,0])
normals[:,-1,0,0] = normals[:,0,0,0]
# y normals = -(x_ip1 - x_i)
normals[:,1:-1,0,1] =-0.5*( x[:,2::,0] - x[:,0:-2,0])
normals[:,0,0,1] = -0.5*( x[:,1,0] - x[:,-2,0])
normals[:,-1,0,1] = normals[:,0,0,1]


normMag[:,:,:] = (abs(normals[:,:,:,0])**2  + abs(normals[:,:,:,1])**2)**0.5
normals[:,:,:,0] = normals[:,:,:,0]/normMag[:,:,:]
normals[:,:,:,1] = normals[:,:,:,1]/normMag[:,:,:]
# no z normals, 2D grid



