from numpy import *
from gas import *
from importModule2d import d_dx
from importModule2d import d_dy
def inviscidJacobians(rho,rhoU,rhoV,rhoE,grid):
  nx,ny = shape(rho)
  JF = zeros((nx*ny,4,4))
  JG = zeros((nx*ny,4,4))
  u = rhoU/rho
  v = rhoV/rho
  phi_sqrd = 0.5*(gam - 1.)*(u**2 + v**2)
  theta = grid.xi_x*u + grid.xi_y*v
  a1 = gam*(rhoE/rho) - phi_sqrd
  realJ = zeros((nx*ny*4,nx*ny*4))
  realG = zeros((nx*ny*4,nx*ny*4))
  JF[:,0,0] = (d_dx(grid.xi_t,grid)).flatten()
  JF[:,1,0] = (d_dx(-u*theta + grid.xi_x*phi_sqrd,grid)).flatten()
  JF[:,2,0] = (d_dx(-v*theta + grid.xi_y*phi_sqrd,grid)).flatten()
  JF[:,3,0] = (d_dx(theta*(phi_sqrd - a1),grid)).flatten()
  
  JF[:,0,1] = (d_dx(grid.xi_x,grid)).flatten()
  JF[:,1,1] = (d_dx(grid.xi_t + theta - (gam - 2.)*grid.xi_x*u,grid)).flatten()
  JF[:,2,1] = (d_dx(grid.xi_x*v - (gam - 1.)*grid.xi_y*u,grid)).flatten()
  JF[:,3,1] = (d_dx(grid.xi_x*a1 - (gam - 1.)*u*theta,grid)).flatten()

  JF[:,0,2] = (d_dx(grid.xi_y,grid)).flatten()
  JF[:,1,2] = (d_dx(grid.xi_y*u - (gam - 1.)*grid.xi_x*v,grid)).flatten()
  JF[:,2,2] = (d_dx(grid.xi_t + theta - (gam - 2.)*grid.xi_y*v,grid)).flatten()
  JF[:,3,2] = (d_dx(grid.xi_y*a1 - (gam - 1.)*v*theta,grid)).flatten()

  JF[:,0,3] = zeros(nx*ny)
  JF[:,1,3] = (d_dx((gam - 1.)*grid.xi_x,grid)).flatten()
  JF[:,2,3] = (d_dx((gam - 1.)*grid.xi_y,grid)).flatten()
  JF[:,3,3] = (d_dx(gam*theta + grid.xi_t,grid)).flatten()
  
  fill_diagonal(realJ[0::4,0::4] , JF[:,0,0])
  fill_diagonal(realJ[1::4,0::4] , JF[:,1,0])
  fill_diagonal(realJ[2::4,0::4] , JF[:,2,0])
  fill_diagonal(realJ[3::4,0::4] , JF[:,3,0])
  fill_diagonal(realJ[0::4,1::4] , JF[:,0,1])
  fill_diagonal(realJ[1::4,1::4] , JF[:,1,1])
  fill_diagonal(realJ[2::4,1::4] , JF[:,2,1])
  fill_diagonal(realJ[3::4,1::4] , JF[:,3,1])
  fill_diagonal(realJ[0::4,2::4] , JF[:,0,2])
  fill_diagonal(realJ[1::4,2::4] , JF[:,1,2])
  fill_diagonal(realJ[2::4,2::4] , JF[:,2,2])
  fill_diagonal(realJ[3::4,2::4] , JF[:,3,2])
  fill_diagonal(realJ[0::4,3::4] , JF[:,0,3])
  fill_diagonal(realJ[1::4,3::4] ,JF[:,1,3])
  fill_diagonal(realJ[2::4,3::4], JF[:,2,3])
  fill_diagonal(realJ[3::4,3::4] ,JF[:,3,3])
  
  JG[:,0,0] = (d_dy(grid.eta_t,grid)).flatten()
  JG[:,1,0] = (d_dy(-u*theta + grid.eta_x*phi_sqrd,grid)).flatten()
  JG[:,2,0] = (d_dy(-v*theta + grid.eta_y*phi_sqrd,grid)).flatten()
  JG[:,3,0] = (d_dy(theta*(phi_sqrd - a1),grid)).flatten()

  JG[:,0,1] = (d_dy(grid.eta_x,grid)).flatten()
  JG[:,1,1] = (d_dy(grid.eta_t + theta - (gam - 2.)*grid.eta_x*u,grid)).flatten()
  JG[:,2,1] = (d_dy(grid.eta_x*v - (gam - 1.)*grid.eta_y*u,grid)).flatten()
  JG[:,3,1] = (d_dy(grid.eta_x*a1 - (gam - 1.)*u*theta,grid)).flatten()

  JG[:,0,2] = (d_dy(grid.eta_y,grid)).flatten()
  JG[:,1,2] = (d_dy(grid.eta_y*u - (gam - 1.)*grid.eta_x*v,grid)).flatten()
  JG[:,2,2] = (d_dy(grid.eta_t + theta - (gam - 2.)*grid.eta_y*v,grid)).flatten()
  JG[:,3,2] = (d_dy(grid.eta_y*a1 - (gam - 1.)*v*theta,grid)).flatten()

  JG[:,0,3] = 0.
  JG[:,1,3] = (d_dy((gam - 1.)*grid.eta_x,grid)).flatten()
  JG[:,2,3] = (d_dy((gam - 1.)*grid.eta_y,grid)).flatten()
  JG[:,3,3] = (d_dy(gam*theta + grid.eta_t,grid)).flatten()
  
  
  fill_diagonal(realG[0::4,0::4],JG[:,0,0])
  fill_diagonal(realG[1::4,0::4],JG[:,1,0])
  fill_diagonal(realG[2::4,0::4],JG[:,2,0])
  fill_diagonal(realG[3::4,0::4],JG[:,3,0])
  fill_diagonal(realG[0::4,1::4],JG[:,0,1])
  fill_diagonal(realG[1::4,1::4],JG[:,1,1])
  fill_diagonal(realG[2::4,1::4],JG[:,2,1])
  fill_diagonal(realG[3::4,1::4],JG[:,3,1])
  fill_diagonal(realG[0::4,2::4],JG[:,0,2])
  fill_diagonal(realG[1::4,2::4],JG[:,1,2])
  fill_diagonal(realG[2::4,2::4],JG[:,2,2])
  fill_diagonal(realG[3::4,2::4],JG[:,3,2])
  fill_diagonal(realG[0::4,3::4],JG[:,0,3])
  fill_diagonal(realG[1::4,3::4],JG[:,1,3])
  fill_diagonal(realG[2::4,3::4],JG[:,2,3])
  fill_diagonal(realG[3::4,3::4],JG[:,3,3])
  return -realJ,-realG