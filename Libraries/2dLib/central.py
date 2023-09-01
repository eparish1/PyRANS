from pylab import *
#from libraries2d import d_dx
#from libraries2d import d_dy 
from importModule2d import d_dx
from importModule2d import d_dy
def centralFlux(rho,rhoU,rhoV,rhoE,grid,p,a):
        ## compute contrivariant velocities
        nx,ny = shape(rho)
        u = rhoU/rho
        v = rhoV/rho
        uC = grid.xi_x*u + grid.xi_y*v
        vC = grid.eta_x*u + grid.eta_y*v
        ## Compute fluxes
        F = zeros((nx,ny,4))
        G = zeros((nx,ny,4))
        F[:,:,0] = (rho*uC)/grid.J
        F[:,:,1] = (rho*u*uC + grid.xi_x*p)/grid.J
        F[:,:,2] = (rho*v*uC + grid.xi_y*p)/grid.J
        F[:,:,3] = (rhoE*uC + p*uC)/grid.J
        
        G[:,:,0] = (rho*vC)/grid.J
        G[:,:,1] = (rho*u*vC + grid.eta_x*p)/grid.J
        G[:,:,2] = (rho*v*vC + grid.eta_y*p)/grid.J
        G[:,:,3] = (rhoE*vC + p*vC)/grid.J
        
        Fflux = zeros((nx,ny,4))
        Gflux = zeros((nx,ny,4))
        Fflux[:,:,0] = d_dx(F[:,:,0],grid)
        Fflux[:,:,1] = d_dx(F[:,:,1],grid)
        Fflux[:,:,2] = d_dx(F[:,:,2],grid)
        Fflux[:,:,3] = d_dx(F[:,:,3],grid)
        
        Gflux[:,:,0] = d_dy(G[:,:,0],grid)
        Gflux[:,:,1] = d_dy(G[:,:,1],grid)
        Gflux[:,:,2] = d_dy(G[:,:,2],grid)
        Gflux[:,:,3] = d_dy(G[:,:,3],grid)
        
        return Fflux,Gflux
