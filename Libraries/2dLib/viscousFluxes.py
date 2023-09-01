from pylab import *
from gas import *
from importModule2d import d_dy,d_dx
def viscFlux(rho,rhoU,rhoV,rhoE,grid,mu,p,T):
    nx,ny = shape(rho)
    u = rhoU/rho
    v = rhoV/rho
    p[:,:] = (gam -1.)*(rhoE - 0.5*(rhoU[:,:]**2/rho[:,:] + rhoV[:,:]**2/rho[:,:]))
    T[:,:] = p[:,:]/(rho[:,:]*R)
    u_xi = d_dx(u,grid)
    u_eta = d_dy(u,grid)
    v_xi = d_dx(v,grid)
    v_eta = d_dy(v,grid)
    T_eta = d_dy(T,grid)
    T_xi = d_dx(T,grid)
    Ev = zeros((nx,ny,4))
    Fv = zeros((nx,ny,4))
    Ev[:,:,1] = mu*(4.*(grid.xi_x*u_xi + grid.eta_x*u_eta) - \
      2.*(grid.xi_y*v_xi + grid.eta_y*v_eta))/3.
    Ev[:,:,2] = mu*(grid.xi_y*u_xi + grid.eta_y*u_eta + grid.xi_x*v_xi + grid.eta_x*v_eta)
    Fv[:,:,1] = Ev[:,:,2]
    Fv[:,:,2] = mu*(-2.*(grid.xi_x*u_xi + grid.eta_x*u_eta) + \
      4.*(grid.xi_y*v_xi + grid.eta_y*v_eta))/3.
    Ev[:,:,3] = u*Ev[:,:,1] + v*Ev[:,:,2] + cp*mu/Pr*(grid.xi_x*T_xi + grid.eta_x*T_eta) 
    Fv[:,:,3] = u*Ev[:,:,2] + v*Fv[:,:,2] + cp*mu/Pr*(grid.xi_y*T_xi + grid.eta_y*T_eta)
    
    
    EvhF = zeros((nx,ny,4))
    FvhF = zeros((nx,ny,4))
    EvhF[:,:,0] = d_dx(1./grid.J*(grid.xi_x*Ev[:,:,0] + grid.xi_y*Fv[:,:,0]),grid)
    EvhF[:,:,1] = d_dx(1./grid.J*(grid.xi_x*Ev[:,:,1] + grid.xi_y*Fv[:,:,1]),grid)
    EvhF[:,:,2] = d_dx(1./grid.J*(grid.xi_x*Ev[:,:,2] + grid.xi_y*Fv[:,:,2]),grid)
    EvhF[:,:,3] = d_dx(1./grid.J*(grid.xi_x*Ev[:,:,3] + grid.xi_y*Fv[:,:,3]),grid)
    FvhF[:,:,0] = d_dy(1./grid.J*(grid.eta_x*Ev[:,:,0] + grid.eta_y*Fv[:,:,0]),grid)
    FvhF[:,:,1] = d_dy(1./grid.J*(grid.eta_x*Ev[:,:,1] + grid.eta_y*Fv[:,:,1]),grid)
    FvhF[:,:,2] = d_dy(1./grid.J*(grid.eta_x*Ev[:,:,2] + grid.eta_y*Fv[:,:,2]),grid)
    FvhF[:,:,3] = d_dy(1./grid.J*(grid.eta_x*Ev[:,:,3] + grid.eta_y*Fv[:,:,3]),grid)
    return EvhF,FvhF
  
  
