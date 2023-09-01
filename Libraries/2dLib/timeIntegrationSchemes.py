from pylab import *
import scipy.sparse.linalg
from importModule2d import computeFlux
from BoundaryConditions import applyBCs
from mesh import normals
from gas import *
from Jacobians import inviscidJacobians
import time
rk4const = array([1./4,1./3,1./2,1.])
##RK4 integration variables

def rk4(dt,rho,rhoU,rhoV,rhoE,grid,p,a):
    nx,ny = shape(rho)
    rho1 = zeros((nx,ny))
    rhoU1 = zeros((nx,ny))
    rhoV1 = zeros((nx,ny))
    rhoE1 = zeros((nx,ny))
    rho1[:,:] = rho[:,:]
    rhoU1[:,:] = rhoU[:,:]
    rhoV1[:,:] = rhoV[:,:]
    rhoE1[:,:] = rhoE[:,:]
    for i in range(0,4):
        Fflux,Gflux = computeFlux(rho,rhoU,rhoV,rhoE,grid,p,a)
        ## continuity equation
        rhoResid = Fflux[:,:,0] + Gflux[:,:,0] 
        rho[:,:] = rho1[:,:] - dt*rk4const[i]*grid.J*( rhoResid )
        
        ## U Momentum Equation
        rhoUResid = Fflux[:,:,1] + Gflux[:,:,1] 
        rhoU[:,:] = rhoU1[:,:] - dt*rk4const[i]*grid.J*( rhoUResid )
        
        ## V Momentum Equation
        rhoVResid = Fflux[:,:,2] + Gflux[:,:,2] 
        rhoV[:,:] = rhoV1[:,:] - dt*rk4const[i]*grid.J*( rhoVResid )
        
        ## Energy Equation
        rhoEResid = Fflux[:,:,3] + Gflux[:,:,3] 
        rhoE[:,:] = rhoE1[:,:] - dt*rk4const[i]*grid.J*( rhoEResid )
        rho,rhoU,rhoV,rhoE = applyBCs(rho,rhoU,rhoV,rhoE,normals)
        p[:,:] = (gam -1.)*(rhoE - 0.5*(rhoU[:,:]**2/rho[:,:] + rhoV[:,:]**2/rho[:,:]))
        a[:,:] = (gam*p[:,:]/rho[:,:])**0.5
    return rho,rhoU,rhoV,rhoE,rhoUResid
    
def implicitEuler(dt,rho,rhoU,rhoV,rhoE,grid,p,a):
    nx,ny = shape(rho)
    JF,JG = inviscidJacobians(rho,rhoU,rhoV,rhoE,grid)
    Fflux,Gflux = computeFlux(rho,rhoU,rhoV,rhoE,grid,p,a)
    ## continuity equation
    rhoResid = Fflux[:,:,0] + Gflux[:,:,0] 
    ## U Momentum Equation
    rhoUResid = Fflux[:,:,1] + Gflux[:,:,1] 
    ## V Momentum Equation
    rhoVResid = Fflux[:,:,2] + Gflux[:,:,2] 
    ## Energy Equation
    rhoEResid = Fflux[:,:,3] + Gflux[:,:,3]
    RHS = zeros((nx*ny*4))
    RHS[0::4]= (-Fflux[:,:,0] - Gflux[:,:,0]).flatten()
    RHS[1::4]= (-Fflux[:,:,1] - Gflux[:,:,1]).flatten()
    RHS[2::4]= (-Fflux[:,:,2] - Gflux[:,:,2]).flatten()
    RHS[3::4]= (-Fflux[:,:,3] - Gflux[:,:,3]).flatten()
    Am = eye(nx*ny*4)/dt - JF - JG
    t1 = time.time()
    dU = scipy.sparse.linalg.spsolve(Am,RHS)
    t2 = time.time()
    print(t2 - t1)
    rho = rho + grid.J*reshape(dU[0::4],(nx,ny))
    rhoU = rhoU + grid.J*reshape(dU[1::4],(nx,ny))
    rhoV = rhoV + grid.J*reshape(dU[2::4],(nx,ny))
    rhoE = rhoE + grid.J*reshape(dU[3::4],(nx,ny))
    rho,rhoU,rhoV,rhoE = applyBCs(rho,rhoU,rhoV,rhoE,normals)
    p[:,:] = (gam -1.)*(rhoE - 0.5*(rhoU[:,:]**2/rho[:,:] + rhoV[:,:]**2/rho[:,:]))
    a[:,:] = (gam*p[:,:]/rho[:,:])**0.5
    return rho,rhoU,rhoV,rhoE
