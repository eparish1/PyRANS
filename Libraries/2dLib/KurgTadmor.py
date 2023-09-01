from pylab import *

def buildFlux(rho,rhoU,rhoV,rhoE):
    nx,ny = shape(rho)
    Fflux = zeros((nx,ny,4))
    Gflux = zeros((nx,ny,4))
    Fflux[:,:,0] = rho*uC
def kurgTadmorFlux(rho,rhoU,rhoV,rhoE,grid,p,a):
    Fm,Fp,Gm,Gp = kurgTadmorSplitting(rho,rhoU,rhoV,rhoE,grid,p,a)
    nx,ny = shape(rho)
    Fflux = zeros((nx,ny,4))
    Gflux = zeros((nx,ny,4))
    Fflux[:,:,0] = (Fp[:,:,0] - Fm[:,:,0])
    Fflux[:,:,1] = (Fp[:,:,1] - Fm[:,:,1])
    Fflux[:,:,2] = (Fp[:,:,2] - Fm[:,:,2])
    Fflux[:,:,3] = (Fp[:,:,3] - Fm[:,:,3])
    
    Gflux[:,:,0] = (Gp[:,:,0] - Gm[:,:,0])
    Gflux[:,:,1] = (Gp[:,:,1] - Gm[:,:,1])
    Gflux[:,:,2] = (Gp[:,:,2] - Gm[:,:,2])
    Gflux[:,:,3] = (Gp[:,:,3] - Gm[:,:,3])
    return Fflux,Gflux
def kurgTadmorSplitting(rho,rhoU,rhoV,rhoE,grid,p,a):
    u = rhoU/rho
    v = rhoV/rho
    J = grid.J
    uC = grid.xi_x*u + grid.xi_y*v
    vC = grid.eta_x*u + grid.eta_y*v

    rhoLm_x,rhoLp_x,rhoRm_x,rhoRp_x = linearReconstruct_x2(rho/J,grid) 
    rhoLm_y,rhoLp_y,rhoRm_y,rhoRp_y = linearReconstruct_y2(rho/J,grid)
    
    rhoULm_x,rhoULp_x,rhoURm_x,rhoURp_x = linearReconstruct_x2(rhoU/J,grid)
    rhoULm_y,rhoULp_y,rhoURm_y,rhoURp_y = linearReconstruct_y2(rhoU/J,grid)
    
    rhoVLm_x,rhoVLp_x,rhoVRm_x,rhoVRp_x = linearReconstruct_x2(rhoV/J,grid)       
    rhoVLm_y,rhoVLp_y,rhoVRm_y,rhoVRp_y = linearReconstruct_y2(rhoV/J,grid)
    
    uCLm_x,uCLp_x,uCRm_x,uCRp_x = linearReconstruct_x2(uC,grid)
    uCLm_y,uCLp_y,uCRm_y,uCRp_y = linearReconstruct_y2(uC,grid)

    vCLm_x,vCLp_x,vCRm_x,vCRp_x = linearReconstruct_x2(vC,grid)
    vCLm_y,vCLp_y,vCRm_y,vCRp_y = linearReconstruct_y2(vC,grid)


    rhoELm_x,rhoELp_x,rhoERm_x,rhoERp_x = linearReconstruct_x2(rhoE/J,grid)
    rhoELm_y,rhoELp_y,rhoERm_y,rhoERp_y = linearReconstruct_y2(rhoE/J,grid)
    pLm_x,pLp_x,pRm_x,pRp_x = linearReconstruct_x2(p/J,grid)
    pLm_y,pLp_y,pRm_y,pRp_y = linearReconstruct_y2(p/J,grid)

    nx,ny = shape(u)
    F = zeros((nx,ny,4))
    Fm = zeros((nx,ny,4))
    Fp = zeros((nx,ny,4))
    G = zeros((nx,ny,4))
    Gm = zeros((nx,ny,4))
    Gp = zeros((nx,ny,4))
    lam1Rp = (-a*grid.XI+rhoVRp_x/rhoRp_x*grid.xi_y + rhoURp_x/rhoRp_x*grid.xi_x+grid.xi_t)
    lam2Rp = ( a*grid.XI+rhoVRp_x/rhoRp_x*grid.xi_y + rhoURp_x/rhoRp_x*grid.xi_x+grid.xi_t)
    lam3Rp = ( rhoVRp_x/rhoRp_x*grid.xi_y+rhoURp_x/rhoRp_x*grid.xi_x+grid.xi_t)
    
    lam1Lp = (-a*grid.XI+rhoVLp_x/rhoLp_x*grid.xi_y + rhoULp_x/rhoLp_x*grid.xi_x+grid.xi_t)
    lam2Lp = ( a*grid.XI+rhoVLp_x/rhoLp_x*grid.xi_y + rhoULp_x/rhoLp_x*grid.xi_x+grid.xi_t)
    lam3Lp = ( rhoVLp_x/rhoLp_x*grid.xi_y+rhoULp_x/rhoLp_x*grid.xi_x+grid.xi_t)

    lam1Rm = (-a*grid.XI+rhoVRm_x/rhoRm_x*grid.xi_y + rhoURm_x/rhoRm_x*grid.xi_x+grid.xi_t)
    lam2Rm = ( a*grid.XI+rhoVRm_x/rhoRm_x*grid.xi_y + rhoURm_x/rhoRm_x*grid.xi_x+grid.xi_t)
    lam3Rm = ( rhoVRm_x/rhoRm_x*grid.xi_y+rhoURm_x/rhoRm_x*grid.xi_x+grid.xi_t)
    
    lam1Lm = (-a*grid.XI+rhoVLm_x/rhoLm_x*grid.xi_y + rhoULm_x/rhoLm_x*grid.xi_x+grid.xi_t)
    lam2Lm = ( a*grid.XI+rhoVLm_x/rhoLm_x*grid.xi_y + rhoULm_x/rhoLm_x*grid.xi_x+grid.xi_t)
    lam3Lm = ( rhoVLm_x/rhoLm_x*grid.xi_y+rhoULm_x/rhoLm_x*grid.xi_x+grid.xi_t)
    
    amaxL_x = zeros(shape(u))
    amaxL_x[:,:] = fmax(lam1Lm,lam1Rm)

    amaxR_x = zeros(shape(u))
    amaxR_x[:,:] = fmax(lam1Lp,lam1Rp)
    
    xi_x = grid.xi_x
    xi_y = grid.xi_y
    
    Fm[:,:,0] = 0.5*( rhoRm_x*uCRm_x + rhoLm_x*uCLm_x ) - amaxL_x*(rhoRm_x - rhoLm_x)
    Fm[:,:,1] = 0.5*( (rhoURm_x*uCRm_x + xi_x*pRm_x) + (rhoULm_x*uCLm_x + xi_x*pLm_x) ) - amaxL_x*(rhoURm_x - rhoULm_x)
    Fm[:,:,2] = 0.5*((rhoVRm_x*uCRm_x + xi_y*pRm_x) + (rhoVLm_x*uCLm_x + xi_y*pLm_x) ) - amaxL_x*(rhoVRm_x - rhoVLm_x) 
    Fm[:,:,3] = 0.5*((rhoERm_x*uCRm_x + pRm_x*uCRm_x) + (rhoELm_x*uCLm_x + pLm_x*uCLm_x)) - amaxL_x*(rhoERm_x - rhoELm_x)
    
    Fp[:,:,0] = 0.5*( (rhoRp_x*uCRp_x) + (rhoLp_x*uCLp_x) ) - amaxR_x*(rhoRp_x - rhoLp_x)
    Fp[:,:,1] = 0.5*((rhoURp_x*uCRp_x + xi_x*pRp_x) + (rhoULp_x*uCLp_x + xi_x*pLp_x) ) - amaxR_x*(rhoURp_x - rhoULp_x)
    Fp[:,:,2] = 0.5*((rhoVRp_x*uCRp_x + xi_y*pRp_x) + (rhoVLp_x*uCLp_x + xi_y*pLp_x) ) - amaxR_x*(rhoVRp_x - rhoVLp_x)
    Fp[:,:,3] = 0.5*((rhoERp_x*uCRp_x + pRp_x*uCRp_x) + (rhoELp_x*uCLp_x + pLp_x*uCLp_x)) - amaxR_x*(rhoERp_x - rhoELp_x)

    eta_x = grid.eta_x
    eta_y = grid.eta_y
        
    lam1Rp = (-a*grid.ETA+rhoVRp_y/rhoRp_y*grid.eta_y+rhoURp_y/rhoRp_y*grid.eta_x+grid.eta_t)
    lam2Rp = ( a*grid.ETA+rhoVRp_y/rhoRp_y*grid.eta_y+rhoURp_y/rhoRp_y*grid.eta_x+grid.eta_t)
    lam3Rp = ( rhoVRp_y/rhoRp_y*grid.eta_y+rhoURp_y/rhoRp_y*grid.eta_x+grid.eta_t  )
 
    lam1Rm = (-a*grid.ETA+rhoVRm_y/rhoRm_y*grid.eta_y+rhoURm_y/rhoRm_y*grid.eta_x+grid.eta_t)
    lam2Rm = ( a*grid.ETA+rhoVRm_y/rhoRm_y*grid.eta_y+rhoURm_y/rhoRm_y*grid.eta_x+grid.eta_t)
    lam3Rm = ( rhoVRm_y/rhoRm_y*grid.eta_y+rhoURm_y/rhoRm_y*grid.eta_x+grid.eta_t  )
    
    lam1Lp = (-a*grid.ETA+rhoVLp_y/rhoLp_y*grid.eta_y+rhoULp_y/rhoLp_y*grid.eta_x+grid.eta_t)
    lam2Lp = ( a*grid.ETA+rhoVLp_y/rhoLp_y*grid.eta_y+rhoULp_y/rhoLp_y*grid.eta_x+grid.eta_t)
    lam3Lp = ( rhoVLp_y/rhoLp_y*grid.eta_y+rhoULp_y/rhoLp_y*grid.eta_x+grid.eta_t  )
    
    lam1Lm = (-a*grid.ETA+rhoVRp_y/rhoRp_y*grid.eta_y+rhoURp_y/rhoRp_y*grid.eta_x+grid.eta_t)
    lam2Lm = ( a*grid.ETA+rhoVRp_y/rhoRp_y*grid.eta_y+rhoURp_y/rhoRp_y*grid.eta_x+grid.eta_t)
    lam3Lm = ( rhoVRp_y/rhoRp_y*grid.eta_y+rhoURp_y/rhoRp_y*grid.eta_x+grid.eta_t  )

    amaxL_y = zeros(shape(u))
    amaxL_y[:,:] = fmax(lam1Lm,lam1Rm)

    amaxR_y = zeros(shape(u))
    amaxR_y[:,:] = fmax(lam1Lp,lam1Rp)
        
    
    Gm[:,:,0] = 0.5*((rhoRm_y*vCRm_y) + (rhoLm_y*vCLm_y) ) - amaxL_y*(rhoRm_y - rhoLm_y)
    Gm[:,:,1] = 0.5*((rhoURm_y*vCRm_y + eta_x*pRm_y) + (rhoULm_y*vCLm_y + eta_x*pLm_y) ) - amaxL_y*(rhoURm_y - rhoULm_y)
    Gm[:,:,2] = 0.5*((rhoVRm_y*vCRm_y + eta_y*pRm_y) + (rhoVLm_y*vCLm_y + eta_y*pLm_y) ) - amaxL_y*(rhoVRm_y - rhoVLm_y)
    Gm[:,:,3] = 0.5*((rhoERm_y*vCRm_y + pRm_y*vCRm_y) + (rhoELm_y*vCLm_y + pLm_y*vCLm_y)) - amaxL_y*(rhoERm_y - rhoELm_y)
    
    Gp[:,:,0] = 0.5*((rhoRp_y*vCRp_y) + (rhoLp_y*vCLp_y) ) - amaxR_y*(rhoRp_y - rhoLp_y)
    Gp[:,:,1] = 0.5*((rhoURp_y*vCRp_y + eta_x*pRp_y) + (rhoULp_y*vCLp_y + eta_x*pLp_y) ) - amaxR_y*(rhoURp_y - rhoULp_y)
    Gp[:,:,2] = 0.5*((rhoVRp_y*vCRp_y + eta_y*pRp_y) + (rhoVLp_y*vCLp_y + eta_y*pLp_y) ) - amaxR_y*(rhoVRp_y - rhoVLp_y)
    Gp[:,:,3] = 0.5*((rhoERp_y*vCRp_y + pRp_y*vCRp_y) + (rhoELp_y*vCLp_y + pLp_y*vCLp_y)) - amaxR_y*(rhoERp_y - rhoELp_y)
    
    return Fm,Fp,Gm,Gp
def fluxLimiter(r):
  phi = zeros(shape(r))
  phi[:] = 1.5*(r[:]**2 + r[:])/(r[:]**2 + r[:] + 1. + 1e-20)
  #phi[:] = fmax(zeros(nx),fmax(fmin(2*r,ones(nx)),fmin(r,ones(nx)*2)))
  #phi[:] = (2.*r)/(r**2 + 1.)
  #phi[:,:] = (r**2 + r)/(r**2 + 1)
  #phi = fmax(0,fmin(fmin(2.*r,(2. + r)/3.),2.))
  return phi

def linearReconstruct_x(u,grid):
    rx = zeros(shape(u))
    rx[1:-1,:] = (u[1:-1,:] - u[0:-2,:])/(u[2::,:] - u[1:-1,:] + 1.e-30)
    rx[0,:] = rx[1,:]
    rx[-1,:] = rx[-2,:]
    phi = fluxLimiter(rx)
    uLp = zeros(shape(u))
    uRp = zeros(shape(u))
    uLm = zeros(shape(u))
    uRm = zeros(shape(u))
    #uLp[1::,:] = u[1::,:] + 0.5*phi[1::,:]*(u[1::,:] - u[0:-1,:])
    #uLm[1::,:] = u[0:-1,:] + 0.5*phi[0:-1,:]*(u[1::,:] - u[0:-1,:])
    #uRp[0:-1,:] = u[1::,:] - 0.5*phi[0:-1,:]*(u[1::,:] - u[0:-1,:])
    #uRm[0:-1,:] = u[0:-1,:] - 0.5*phi[0:-1,:]*(u[1::,:] - u[0:-1,:])
    uLp[0:-1,:] = u[0:-1,:] + 0.5*phi[0:-1,:]*(u[1::,:] - u[0:-1,:])
    uLm[1::,:] = u[0:-1,:] + 0.5*phi[0:-1,:]*(u[1::,:] - u[0:-1,:])
    uRp[0:-2,:] = u[1:-1,:] - 0.5*phi[1:-1,:]*(u[2::,:] - u[1:-1,:])
    uRm[0:-1,:] = u[0:-1,:] - 0.5*phi[0:-1,:]*(u[1::,:] - u[0:-1,:])
    if (grid.xperd == 1):
        uLm[0,:] = u[-1,:] + 0.5*phi[-1,:]*(u[0,:] - u[-1,:])
        uLp[-1,:] = u[-1,:] + 0.5*phi[-1,:]*(u[0,:] - u[-1,:])
    
        uRp[-2,:] = u[-1,:] - 0.5*phi[-1,:]*(u[0,:] - u[-1,:])
        uRp[-1,:] = u[0,:] - 0.5*phi[0,:]*(u[1,:] - u[0,:])
        uRm[-1,:] = u[-1,:] - 0.5*phi[-1,:]*(u[0,:] - u[-1,:])
    else:
      uLm[0,:] = u[0,:]
      uLp[-1,:] = u[-1,:]
      uRp[-2,:] = u[-1,:] - 0.5*phi[-1,:]*(u[-2,:] - u[-2,:])
      uRp[-1,:] = u[-1,:]
      uRm[-1,:] = u[-1,:]
    return uLm,uLp,uRm,uRp
    
def linearReconstruct_y(u,grid):
    ry = zeros(shape(u))
    ry[:,1:-1] = (u[:,1:-1] - u[:,0:-2])/(u[:,2::] - u[:,1:-1]  + 1.e-10)
    ry[:,0] = ry[:,1]
    ry[:,-1] = ry[:,-2]
    phi = fluxLimiter(ry)
    uLp = zeros(shape(u))
    uRp = zeros(shape(u))
    uLm = zeros(shape(u))
    uRm = zeros(shape(u))
    #uLp[:,1::] = u[:,1::] + 0.5*phi[:,1::]*(u[:,1::] - u[:,0:-1])
    #uLm[:,1::] = u[:,0:-1] + 0.5*phi[:,0:-1]*(u[:,1::] - u[:,0:-1])
    #uRp[:,0:-1] = u[:,1::] - 0.5*phi[:,0:-1]*(u[:,1::] - u[:,0:-1])
    #uRm[:,0:-1] = u[:,0:-1] - 0.5*phi[:,0:-1]*(u[:,1::] - u[:,0:-1])
    uLp[:,0:-1] = u[:,0:-1] + 0.5*phi[:,0:-1]*(u[:,1::] - u[:,0:-1])
    uLm[:,1::] = u[:,0:-1] + 0.5*phi[:,0:-1]*(u[:,1::] - u[:,0:-1])
    uRp[:,0:-2] = u[:,1:-1] - 0.5*phi[:,1:-1]*(u[:,2::] - u[:,1:-1])
    uRm[:,0:-1] = u[:,0:-1] - 0.5*phi[:,0:-1]*(u[:,1::] - u[:,0:-1])
    if (grid.yperd == 1):
        uLm[:,0] = u[:,-1] + 0.5*phi[:,-1]*(u[:,0] - u[:,-1])
        uLp[:,-1] = u[:,-1] + 0.5*phi[:,-1]*(u[:,0] - u[:,-1])
    
        uRp[:,-2] = u[:,-1] - 0.5*phi[:,-1]*(u[:,0] - u[:,-1])
        uRp[:,-1] = u[:,0] - 0.5*phi[:,0]*(u[:,1] - u[:,0])
        uRm[:,-1] = u[:,-1] - 0.5*phi[:,-1]*(u[:,0] - u[:,-1])
    else:
      uLm[:,0] = uLp[:,0]
      uLp[:,-1] = u[:,-1]
      uRp[:,-2] = u[:,-1] - 0.5*phi[:,-1]*(u[:,-2] - u[:,-2])
      uRp[:,-1] = u[:,-1]
      uRm[:,-1] = u[:,-1]
    return uLm,uLp,uRm,uRp
    
    
def linearReconstruct_x2(u,grid):
    rx = zeros(shape(u))
    rx[1:-1,:] = (u[1:-1,:] - u[0:-2,:])/(u[2::,:] - u[1:-1,:] + 1.e-30)
    rx[0,:] = rx[1,:]
    rx[-1,:] = rx[-2,:]
    phi = fluxLimiter(rx)
    uLp = zeros(shape(u))
    uRp = zeros(shape(u))
    uLm = zeros(shape(u))
    uRm = zeros(shape(u))
    uLp[1::,:] = u[1::,:] + 0.5*phi[1::,:]*(u[1::,:] - u[0:-1,:])
    uLm[1::,:] = u[0:-1,:] + 0.5*phi[0:-1,:]*(u[1::,:] - u[0:-1,:])
    uRp[0:-1,:] = u[1::,:] - 0.5*phi[0:-1,:]*(u[1::,:] - u[0:-1,:])
    uRm[0:-1,:] = u[0:-1,:] - 0.5*phi[0:-1,:]*(u[1::,:] - u[0:-1,:])
    if (grid.xperd == 1):
      uLp[0,:] = u[0,:] + 0.5*phi[0,:]*(u[0,:] - u[-1,:])
      uLm[0,:] = u[-1,:] + 0.5*phi[-1,:]*(u[0,:] - u[-1,:])
      uRp[-1,:] = u[0,:] - 0.5*phi[-1,:]*(u[0,:] - u[-1,:])
      uRm[-1,:] = u[-1,:] - 0.5*phi[-1,:]*(u[0,:] - u[-1,:])
    else:
      uLm[0,:] = u[0,:]
      uLp[-1,:] = u[-1,:]
      uRp[-2,:] = u[-1,:] - 0.5*phi[-1,:]*(u[-2,:] - u[-2,:])
      uRp[-1,:] = u[-1,:]
      uRm[-1,:] = u[-1,:]
    return uLm,uLp,uRm,uRp
    
def linearReconstruct_y2(u,grid):
    ry = zeros(shape(u))
    ry[:,1:-1] = (u[:,1:-1] - u[:,0:-2])/(u[:,2::] - u[:,1:-1]  + 1.e-10)
    ry[:,0] = ry[:,1]
    ry[:,-1] = ry[:,-2]
    phi = fluxLimiter(ry)
    uLp = zeros(shape(u))
    uRp = zeros(shape(u))
    uLm = zeros(shape(u))
    uRm = zeros(shape(u))
    uLp[:,1::] = u[:,1::] + 0.5*phi[:,1::]*(u[:,1::] - u[:,0:-1])
    uLm[:,1::] = u[:,0:-1] + 0.5*phi[:,0:-1]*(u[:,1::] - u[:,0:-1])
    uRp[:,0:-1] = u[:,1::] - 0.5*phi[:,0:-1]*(u[:,1::] - u[:,0:-1])
    uRm[:,0:-1] = u[:,0:-1] - 0.5*phi[:,0:-1]*(u[:,1::] - u[:,0:-1])
    if (grid.yperd == 1):
      uLp[:,0] = u[:,0] + 0.5*phi[:,0]*(u[:,0] - u[:,-1])
      uLm[:,0] = u[:,-1] + 0.5*phi[:,-1]*(u[:,0] - u[:,-1])
      uRp[:,-1] = u[:,0] - 0.5*phi[:,-1]*(u[:,0] - u[:,-1])
      uRm[:,-1] = u[:,-1] - 0.5*phi[:,-1]*(u[:,0] - u[:,-1])
    else:
      uLm[:,0] = uLp[:,0]
      uLp[:,-1] = u[:,-1]
      uRp[:,-2] = u[:,-1] - 0.5*phi[:,-1]*(u[:,-2] - u[:,-2])
      uRp[:,-1] = u[:,-1]
      uRm[:,-1] = u[:,-1]
    return uLm,uLp,uRm,uRp