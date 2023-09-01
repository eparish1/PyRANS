from numpy import *
from numpy.linalg import*
######################## FILTERING #############################################################
################################################################################################
##### Function to get the filter matrices
### Filter procedure to remove high frequency content (stabilizes central differences)
### Implemented as outlined in
##   High-Order Schemes for Navier-Stokes Equations:
#    Algorithm and Implementation into FDL3DI
def getFiltMatrices(u,xperd,yperd):
    # filter coefficient alpha->[0,0.5]
    # alpha = 0 is very diffusive, alpha = 0.5 is not diffusive
    from schemes import alpha
    nx,ny = shape(u)
    Afx = zeros((nx,nx))
    Afy = zeros((ny,ny))
    
    
    ### Second order filter
    fill_diagonal(Afx[0:nx,0:nx],1)
    #Afx[0,1] = alpha
    fill_diagonal(Afx[1:nx-1,0:nx-2],alpha)    
    fill_diagonal(Afx[1:nx-1,2:nx],alpha)
    if (xperd == 1):
        Afx[0,1] = alpha
        Afx[0,-1] = alpha
        Afx[-1,0] = alpha
        Afx[-1,-2] = alpha
    #Afx[-1,-2] = alpha
    Afx = inv(Afx)
    
    # If it isn't 1D then apply in y
    if (ny > 3):
      fill_diagonal(Afy[:,:],1)
      #Afy[0,1] = alpha
      fill_diagonal(Afy[1:-1,0:-2],alpha)    
      fill_diagonal(Afy[1:-1,2::],alpha)
      if (yperd == 1):
          Afy[0,-1] = alpha
          Afy[-1,0] = alpha
          Afy[0,1] = alpha
          Afy[-1,-2] = alpha
      #Afy[-1,-2] = alpha
      Afy = inv(Afy) #A is constant so we just invert it and use the same at every tstep
    else:
      fill_diagonal(Afy,1)
    return Afx,Afy


### Function to apply 2nd order filtering
def filt2(u,grid):
    from schemes import alpha
    nx,ny = shape(u)
    rhs1 = zeros((nx,ny))
    rhs2 = zeros((nx,ny))
    rhs3 = zeros((nx,ny))
    a0 = 0.5 + alpha
    a1 = 0.5 + alpha
    
    ### Now filter in the three principle directions
    ## x direction
    rhs1[1:nx-1,:] = a0*u[1:nx-1,:] + a1/2.*(u[2:nx,:] + u[0:nx-2,:])
    if (grid.xperd == 1):
        rhs1[0,:] = a0*u[0,:] + a1/2.*(u[1,:] + u[-1,:])
        rhs1[-1,:] = a0*u[-1,:] + a1/2.*(u[0,:] + u[-2,:]) 
    else:
        rhs1[0,:] =u[0,:]
        rhs1[-1,:] = u[-1,:]
    rhs1 = reshape(rhs1,(nx,ny))
    un1 = dot(grid.Afx,rhs1)  #this is the linear solve, but A is constant so we inverted
    #u = einsum("ij,jlm -> ilm",grid.Afx,rhs1)
    u[:,:] = reshape(un1,(nx,ny))[:,:]
    ## y direction
    if (ny > 3):
      rhs2[:,1:-1] = a0*u[:,1:-1] + a1/2.*(u[:,2::] + u[:,0:-2])
      if (grid.yperd == 1):
        rhs2[:,0] = a0*u[:,0] + a1/2.*(u[:,1] + u[:,-1])
        rhs2[:,-1] = a0*u[:,-1] + a1/2.*(u[:,0] + u[:,-2])
      else:
        rhs2[:,0] = u[:,0]
        rhs2[:,-1] = u[:,-1]
    rhs2= rollaxis(rhs2,1)
    rhs2 = reshape(rhs2,(ny,nx))
    utemp = dot(grid.Afy,rhs2)
    utemp = rollaxis(utemp,1)
    u[:,:] = reshape(utemp,(nx,ny))
    return u


### Function to apply 4th order filtering
def filt4(u,grid):
    from schemes import alpha
    nx,ny = shape(u)
    rhs1 = zeros((nx,ny))
    rhs2 = zeros((nx,ny))
    rhs3 = zeros((nx,ny))
    ## Filter coeffs for interior stencil
    a0 = 5./8 + 3./4.*alpha
    a1 = 0.5 + alpha
    a2 = -1./8. + 0.25*alpha
    
    ## Filter coeffs for second point in stencils
    ab2 = 1./16 + 7./8.*alpha
    bb2 = 0.75 + 0.5*alpha
    cb2 = 3./8. + 0.25*alpha
    db2 = -1./4. + 0.5*alpha
    eb2 = 1./16. - alpha/8.
    
    ## Filter coeffs for edge point in stencis
    ab1 = 15./16. + alpha/16.
    bb1 = 0.25 + 0.75*alpha
    cb1 = -3./8. + 3./8.*alpha
    db1 = 1./4. - alpha/4.
    eb1 = -1./16. + alpha/16.
    
    a0p = 0.5 + alpha
    a1p = 0.5 + alpha
    
    ### Now filter in the three principle directions
    ## x direction

    rhs1[2:-2,:] = a0*u[2:-2,:] + 0.5*a1*(u[3:-1,:] + u[1:-3,:]) + 0.5*a2*(u[4::,:] + u[0:-4,:])
    if (grid.xperd == 1):
        rhs1[1,:] = a0*u[1,:] + 0.5*a1*(u[2,:] + u[0,:]) + 0.5*a2*(u[3,:] + u[-1,:])
        rhs1[0,:] = a0*u[0,:] + 0.5*a1*(u[1,:] + u[-1,:]) + 0.5*a2*(u[2,:] + u[-2,:])
        rhs1[-2,:] = a0*u[-2,:] + 0.5*a1*(u[-1,:] + u[-3,:]) + 0.5*a2*(u[0,:] + u[-4,:])
        rhs1[-1,:] = a0*u[-1,:] + 0.5*a1*(u[0,:] + u[-2,:]) + 0.5*a2*(u[1,:] + u[-3,:])
    else: 
        rhs1[0,:] = u[0,:]
        rhs1[1,:] = a0p*u[1,:] + a1p/2.*(u[2,:] + u[0,:])
        rhs1[-2,:] = a0p*u[-2,:] + a1p/2.*(u[-1,:] + u[-3,:])
        rhs1[-1,:] = u[-1,:]   
    rhs1 = reshape(rhs1,(nx,ny))
    un1 = dot(grid.Afx,rhs1)  #this is the linear solve, but A is constant so we inverted
    #u = einsum("ij,ilm -> jlm",grid.Afx,rhs1)
    u[:,:] = reshape(un1,(nx,ny))[:,:]
    
    ## y direction
    if (ny > 3):
        rhs2[:,2:-2] = a0*u[:,2:-2] + 0.5*a1*(u[:,3:-1] + u[:,1:-3]) + 0.5*a2*(u[:,4::] + u[:,0:-4])
        if (grid.yperd == 1):
            rhs2[:,1] = a0*u[:,1] + 0.5*a1*(u[:,2] + u[:,0]) + 0.5*a2*(u[:,3] + u[:,-1])
            rhs2[:,0] = a0*u[:,0] + 0.5*a1*(u[:,1] + u[:,-1]) + 0.5*a2*(u[:,2] + u[:,-2])
            rhs2[:,-2] = a0*u[:,-2] + 0.5*a1*(u[:,-1] + u[:,-3]) + 0.5*a2*(u[:,0] + u[:,-4])
            rhs2[:,-1] = a0*u[:,-1] + 0.5*a1*(u[:,0] + u[:,-2]) + 0.5*a2*(u[:,1] + u[:,-3])
        else:
          rhs2[:,0] = u[:,0]
          rhs2[:,1] = a0p*u[:,1] + a1p/2.*(u[:,2] + u[:,0])
          rhs2[:,-2] = a0p*u[:,-2] + a1p/2.*(u[:,-1] + u[:,-3])
          rhs2[:,-1] = u[:,-1]
        
#        rhs2[:,0] =    a0*u[:,0] +    0.5*a1*(u[:,1]    + u[:,-1])     + 0.5*a2*(u[:,2] + u[:,-2])   
#        rhs2[:,1] =    a0*u[:,1] +    0.5*a1*(u[:,2]    + u[:,0])     + 0.5*a2*(u[:,3] + u[:,-1])  
#        rhs2[:,-2] =    a0*u[:,-2] +    0.5*a1*(u[:,-1]    + u[:,-3])     + 0.5*a2*(u[:,0] + u[:,-4])   
#        rhs2[:,-1] =    a0*u[:,-1] +    0.5*a1*(u[:,0]    + u[:,-2])     + 0.5*a2*(u[:,1] + u[:,-3]) 
        rhs2= rollaxis(rhs2,1)
#        rhs2 = reshape(rhs2,(ny,nx))
        utemp = dot(grid.Afy,rhs2)
        utemp = rollaxis(utemp,1)
        u[:,:] = reshape(utemp,(nx,ny))
    return u

############################################################################################
    


###########################################################################################
########################### COMPACT FINITE DIFFERENCING ###################################
### Implemented as outlined in
##   High-Order Schemes for Navier-Stokes Equations:
#    Algorithm and Implementation into FDL3DI
### Added by Eric 7/20/15
######### Get matrices
def getCompactMat(nx,ny,xperd,yperd):
  alpha = 0.25
  alpha1 = 1.
  Ax = zeros((nx,nx))
  Ay = zeros((ny,ny))
  # Build matrix 
  fill_diagonal(Ax,1.)
  fill_diagonal(Ax[1:-1,2::],alpha)
  fill_diagonal(Ax[1:-1,0:-2],alpha)
  if (xperd == 1):
      Ax[0,1] = alpha
      Ax[0,-1] = alpha
      Ax[-1,-2] = alpha
      Ax[-1,0] = alpha
  else:
      Ax[0,1] = alpha1
      Ax[-1,-2] = alpha1

  fill_diagonal(Ay,1.)
  fill_diagonal(Ay[1:-1,2::],alpha)
  fill_diagonal(Ay[1:-1,0:-2],alpha)
  Ay[0,1] = alpha1
  Ay[-1,-2] = alpha1
  if (yperd == 1):
      Ay[0,1] = alpha
      Ay[0,-1] = alpha
      Ay[-1,-2] = alpha
      Ay[-1,0] = alpha
  else:
      Ay[0,1] = alpha1
      Ay[-1,-2] = alpha1
  Ax = inv(Ax)
  Ay = inv(Ay)
  return Ax,Ay

def d_dx_compact4(u,grid):
      ### Algorithm to compute 4th order compact differences
      ### Added by Eric 7/20/15
      nx,ny = shape(u)
      d_dx = zeros((nx,ny))
      RHS = zeros((nx,ny))
      a = 1.5
      a1 = -2.
      b1 = 2
      c1 = 0.
      d1 = 0.

      # Compute the RHS for the x derivatives
      RHS[1:-1,:] = 0.5*a*(u[2::,:] - u[0:-2,:])
      if (grid.xperd == 1):
        RHS[0,:] = 0.5*a*(u[1,:] - u[-1,:])
        RHS[-1,:] = 0.5*a*(u[0,:] - u[-2,:])
      else:
        RHS[0,:] = a1*u[0,:] + b1*u[1,:] + c1*u[2,:] + d1*u[3,:]
        RHS[-1,:] = -a1*u[-1,:] - b1*u[-2,:] - c1*u[-3,:] - d1*u[-4,:]
      RHS = reshape(RHS,(nx,ny))
      #d_dx[:,:] = einsum("ii,ij -> ij",grid.Cx,RHS)
      d_dx[:,:] = dot(grid.Cx,RHS)
      return d_dx
def d_dy_compact4(u,grid):
      nx,ny = shape(u)
      d_dy = zeros((nx,ny))
      RHS = zeros((nx,ny))
      a = 1.5
      a1 = -2.
      b1 = 2
      c1 = 0.
      d1 = 0.
      # Compute the RHS for the x derivatives
      RHS[:,1:-1] = 0.5*a*(u[:,2::] - u[:,0:-2])
      if (grid.yperd == 1):
        RHS[:,0] = 0.5*a*(u[:,1] - u[:,-1])
        RHS[:,-1] = 0.5*a*(u[:,0] - u[:,-2])
      else:
        RHS[:,0] = a1*u[:,0] + b1*u[:,1] + c1*u[:,2] + d1*u[:,3]
        RHS[:,-1] = -a1*u[:,-1] - b1*u[:,-2] - c1*u[:,-3] - d1*u[:,-4]
      temp1 = rollaxis(RHS,1)
      temp = dot(grid.Cy,temp1)
      #d_dy[:,:] = reshape(einsum("jj,ji -> ji",grid.Cy,rollaxis(RHS,1)),(nx,ny))
      d_dy = rollaxis(temp,1)[:,:]
      return d_dy
          
def d_dXC(u,grid):
  ### Algorithm to compute 4th order compact differences
  ### Added by Eric 7/20/15
  nx,ny,nz = shape(u)
  d_dX = zeros((nx,ny,nz,3))
  a = 1.5
  a1 = -17./6.
  b1 = 3./2.
  c1 = 3./2.
  d1 = -1./6.
# Compute the RHS for the x derivatives
  RHS = zeros((nx,ny,nz))
  RHS[1:-1,:,:] = 0.5*a*(u[2::,:,:] - u[0:-2,:,:])
  RHS[0,:,:] = a1*u[0,:,:] + b1*u[1,:,:] + c1*u[2,:,:] + d1*u[3,:,:]
  RHS[-1,:,:] = -a1*u[-1,:,:] - b1*u[-2,:,:] - c1*u[-3,:,:] - d1*u[-4,:,:]
  #d_dX[:,:,:,0] = einsum("ij,jkm -> ikm",Cx,RHS)
  RHS = reshape(RHS,(nx,ny*nz))
  deriv = dot(grid.Cx,RHS)
  d_dX[:,:,:,0] = reshape(deriv,(nx,ny,nz))

# Compute the RHS for the y derivatives
  if (ny > 3):
    RHS = zeros((nx,ny,nz))
    RHS[:,1:-1,:] = 0.5*a*(u[:,2::,:] - u[:,0:-2,:])
    RHS[:,0,:] = a1*u[:,0,:] + b1*u[:,1,:] + c1*u[:,2,:] + d1*u[:,3,:]
    RHS[:,-1,:] = -a1*u[:,-1,:] - b1*u[:,-2,:] - c1*u[:,-3,:] - d1*u[:,-4,:]
  #d_dX[:,:,:,1] = einsum("jk,ikm -> ijm",Cy,RHS)
    RHS = rollaxis(RHS,1)
    RHS = reshape(RHS,(ny,nx*nz))
    deriv = dot(grid.Cy,RHS)
    deriv = rollaxis(deriv,1)
    d_dX[:,:,:,1] = reshape(deriv,(nx,ny,nz))
  else:
    d_dX[:,:,:,1] = 0.
  if (nz > 3):
    # Compute the RHS for the z derivatives
    RHS = zeros((nx,ny,nz))
    RHS[:,:,1:-1] = 0.5*a*(u[:,:,2::] - u[:,:,0:-2])
    RHS[:,:,0] = a1*u[:,:,0] + b1*u[:,:,1] + c1*u[:,:,2] + d1*u[:,:,3]
    RHS[:,:,-1] = -a1*u[:,:,-1] - b1*u[:,:,-2] - c1*u[:,:,-3] - d1*u[:,:,-4]
    d_dX[:,:,:,2] = einsum("kl,ijm -> ijl",grid.Cy,RHS)
  else:
    d_dX[:,:,:,2] = 0.

  ## Apply transformation from curvilinear coords.
    d_dX = einsum("ijklm,ijkm -> ijkl",grid.A,d_dX)
#  d_dXR = zeros((nx,ny,nz,3))
#  d_dXR[:,:,:,0] = A[:,:,:,0,0]*d_dX[:,:,:,0] + A[:,:,:,0,1]*d_dX[:,:,:,1] + A[:,:,:,0,2]*d_dX[:,:,:,2] 
#  d_dXR[:,:,:,1] = A[:,:,:,1,0]*d_dX[:,:,:,0] + A[:,:,:,1,1]*d_dX[:,:,:,1] + A[:,:,:,2,2]*d_dX[:,:,:,2] 
#d_dXR[:,:,:,2] = A[:,:,:,2,0]*d_dX[:,:,:,0] + A[:,:,:,2,1]*d_dX[:,:,:,1] + A[:,:,:,2,2]*d_dX[:,:,:,2] 
  return d_dX


## Function to compute the strain tensor on a YZ plane
def div(u,A):
    nx,ny,nz = shape(u)
    div_u = zeros((nx,ny,nz))
    div_u = sum(d_dX(u,A),axis=3)
    return div_u


######################### DIFFERENCING #######################################################
def d_dX(u,grid):
    ## Compute the derivative of a scalar in x,y,z
    nx,ny,nz = shape(u)
    # derivitives computed by,
    # d/dx = A[0,0]*d/chi + A[0,1]*d/eta + A[0,2]*d/dzet
    # d/dy = A[1,0]*d/chi + A[1,1]*d/eta + A[1,2]*d/dzet
    # d/dz = A[2,0]*d/chi + A[2,1]*d/eta + A[2,2]*d/dzet
    ## If no grid transform, use A  = eye(
    d_dX = zeros((nx,ny,nz,3))
    
    # compute x derivative
    if (grid.xperd == 1):
      d_dX[0,:,:,0] = 0.5*(u[1,:,:] - u[-1,:,:])
      d_dX[1:-1,:,:,0] = 0.5*(u[2::,:,:] - u[0:-2,:,:])
      d_dX[-1,:,:,0] = 0.5*(u[0,:,:] - u[-2,:,:])
    else:
      d_dX[0,:,:,0] = -3./2*u[0,:,:] + 2.*u[1,:,:] - 0.5*u[2,:,:]
      d_dX[1:-1,:,:,0] = 0.5*(u[2::,:,:] - u[0:-2,:,:])
      d_dX[-1,:,:,0] = 3./2*u[-1,:,:] - 2.*u[-2,:,:] + 0.5*u[-3,:,:]
   

    # compute y derivative
    if (ny > 3):
        ## If periodic in y
        if (grid.yperd == 1):
          #d_dX[:,0,:,1] = -3./2*u[:,0,:] + 2.*u[:,1,:] - 0.5*u[:,2,:]
          d_dX[:,0,:,1] = 0.5*( u[:,1,:] - u[:,-1,:] )
          d_dX[:,1:-1,:,1] = 0.5*(u[:,2::,:] - u[:,0:-2,:])
          #d_dX[:,-1,:,1] = 3./2*u[:,-1,:] - 2.*u[:,-2,:] + 0.5*u[:,-3,:]
          d_dX[:,-1,:,1] = 0.5*( u[:,0,:] - u[:,-2,:] )
        else:
          d_dX[:,0,:,1] = -3./2*u[:,0,:] + 2.*u[:,1,:] - 0.5*u[:,2,:]
          d_dX[:,1:-1,:,1] = 0.5*(u[:,2::,:] - u[:,0:-2,:])
          d_dX[:,-1,:,1] = 3./2*u[:,-1,:] - 2.*u[:,-2,:] + 0.5*u[:,-3,:]
            
    else:
        d_dX[:,:,:,1] = 0
    # compute z derivative
    ## only if the problem isn't 2D
    if (nz > 3):
        d_dX[:,:,0,2] = -3./2*u[:,:,0] + 2.*u[:,:,1] - 0.5*u[:,:,2]
        d_dX[:,:,1:-1,2] = 0.5*(u[:,:,2::] - u[:,:,0:-2])
        d_dX[:,:,-1,2] = 3./2*u[:,:,-1] - 2.*u[:,:,-2] + 0.5*u[:,:,-3]
    else:
        d_dX[:,:,:,2] = 0
    # now apply transformation
    # use einsum to perform multidimensional matrix multiplication
    #newVal = zeros((nx,ny,nz,3))
        #for i in range(0,nx):
        #for j in range(0,ny):
        #       for k in range(0,nz):
        #       newVal[i,j,k,:] = dot(A[i,j,k,:,:],d_dX[i,j,k,:])
        
    d_dX = einsum("ijklm,ijkm -> ijkl",grid.A,d_dX)
    return d_dX


######################### DIFFERENCING #######################################################
def bd_dX(u,xperd,yperd):
    ## Compute the derivative of a scalar in x,y
    nx,ny = shape(u)
    d_dX = zeros((nx,ny,2))
    # compute x derivative
    d_dX[1:-1,:,0] = 0.5*(u[2::,:] - u[0:-2,:])
    if (xperd == 0):
      d_dX[-1,:,0] = 3./2*u[-1,:] - 2.*u[-2,:] + 0.5*u[-3,:]
      d_dX[0,:,0] = -3./2*u[0,:] + 2.*u[1,:] - 0.5*u[2,:]
    else:    
      d_dX[-1,:,0] = 0.5*(u[0,:] - u[-2,:]) 
      d_dX[0,:,0] = 0.5*(u[1,:] - u[-1,:]) 

    # compute y derivative
    if (ny > 3):
        d_dX[:,1:-1,1] = 0.5*(u[:,2::] - u[:,0:-2])
        if (yperd == 0):
          d_dX[:,0,1] = -3./2*u[:,0] + 2.*u[:,1] - 0.5*u[:,2]
          d_dX[:,-1,1] = 3./2*u[:,-1] - 2.*u[:,-2] + 0.5*u[:,-3]
        else:
          d_dX[:,0,1] = 0.5*(u[:,1] - u[:,-1])
          d_dX[:,-1,1] = 0.5*(u[:,0] - u[:,-2]) 
    else:
        d_dX[:,:,1] = 0
    return d_dX


def d_dx_central2(u,grid):
    ## Compute the derivative of a scalar in x,y
    nx,ny = shape(u)
    d_dx = zeros((nx,ny))
    # compute x derivative
    d_dx[1:-1,:] = 0.5*(u[2::,:] - u[0:-2,:])
    if (grid.xperd == 0):
      d_dx[0,:] = -3./2*u[0,:] + 2.*u[1,:] - 0.5*u[2,:]
      d_dx[-1,:] = 3./2*u[-1,:] - 2.*u[-2,:] + 0.5*u[-3,:]
    
    else:
      d_dx[0,:] =  0.5*(u[1,:] - u[-1,:])
      d_dx[-1,:] = 0.5*(u[0,:] - u[-2,:])   
    return d_dx

 
def d_dxF1(u,grid):
    ## Compute the derivative of a scalar in x,y
    nx,ny = shape(u)
    d_dx = zeros((nx,ny))
    # compute x derivative
    d_dx[0:-1,:] = (u[1::,:] - u[0:-1,:])
    if (grid.xperd == 1):
      d_dx[-1,:] = (u[0,:] - u[-1,:])
    
    else:
      d_dx[-1,:] =d_dx[-2,:]
    return d_dx
    


def d_dxB1(u,grid):
    ## Compute the derivative of a scalar in x,y
    nx,ny = shape(u)
    d_dx = zeros((nx,ny))
    # compute x derivative
    d_dx[1::,:] = (u[1::,:] - u[0:-1,:])
    if (grid.xperd == 1):
      d_dx[0,:] = (u[0,:] - u[-1,:])
    
    else:
      d_dx[0,:] =d_dx[1,:]
    return d_dx


def d_dyF1(u,grid):
    ## Compute the derivative of a scalar in x,y
    nx,ny = shape(u)
    d_dy = zeros((nx,ny))
    # compute x derivative
    d_dy[:,0:-1] = (u[:,1::] - u[:,0:-1])
    if (grid.yperd == 1):
      d_dy[:,-1] = (u[:,0] - u[:,-1])
    
    else:
      d_dy[:,-1] =d_dy[:,-2]
    return d_dy
    
    
def d_dyB1(u,grid):
    ## Compute the derivative of a scalar in x,y
    nx,ny = shape(u)
    d_dy = zeros((nx,ny))
    # compute x derivative
    d_dy[:,1::] = (u[:,1::] - u[:,0:-1])
    if (grid.yperd == 1):
      d_dy[:,0] = (u[:,0] - u[:,-1])
    
    else:
      d_dy[:,0] =d_dy[:,1]
    return d_dy

######### Second Order Forward and Backward Schemes
def d_dxF2(u,grid):
    nx,ny = shape(u)
    d_dx = zeros((nx,ny))
    d_dx[0:-2,:] = -1.5*u[0:-2,:] + 2.*u[1:-1,:] - 0.5*u[2::,:]
    if (grid.xperd == 1):
      d_dx[-2,:] =  -1.5*u[-2,:] + 2.*u[-1,:] - 0.5*u[0,:]
      d_dx[-1,:] =  -1.5*u[-1,:] + 2.*u[0,:] - 0.5*u[1,:]
    else:
      d_dx[-2,:] = u[-1,:] - u[-2,:]
      d_dx[-1,:] = d_dx[-2,:]
    return d_dx

def d_dyF2(u,grid):
    nx,ny = shape(u)
    d_dy = zeros((nx,ny))
    d_dy[:,0:-2] = -1.5*u[:,0:-2] + 2.*u[:,1:-1] - 0.5*u[:,2::]
    if (grid.xperd == 1):
      d_dy[:,-2] =  -1.5*u[:,-2] + 2.*u[:,1] - 0.5*u[:,0]
      d_dy[:,-1] =  -1.5*u[:,-1] + 2.*u[:,0] - 0.5*u[:,1]
    else:
      d_dy[:,-2] = u[:,-1] - u[:,-2]
      d_dy[:,-1] = d_dy[:,-2]
    return d_dy


def d_dxB2(u,grid):
    ## Compute the derivative of a scalar in x,y
    nx,ny = shape(u)
    d_dx = zeros((nx,ny))
    # compute x derivative
    d_dx[2::,:] = 1.5*u[2::,:] - 2.*u[1:-1,:] +  0.5*u[0:-2,:]
    if (grid.xperd == 1):
      d_dx[1,:] = 1.5*u[1,:] - 2.*u[0,:] +  0.5*u[-1,:]
      d_dx[0,:] = 1.5*u[0,:] - 2.*u[-1,:] +  0.5*u[-2,:]
    else:
      d_dx[1:,] = u[1,:] - u[0,:]
      d_dx[0,:] = d_dx[1,:]
    return d_dx
    
    
def d_dyB2(u,grid):
    ## Compute the derivative of a scalar in x,y
    nx,ny = shape(u)
    d_dy = zeros((nx,ny))
    # compute x derivative
    d_dy[:,2::] = 1.5*u[:,2::] -2.*u[:,1:-1] +  0.5*u[:,0:-2]
    if (grid.yperd == 1):
      d_dy[:,1] = 1.5*u[:,1] - 2.*u[:,0] +  0.5*u[:,-1]
      d_dy[:,0] = 1.5*u[:,0] - 2.*u[:,-1] +  0.5*u[:,-2]
    
    else:
      d_dy[:,1] = u[:,1] - u[:,0]
      d_dy[:,0] = d_dy[:,1]
    return d_dy



######### Third Order Forward and Backward Schemes  
def d_dxF3(u,grid):
    nx,ny = shape(u)
    d_dx = zeros((nx,ny))
    d_dx[0:-3,:] = -11./6.*u[0:-3,:] + 3.*u[1:-2,:] - 1.5*u[2:-1,:] + 1./3.*u[3::,:] 
    if (grid.xperd == 1):
        d_dx[-3,:] = -11./6.*u[-3,:] + 3.*u[-2,:] - 1.5*u[-1,:] + 1./3.*u[0,:] 
        d_dx[-2,:] = -11./6.*u[-2,:] + 3.*u[-1,:] - 1.5*u[0,:] + 1./3.*u[1,:] 
        d_dx[-1,:] = -11./6.*u[-1,:] + 3.*u[0,:] - 1.5*u[1,:] + 1./3.*u[2,:] 
    else:
        d_dx[-3,:] = -3./2.*u[-3,:] + 2.*u[-2,:] - 0.5*u[-1,:]
        d_dx[-2,:] = -1.*u[-2,:] + u[-1,:]
        d_dx[-1,:] = d_dx[-2,:]
    return d_dx

def d_dyF3(u,grid):
    nx,ny = shape(u)
    d_dy = zeros((nx,ny))
    d_dy[:,0:-3] = -11./6.*u[:,0:-3] + 3.*u[:,1:-2] - 1.5*u[:,2:-1] + 1./3.*u[:,3::] 
    if (grid.yperd == 1):
        d_dy[:,-3] = -11./6.*u[:,-3] + 3.*u[:,-2] - 1.5*u[:,-1] + 1./3.*u[:,0] 
        d_dy[:,-2] = -11./6.*u[:,-2] + 3.*u[:,-1] - 1.5*u[:,0] + 1./3.*u[:,1] 
        d_dy[:,-1] = -11./6.*u[:,-1] + 3.*u[:,0] - 1.5*u[:,1] + 1./3.*u[:,2] 
    else:
        d_dy[:,-3] = -3./2.*u[-3,:] + 2.*u[-2,:] - 0.5*u[-1,:]
        d_dy[:,-2] = -1.*u[-2,:] + u[-1,:]
        d_dy[:,-1] = d_dy[-2,:]
    return d_dy

def d_dxB3(u,grid):
    nx,ny = shape(u)
    d_dx = zeros((nx,ny))
    d_dx[3::,:] = 11./6.*u[3::,:] - 3.*u[2:-1,:] + 3./2.*u[1:-2,:] - 1./3.*u[0:-3,:]
    if (grid.xperd == 1):
      d_dx[2,:] = 11./6.*u[2,:] - 3.*u[1,:] + 1.5*u[0,:] - 1./3.*u[-1,:]
      d_dx[1,:] = 11./6.*u[1,:] - 3.*u[0,:] + 1.5*u[-1,:] - 1./3.*u[-2,:]
      d_dx[0,:] = 11./6.*u[0,:] - 3.*u[-1,:] + 1.5*u[-2,:] - 1./3.*u[-3,:]
    else:
      d_dx[2,:] = 3./2.*u[2,:] - 2.*u[1,:] + 1./2.*u[0,:]
      d_dx[1,:] = u[1,:] - u[0,:]
      d_dx[0,:] = d_dx[1,:]
    return d_dx

def d_dyB3(u,grid):
    nx,ny = shape(u)
    d_dy = zeros((nx,ny))
    d_dy[:,3::] = 11./6.*u[:,3::] - 3.*u[:,2:-1] + 3./2.*u[:,1:-2] - 1./3.*u[:,0:-3]
    if (grid.yperd == 1):
      d_dy[:,2] = 11./6.*u[:,2] - 3.*u[:,1] + 1.5*u[:,0] - 1./3.*u[:,-1]
      d_dy[:,1] = 11./6.*u[:,1] - 3.*u[:,0] + 1.5*u[:,-1] - 1./3.*u[:,-2]
      d_dy[:,0] = 11./6.*u[:,0] - 3.*u[:,-1] + 1.5*u[:,-2] - 1./3.*u[:,-3]
    else:
      d_dy[:,2] = 3./2.*u[:,2] - 2.*u[:,1] + 1./2.*u[:,0]
      d_dy[:,1] = u[:,1] - u[:,0]
      d_dy[:,0] = d_dy[:,1]
    return d_dy
######################################################################################




def d_dy_central2(u,grid):
    ## Compute the derivative of a scalar in x,y
    nx,ny = shape(u)
    d_dy = zeros((nx,ny))
    d_dy[:,1:-1] = 0.5*(u[:,2::] - u[:,0:-2])
    if (grid.yperd == 0):
      d_dy[:,0] = -3./2*u[:,0] + 2.*u[:,1] - 0.5*u[:,2]
      d_dy[:,-1] = 3./2*u[:,-1] - 2.*u[:,-2] + 0.5*u[:,-3]
    else:
      d_dy[:,0] = 0.5*(u[:,1] - u[:,-1])
      d_dy[:,-1] = 0.5*(u[:,0] - u[:,-2])
    return d_dy


def d_dxC(u,grid):
  ### Algorithm to compute 4th order compact differences
  ### Added by Eric 7/20/15
  nx,ny = shape(u)
  d_dx = zeros((nx,ny))
  a = 1.5
  a1 = -17./6.
  b1 = 3./2.
  c1 = 3./2.
  d1 = -1./6.
# Compute the RHS for the x derivatives
  RHS = zeros((nx,ny))
  RHS[1:-1,:] = 0.5*a*(u[2::,:] - u[0:-2,:])
  RHS[0,:] = a1*u[0,:] + b1*u[1,:] + c1*u[2,:] + d1*u[3,:]
  RHS[-1,:] = -a1*u[-1,:] - b1*u[-2,:] - c1*u[-3,:] - d1*u[-4,:]
  #d_dX[:,:,:,0] = einsum("ij,jkm -> ikm",Cx,RHS)
  d_dX = dot(grid.Cx,RHS)


def d_dyC(u,grid):
  ### Algorithm to compute 4th order compact differences
  ### Added by Eric 7/20/15
  nx,ny = shape(u)
  d_dy = zeros((nx,ny))
  a = 1.5
  a1 = -17./6.
  b1 = 3./2.
  c1 = 3./2.
  d1 = -1./6.
# Compute the RHS for the x derivatives
  RHS = zeros((nx,ny))
  RHS[:,1:-1] = 0.5*a*(u[:,2::] - u[:,0:-2])
  RHS[:,0] = a1*u[:,0] + b1*u[:,1] + c1*u[:,2] + d1*u[:,3]
  RHS[:,-1] = -a1*u[:,-1] - b1*u[:,-2] - c1*u[:,-3] - d1*u[:,-4]
  #d_dX[:,:,:,0] = einsum("ij,jkm -> ikm",Cx,RHS)
  d_dX = dot(grid.Cx,RHS,axis=1)


# Function to compute the metrics for 2D curvilinear grids
def gridMetrics_central2(x,y):
    from mesh import xperd as xperd1
    from mesh import yperd as yperd1
    from mesh import xperd_grid 
    from mesh import yperd_grid
    #Get mesh sizes
    nx,ny = shape(x)
    #Use uniform spacing in computational domain
    dchi = 1.
    deta = 1.
    #Compute the derivatives of x,y,z w.r.p to chi,eta,zet
    #dx_dCOMP = bd_dX(x,xperd_grid,yperd_grid)
    #dy_dCOMP = bd_dX(y,xperd_grid,yperd_grid)
    dx_dCOMP = bd_dX(x,xperd_grid,yperd_grid)
    dy_dCOMP = bd_dX(y,xperd_grid,yperd_grid)
    #Put the data into a transformation matrix
    #this is the transformation that results from chain rule
    #i.e d/dchi = d/dx*dx/dchi + d/dy*dy/dchi + d/dz*d/dchi
    #    d/deta = d/dx*dx/deta + d/dy*dy/deta + d/dz*d/deta
    #    d/dzet = d/dx*dx/dzet + d/dy*dy/dzet + d/dz*d/dzet
    J = zeros((nx,ny))
    J[:,:] = 1./(dx_dCOMP[:,:,0]*dy_dCOMP[:,:,1] - dx_dCOMP[:,:,1]*dy_dCOMP[:,:,0] )
    xi_x = J*dy_dCOMP[:,:,1]
    eta_x = J*-dy_dCOMP[:,:,0]
    xi_y = J*-dx_dCOMP[:,:,1]
    eta_y = J*dx_dCOMP[:,:,0]
    xit = zeros(shape(x))
    XI = sqrt(xi_x**2 + xi_y**2)
    ETA = sqrt(eta_x**2 + eta_y**2)
    etat = zeros(shape(x))
    return J,xi_x,eta_x,xi_y,eta_y,XI,ETA,xit,etat

def gridMetrics_compact4(x,y):
    from mesh import xperd 
    from mesh import yperd 
    from mesh import xperd_grid 
    from mesh import yperd_grid
    #Get mesh sizes
    nx,ny = shape(x)
    #Use uniform spacing in computational domain
    dchi = 1.
    deta = 1.
    #Compute the derivatives of x,y,z w.r.p to chi,eta,zet
    Cxgrid,Cygrid = getCompactMat(nx,ny,xperd_grid,yperd_grid)
    dx_dxi = d_dxc_grid(x,xperd_grid,Cxgrid)
    dx_deta = d_dyc_grid(x,yperd_grid,Cygrid)
    dy_dxi = d_dxc_grid(y,xperd_grid,Cxgrid)
    dy_deta = d_dyc_grid(y,yperd_grid,Cygrid)

    J = zeros((nx,ny))
    J[:,:] = 1./(dx_dxi[:,:]*dy_deta[:,:] - dx_deta[:,:]*dy_dxi[:,:] )
    xi_x = J*dy_deta[:,:]
    eta_x = J*-dy_dxi[:,:]
    xi_y = J*-dx_deta[:,:]
    eta_y = J*dx_dxi[:,:]
    xit = zeros(shape(x))
    XI = sqrt(xi_x**2 + xi_y**2)
    ETA = sqrt(eta_x**2 + eta_y**2)
    etat = zeros(shape(x))
    return J,xi_x,eta_x,xi_y,eta_y,XI,ETA,xit,etat
    
def d_dxc_grid(u,xperd_grid,Cx):
      ### Algorithm to compute 4th order compact differences
      ### Added by Eric 7/20/15
      nx,ny = shape(u)
      d_dx = zeros((nx,ny))
      RHS = zeros((nx,ny))
      a = 1.5
      a1 = -2.
      b1 = 2
      c1 = 0.
      d1 = 0.
      # Compute the RHS for the x derivatives
      RHS[1:-1,:] = 0.5*a*(u[2::,:] - u[0:-2,:])
      if (xperd_grid == 1):
        RHS[0,:] = 0.5*a*(u[1,:] - u[-1,:])
        RHS[-1,:] = 0.5*a*(u[0,:] - u[-2,:])
      else:
        RHS[0,:] = a1*u[0,:] + b1*u[1,:] + c1*u[2,:] + d1*u[3,:]
        RHS[-1,:] = -a1*u[-1,:] - b1*u[-2,:] - c1*u[-3,:] - d1*u[-4,:]
      #RHS = reshape(RHS,(nx,ny))
      d_dx[:,:] = dot(Cx,RHS)
      return d_dx
      
def d_dyc_grid(u,yperd_grid,Cy):
      nx,ny = shape(u)
      d_dx = zeros((nx,ny))
      RHS = zeros((nx,ny))
      a = 1.5
      a1 = -2.
      b1 = 2
      c1 = 0.
      d1 = 0.
      # Compute the RHS for the x derivatives
      RHS[:,1:-1] = 0.5*a*(u[:,2::] - u[:,0:-2])
      if (yperd_grid == 1):
        RHS[:,0] = 0.5*a*(u[:,1] - u[:,-1])
        RHS[:,-1] = 0.5*a*(u[:,0] - u[:,-2])
      else:
        RHS[:,0] = a1*u[:,0] + b1*u[:,1] + c1*u[:,2] + d1*u[:,3]
        RHS[:,-1] = -a1*u[:,-1] - b1*u[:,-2] - c1*u[:,-3] - d1*u[:,-4]
      temp1 = rollaxis(RHS,1)
      temp = dot(Cy,temp1)
      d_dy = rollaxis(temp,1)[:,:]
      return d_dy



def gradedSpace(x0,xF,totalPoints,cell_rat):
    totalLength = xF-x0
    nx = totalPoints
    x= zeros(totalPoints)
    dx = zeros(totalPoints-1)
    if (cell_rat == 1.0):
        r = 1
        dx0 = totalLength/(totalPoints-1.)
    else:
        r = cell_rat**(1./(nx-2))
        L = totalLength
        dx0 = L/((1-r**(nx-1))/(1-r))
    x[0] = x0
    for i in range(1,nx):
        dx[i-1] = dx0*r**float(i-1)
        x[i] = x[i-1]+dx[i-1]
        
    dxEnd = x[nx-1] - x[nx-3]
    x[nx-2] = x[nx-3] + 0.5*dxEnd
    return x

### cylinder mesh
#r = linspace(1,3,100)
#theta = linspace(0,2*pi,100)
#x = zeros((100,100,10))
#y = zeros((100,100,10))
#z = zeros((100,100,10))
#for i in range(0,100):
#    for j in range(0,100):
#        for k in range(0,10):
#          x[i,j,k] = r[i]*cos(theta[j])
#          y[i,j,k] = r[i]*sin(theta[j])
#          z[i,j,k] = k

