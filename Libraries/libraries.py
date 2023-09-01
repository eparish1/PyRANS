from pylab import *

######################## FILTERING #############################################################
################################################################################################
##### Function to get the filter matrices
### Filter procedure to remove high frequency content (stabilizes central differences)
### Implemented as outlined in
##   High-Order Schemes for Navier-Stokes Equations:
#    Algorithm and Implementation into FDL3DI
def getFiltMatrices(u):
    # filter coefficient alpha->[0,0.5]
    # alpha = 0 is very diffusive, alpha = 0.5 is not diffusive
    alpha = 0.475
    nx,ny,nz = shape(u)
    Afx = zeros((nx,nx))
    Afy = zeros((ny,ny))
    Afz = zeros((nz,nz))
    
    
    ### Second order filter
    fill_diagonal(Afx[0:nx,0:nx],1)
    #Afx[0,1] = alpha
    fill_diagonal(Afx[1:nx-1,0:nx-2],alpha)    
    fill_diagonal(Afx[1:nx-1,2:nx],alpha)
    #Afx[-1,-2] = alpha
    Afx = inv(Afx)
    
    # If it isn't 1D then apply in y
    if (ny > 3):
      fill_diagonal(Afy[:,:],1)
      #Afy[0,1] = alpha
      fill_diagonal(Afy[1:-1,0:-2],alpha)    
      fill_diagonal(Afy[1:-1,2::],alpha)
      #Afy[-1,-2] = alpha
      Afy = inv(Afy) #A is constant so we just invert it and use the same at every tstep
    else:
      fill_diagonal(Afy,1)

    # if it isn't 2D then apply in z
    if (nz > 3):
      fill_diagonal(Afz[:,:],1)
      #Afz[0,1] = alpha
      fill_diagonal(Afz[1:-1,0:-2],alpha)    
      fill_diagonal(Afz[1:-1,2::],alpha)
      #Afz[-1,-2] = alpha
      Afz = inv(Afz)
    else:
      fill_diagonal(Afz,1)
    return Afx,Afy,Afz


### Function to apply 2nd order filtering
def filt(u,grid):
    alpha = 0.475
    nx,ny,nz = shape(u)
    rhs1 = zeros((nx,ny,nz))
    rhs2 = zeros((nx,ny,nz))
    rhs3 = zeros((nx,ny,nz))
    a0 = 0.5 + alpha
    a1 = 0.5 + alpha
    
    ### Now filter in the three principle directions
    ## x direction
    rhs1[0,:,:] =u[0,:,:]
    rhs1[-1,:,:] = u[-1,:,:]
    rhs1[1:nx-1,:,:] = a0*u[1:nx-1,:,:] + a1/2.*(u[2:nx,:,:] + u[0:nx-2,:,:])
    rhs1 = reshape(rhs1,(nx,ny*nz))
    un1 = dot(grid.Afx,rhs1)  #this is the linear solve, but A is constant so we inverted
    #u = einsum("ij,jlm -> ilm",grid.Afx,rhs1)
    u[:,:,:] = reshape(un1,(nx,ny,nz))[:,:,:]
    ## y direction
    if (ny > 3):
       rhs2[:,0,:] = u[:,0,:]
       rhs2[:,-1,:] = u[:,-1,:]
       rhs2[:,1:-1,:] = a0*u[:,1:-1,:] + a1/2.*(u[:,2::,:] + u[:,0:-2,:])
           #u[:,:,:] = einsum("mj,ijk -> imk",grid.Afy,rhs2)[:,:,:] #this is SLOW. Find a way to reshape and use dot function
       rhs2= rollaxis(rhs2,1)
       rhs2 = reshape(rhs2,(ny,nx*nz))
       utemp = dot(grid.Afy,rhs2)
       utemp = rollaxis(utemp,1)
       u[:,:,:] = reshape(utemp,(nx,ny,nz))
    
    if (nz > 3):
      rhs3[:,:,0] = (0.75 + 0.25*alpha)*u[:,:,0] + (0.5 + 0.5*alpha)*u[:,:,1] + (-0.25 + 0.25*alpha)*u[:,:,2]
      rhs3[:,:,-1] = (0.75 + 0.25*alpha)*u[:,:,-1] + (0.5 + 0.5*alpha)*u[:,:,-2] + (-0.25 + 0.25*alpha)*u[:,:,-3]
      rhs3[:,:,1:-1] = a0*u[:,:,1:-1] + a1/2.*(u[:,:,2::] + u[:,:,0:-2])
      u[:,:,:] = einsum("kl,ijk -> ijl",grid.Afz,rhs3)[:,:,:] #again, slow
    return u


### Function to apply 4th order filtering
def filt4(u,grid):
    alpha = 0.475
    nx,ny,nz = shape(u)
    rhs1 = zeros((nx,ny,nz))
    rhs2 = zeros((nx,ny,nz))
    rhs3 = zeros((nx,ny,nz))
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
    rhs1[0,:,:] = u[0,:,:]
    rhs1[1,:,:] = a0p*u[1,:,:] + a1p/2.*(u[2,:,:] + u[0,:,:])
    rhs1[-2,:,:] = a0p*u[-2,:,:] + a1p/2.*(u[-1,:,:] + u[-3,:,:])
    rhs1[-1,:,:] = u[-1,:,:]
    rhs1[2:-2,:,:] = a0*u[2:-2,:,:] + 0.5*a1*(u[3:-1,:,:] + u[1:-3,:,:]) + 0.5*a2*(u[4::,:,:] + u[0:-4,:,:])
    rhs1 = reshape(rhs1,(nx,ny*nz))
    un1 = dot(grid.Afx,rhs1)  #this is the linear solve, but A is constant so we inverted
    #u = einsum("ij,ilm -> jlm",grid.Afx,rhs1)
    u[:,:,:] = reshape(un1,(nx,ny,nz))[:,:,:]
    ## y direction
    if (ny > 3):
        rhs2[:,0,:] = u[:,0,:]
        rhs2[:,1,:] = a0p*u[:,1,:] + a1p/2.*(u[:,2,:] + u[:,0,:])
        rhs2[:,-2,:] = a0p*u[:,-2,:] + a1p/2.*(u[:,-1,:] + u[:,-3,:])
        rhs2[:,-1,:] = u[:,-1,:]
        rhs2[:,2:-2,:] = a0*u[:,2:-2,:] + 0.5*a1*(u[:,3:-1,:] + u[:,1:-3,:]) + 0.5*a2*(u[:,4::,:] + u[:,0:-4,:])
    #u[:,:,:] = einsum("jk,ijl -> ikl",grid.Afy,rhs2)[:,:,:] #this is SLOW. Find a way to reshape and use dot function
        rhs2= rollaxis(rhs2,1)
        rhs2 = reshape(rhs2,(ny,nx*nz))
        utemp = dot(grid.Afy,rhs2)
        utemp = rollaxis(utemp,1)
        u[:,:,:] = reshape(utemp,(nx,ny,nz))
                                
    if (nz > 3):
        rhs3[:,:,0] = ab1*u[:,:,0] + bb1*u[:,:,1] + cb1*u[:,:,2] + db1*u[:,:,3] + eb1*u[:,:,4]
        rhs3[:,:,1] = ab2*u[:,:,0] + bb2*u[:,:,1] + cb2*u[:,:,2] + db2*u[:,:,3] + eb2*u[:,:,4]
        rhs3[:,:,-2] = ab2*u[:,:,-1] + bb2*u[:,:,-2] + cb2*u[:,:,-3] + db2*u[:,:,-4] + eb2*u[:,:,-5]
        rhs3[:,:,-1] = ab1*u[:,:,-1] + bb1*u[:,:,-2] + cb1*u[:,:,-3] + db1*u[:,:,-4] + eb1*u[:,:,-5]
        rhs3[:,:,2:nx-2] = a0*u[:,:,2:nx-2] + 0.5*a1*(u[:,:,3:nx-1] + u[:,:,1:nx-3]) + 0.5*a2*(u[:,:,4:nx] + u[:,:,0:nx-4])
        u[:,:,:] = einsum("kl,ijk -> ijl",grid.Afz,rhs3)[:,:,:] #again, slow
    return u

############################################################################################
    


###########################################################################################
########################### COMPACT FINITE DIFFERENCING ###################################
### Implemented as outlined in
##   High-Order Schemes for Navier-Stokes Equations:
#    Algorithm and Implementation into FDL3DI
### Added by Eric 7/20/15
######### Get matrices
def getCompactMat(u):
  alpha = 0.25
  alpha1 = 3.
  nx,ny,nz = shape(u)
  Ax = zeros((nx,nx))
  Ay = zeros((ny,ny))
  Az = zeros((nz,nz))
  # Build matrix 
  fill_diagonal(Ax,1.)
  fill_diagonal(Ax[1:-1,2::],alpha)
  fill_diagonal(Ax[1:-1,0:-2],alpha)
  Ax[0,1] = alpha1
  Ax[-1,-2] = alpha1

  fill_diagonal(Ay,1.)
  fill_diagonal(Ay[1:-1,2::],alpha)
  fill_diagonal(Ay[1:-1,0:-2],alpha)
  Ay[0,1] = alpha1
  Ay[-1,-2] = alpha1

  if (nz > 3):
    fill_diagonal(Az,1.)
    fill_diagonal(Az[1:-1,2::],alpha)
    fill_diagonal(Az[1:-1,0:-2],alpha)
    Az[0,1] = alpha1
    Az[-1,-2] = alpha1
  else:
    fill_diagonal(Az,1.)
  Ax = inv(Ax)
  Ay = inv(Ay)
  Az = inv(Az)
  return Ax,Ay,Az

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
def bd_dX(u,A):
    ## Compute the derivative of a scalar in x,y,z
    nx,ny,nz = shape(u)
    # derivitives computed by,
    # d/dx = A[0,0]*d/chi + A[0,1]*d/eta + A[0,2]*d/dzet
    # d/dy = A[1,0]*d/chi + A[1,1]*d/eta + A[1,2]*d/dzet
    # d/dz = A[2,0]*d/chi + A[2,1]*d/eta + A[2,2]*d/dzet
    ## If no grid transform, use A  = eye(
    d_dX = zeros((nx,ny,nz,3))
    
    # compute x derivative
    d_dX[0,:,:,0] = -3./2*u[0,:,:] + 2.*u[1,:,:] - 0.5*u[2,:,:]
    d_dX[1:-1,:,:,0] = 0.5*(u[2::,:,:] - u[0:-2,:,:])
    d_dX[-1,:,:,0] = 3./2*u[-1,:,:] - 2.*u[-2,:,:] + 0.5*u[-3,:,:]
    
    # compute y derivative
    if (ny > 3):
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
    return d_dX


def d_dx1(u,nx,ny,nz,dx):
    du_dx = zeros((ny,nz))
    du_dx[:,:] = (u[2,:,:] - u[0,:,:])/dx
    return du_dx

def d_dx2(u,ny,nz,dy):
    du_dy = zeros((ny,nz))
    du_dy[1:-1,:] = (u[1,2:ny:] - u[1,0:ny-2:])/dy
    return du_dy

def d_dx3(u,ny,nz,dz):
    du_dz = zeros((ny,nz))
    du_dz[:,1:-1] = (u[1,:,2:nz] - u[1,:,0:nz-2])/dz
    return du_dz

def computeStrainTensor(u,v,w,nz,ny,dx,dy,dz):
    S_ij = zeros((3,3,ny,nz)) 
    mu = ones((ny,nz))
    du1_dx1 = d_dx1(u,ny,nz,dx)
    du1_dx2 = d_dx2(u,ny,nz,dy)
    du1_dx3 = d_dx3(u,ny,nz,dz)
    
    du2_dx1 = d_dx1(v,ny,nz,dx)
    du2_dx2 = d_dx2(v,ny,nz,dy)
    du2_dx3 = d_dx3(v,ny,nz,dz)
    
    du3_dx1 = d_dx1(w,ny,nz,dx)
    du3_dx2 = d_dx2(w,ny,nz,dy)
    du3_dx3 = d_dx3(w,ny,nz,dz)
    
    S_ij[0,0,:,:] = mu * du1_dx1
    S_ij[0,1,:,:] = mu * 0.5 * ( du1_dx2 + du2_dx1 )
    S_ij[0,2,:,:] = mu * 0.5 * ( du1_dx3 + du3_dx1 )
    
    S_ij[1,0,:,:] = mu * 0.5 * ( du2_dx1 + du1_dx2 )
    S_ij[1,1,:,:] = mu * du2_dx2
    S_ij[1,2,:,:] = mu * 0.5 * ( du2_dx3 + du3_dx2 ) 
    
    S_ij[2,0,:,:] = mu * 0.5 * ( du3_dx1 + du1_dx3 )
    S_ij[2,1,:,:] = mu * 0.5 * ( du3_dx2 + du2_dx3 )
    S_ij[2,2,:,:] = mu * du3_dx3
    
    ## Add divergence term for compressible flow with stokes hypothese
    divergence = (du1_dx1 + du2_dx2 + du3_dx3)
    lam = -2./3.*mu
    S_ij[0,0,:,:] += lam*divergence
    S_ij[1,1,:,:] += lam*divergence
    S_ij[2,2,:,:] += lam*divergence
    return S_ij

# Function to compute the metrics for 3D curvilinear grids
def gridMetrics(x,y,z):
    #Get mesh sizes
    nx,ny,nz = shape(x)
    #Use uniform spacing in computational domain
    dchi = 1.
    deta = 1.
    dzet = 1.
    #Compute the derivatives of x,y,z w.r.p to chi,eta,zet
    Ag = zeros((nx,ny,nz,3,3))
    Ag[:,:,:] = eye(3)
    dx_dCOMP = bd_dX(x,Ag)
    dy_dCOMP = bd_dX(y,Ag)
    dz_dCOMP = bd_dX(z,Ag)
    #Put the data into a transformation matrix
    #this is the transformation that results from chain rule
    #i.e d/dchi = d/dx*dx/dchi + d/dy*dy/dchi + d/dz*d/dchi
    #    d/deta = d/dx*dx/deta + d/dy*dy/deta + d/dz*d/deta
    #    d/dzet = d/dx*dx/dzet + d/dy*dy/dzet + d/dz*d/dzet
    A = zeros((nx,ny,nz,3,3))
    A[:,:,:,0,0] = dx_dCOMP[:,:,:,0] + 1.e-110  ## add small elements to diagonal for 1d/2d
    A[:,:,:,0,1] = dy_dCOMP[:,:,:,0] #+ 1.e-20
    A[:,:,:,0,2] = dz_dCOMP[:,:,:,0] #+ 1.e-20
    A[:,:,:,1,0] = dx_dCOMP[:,:,:,1] #+ 1.e-20
    A[:,:,:,1,1] = dy_dCOMP[:,:,:,1]  + 1.e-110
    A[:,:,:,1,2] = dz_dCOMP[:,:,:,1] #+ 1.e-20
    A[:,:,:,2,0] = dx_dCOMP[:,:,:,2] #+ 1.e-20
    A[:,:,:,2,1] = dy_dCOMP[:,:,:,2] #+ 1.e-20
    A[:,:,:,2,2] = dz_dCOMP[:,:,:,2] + 1.e-80

    #Invert A to get a general transform for d/d chi,eta,zet -> d/d x,,z
    A = inv(A)
    # derivitives can now be computed by,
    # d/dx = A[0,0]*d/chi + A[0,1]*d/eta + A[0,2]*d/dzet
    # d/dy = A[1,0]*d/chi + A[1,1]*d/eta + A[1,2]*d/dzet
    # d/dz = A[2,0]*d/chi + A[2,1]*d/eta + A[2,2]*d/dzet
    return A


# Function to compute the metrics for 3D curvilinear grids
def gridMetricsC(x,y,z):
    Ax,Ay,Az = getCompactMat(x)
    #Get mesh sizes
    nx,ny,nz = shape(x)
    #Use uniform spacing in computational domain
    dchi = 1.
    deta = 1.
    dzet = 1.
    #Compute the derivatives of x,y,z w.r.p to chi,eta,zet
    Ag = zeros((nx,ny,nz,3,3))
    Ag[:,:,:] = eye(3)
    class dummyGrid():
      Cx = Ax
      Cy = Ay
      Cz = Az
      A = Ag
    dummy = dummyGrid()
    dx_dCOMP = d_dXC(x,dummyGrid)
    dy_dCOMP = d_dXC(y,dummyGrid)
    dz_dCOMP = d_dXC(z,dummyGrid)
    #Put the data into a transformation matrix
    #this is the transformation that results from chain rule
    #i.e d/dchi = d/dx*dx/dchi + d/dy*dy/dchi + d/dz*d/dchi
    #    d/deta = d/dx*dx/deta + d/dy*dy/deta + d/dz*d/deta
    #    d/dzet = d/dx*dx/dzet + d/dy*dy/dzet + d/dz*d/dzet
    A = zeros((nx,ny,nz,3,3))
    A[:,:,:,0,0] = dx_dCOMP[:,:,:,0] + 1.e-110  ## add small elements to diagonal for 1d/2d
    A[:,:,:,0,1] = dy_dCOMP[:,:,:,0] #+ 1.e-20
    A[:,:,:,0,2] = dz_dCOMP[:,:,:,0] #+ 1.e-20
    A[:,:,:,1,0] = dx_dCOMP[:,:,:,1] #+ 1.e-20
    A[:,:,:,1,1] = dy_dCOMP[:,:,:,1]  + 1.e-110
    A[:,:,:,1,2] = dz_dCOMP[:,:,:,1] #+ 1.e-20
    A[:,:,:,2,0] = dx_dCOMP[:,:,:,2] #+ 1.e-20
    A[:,:,:,2,1] = dy_dCOMP[:,:,:,2] #+ 1.e-20
    A[:,:,:,2,2] = dz_dCOMP[:,:,:,2] + 1.e-80
    
    #Invert A to get a general transform for d/d chi,eta,zet -> d/d x,,z
    A = inv(A)
    # derivitives can now be computed by,
    # d/dx = A[0,0]*d/chi + A[0,1]*d/eta + A[0,2]*d/dzet
    # d/dy = A[1,0]*d/chi + A[1,1]*d/eta + A[1,2]*d/dzet
    # d/dz = A[2,0]*d/chi + A[2,1]*d/eta + A[2,2]*d/dzet
    return A



def TbTinv(A):
  det = A[0,0]*(A[1,1]*A[2,2] - A[2,1]*A[1,2]) - \
        A[0,1]*(A[1,0]*A[2,2] - A[2,0]*A[1,2]) + \
        A[0,2]*(A[1,0]*A[2,1] - A[2,0]*A[1,1])
  A = A.transpose()
  
  return det


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

