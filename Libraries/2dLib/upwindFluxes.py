from pylab import *
def upwindFluxes(u,v,rho,a,grid1):
    Fp = zeros((4,grid1.nx,grid1.ny))
    Fm = zeros((4,grid1.nx,grid1.ny))
    Gp = zeros((4,grid1.nx,grid1.ny))
    Gm = zeros((4,grid1.nx,grid1.ny))
    gam = 1.4
    
    
    lam1 = -a*sqrt(grid1.eta_y**2+grid1.eta_x**2)+v*grid1.eta_y+u*grid1.eta_x+grid1.eta_t
    lam2 =  a*sqrt(grid1.eta_y**2+grid1.eta_x**2)+v*grid1.eta_y+u*grid1.eta_x+grid1.eta_t
    lam3 =  v*grid1.eta_y+u*grid1.eta_x+grid1.eta_t
    lam4 =  v*grid1.eta_y+u*grid1.eta_x+grid1.eta_t
    lam1p = 0.5*(lam1 + abs(lam1))
    lam2p = 0.5*(lam2 + abs(lam2))
    lam3p = 0.5*(lam3 + abs(lam3))
    lam4p = 0.5*(lam4 + abs(lam4))
    lam1m = 0.5*(lam1 - abs(lam1))
    lam2m = 0.5*(lam2 - abs(lam2))
    lam3m = 0.5*(lam3 - abs(lam3))
    lam4m = 0.5*(lam4 - abs(lam4))
    Vc = grid1.eta_t + grid1.eta_x*u + grid1.eta_y*v
    ## Positive  U Glux
    Gp[0,:,:] = (rho*(2*lam3p*(gam-1)+lam2p+lam1p))/(2*grid1.J*gam)
      
    Gp[1,:,:] = rho*(u*grid1.ETA*(2.*lam4p*(gam-1.) + (lam2p + lam1p)) + (lam2p - lam1p)*a*grid1.eta_x)/(2.*grid1.J*grid1.ETA*gam)

    Gp[2,:,:] = rho*(gam - 1.)/(grid1.eta_y*grid1.J*gam)*(lam3p*(Vc - grid1.eta_t) - lam4p*u*grid1.eta_x) + rho/(2.*grid1.ETA**2*grid1.J*gam)*\
      (a*grid1.eta_y*grid1.ETA*(lam2p - lam1p) + (v*grid1.eta_y**2 + v*grid1.eta_x**2)*(lam2p + lam1p) )
      
#    t1 = lam1p*(rho*(gam - 1.)*(2.*a*grid1.ETA*(grid1.eta_y*v + grid1.eta_x*u ) - (v**2 + u**2)*grid1.ETA) - 4.*a**2*grid1.ETA**2)/(4.*grid1.ETA**2*grid1.J*(gam - 1.)*gam)
#    t2 = lam2p*rho/(4*grid1.ETA**2*grid1.J*(gam-1.)*gam)*((gam - 1.)*(v*(2.*a*grid1.eta_y*grid1.ETA) + v**2*grid1.ETA  + u*(2.*a*grid1.eta_x*grid1.ETA) + u**2*grid1.ETA ) + 2.*a**2*grid1.ETA**2) 
    t1 = -(lam1p*rho*(2*a*v*grid1.eta_y*grid1.ETA*gam+2*a*u*grid1.eta_x*grid1.ETA*gam-v**2*grid1.eta_y**2*gam-u**2*grid1.eta_y**2*gam-v**2*grid1.eta_x**2*gam-u**2*grid1.eta_x**2*gam-2*a*v*grid1.eta_y*grid1.ETA-2*a*u*grid1.eta_x*grid1.ETA+v**2*grid1.eta_y**2+u**2*grid1.eta_y**2-2*a**2*grid1.eta_y**2+v**2*grid1.eta_x**2+u**2*grid1.eta_x**2-2*a**2*grid1.eta_x**2))/(4*(grid1.eta_y**2+grid1.eta_x**2)*grid1.J*(gam-1)*gam)
    t2 = (lam2p*rho*(2*a*v*grid1.eta_y*grid1.ETA*gam+2*a*u*grid1.eta_x*grid1.ETA*gam+v**2*grid1.eta_y**2*gam+u**2*grid1.eta_y**2*gam+v**2*grid1.eta_x**2*gam+u**2*grid1.eta_x**2*gam-2*a*v*grid1.eta_y*grid1.ETA-2*a*u*grid1.eta_x*grid1.ETA-v**2*grid1.eta_y**2-u**2*grid1.eta_y**2+2*a**2*grid1.eta_y**2-v**2*grid1.eta_x**2-u**2*grid1.eta_x**2+2*a**2*grid1.eta_x**2))/(4*(grid1.eta_y**2+grid1.eta_x**2)*grid1.J*(gam-1)*gam)
    t3 = lam3p*rho*(gam -1.)*(grid1.eta_y*(v**2 - u**2) + 2.*u*v*grid1.eta_x)/(2.*grid1.eta_y*grid1.J*gam)
    t4 = lam4p*rho*u*(gam-1.)*(u*grid1.eta_y - v*grid1.eta_x)/(grid1.eta_y*grid1.J*gam) 
    Gp[3,:,:] = t1 + t2 + t3 + t4

    ## Negative U flux
    Gm[0,:,:] = (rho*(2*lam3m*(gam-1.)+lam2m+lam1m))/(2*grid1.J*gam)
      
    Gm[1,:,:] = rho*(u*grid1.ETA*(2.*lam4m*(gam-1.) + (lam2m + lam1m)) + (lam2m - lam1m)*a*grid1.eta_x)/(2.*grid1.J*grid1.ETA*gam)

    Gm[2,:,:] = rho*(gam - 1.)/(grid1.eta_y*grid1.J*gam)*(lam3m*(Vc - grid1.eta_t) - lam4m*u*grid1.eta_x) + rho/(2.*grid1.ETA**2*grid1.J*gam)*\
      (a*grid1.eta_y*grid1.ETA*(lam2m - lam1m) + (v*grid1.eta_y**2 + v*grid1.eta_x**2)*(lam2m + lam1m) )

    t1 = -(lam1m*rho*(2*a*v*grid1.eta_y*grid1.ETA*gam+2*a*u*grid1.eta_x*grid1.ETA*gam-v**2*grid1.eta_y**2*gam-u**2*grid1.eta_y**2*gam-v**2*grid1.eta_x**2*gam-u**2*grid1.eta_x**2*gam-2*a*v*grid1.eta_y*grid1.ETA-2*a*u*grid1.eta_x*grid1.ETA+v**2*grid1.eta_y**2+u**2*grid1.eta_y**2-2*a**2*grid1.eta_y**2+v**2*grid1.eta_x**2+u**2*grid1.eta_x**2-2*a**2*grid1.eta_x**2))/(4*(grid1.eta_y**2+grid1.eta_x**2)*grid1.J*(gam-1)*gam)
    t2 = (lam2m*rho*(2*a*v*grid1.eta_y*grid1.ETA*gam+2*a*u*grid1.eta_x*grid1.ETA*gam+v**2*grid1.eta_y**2*gam+u**2*grid1.eta_y**2*gam+v**2*grid1.eta_x**2*gam+u**2*grid1.eta_x**2*gam-2*a*v*grid1.eta_y*grid1.ETA-2*a*u*grid1.eta_x*grid1.ETA-v**2*grid1.eta_y**2-u**2*grid1.eta_y**2+2*a**2*grid1.eta_y**2-v**2*grid1.eta_x**2-u**2*grid1.eta_x**2+2*a**2*grid1.eta_x**2))/(4*(grid1.eta_y**2+grid1.eta_x**2)*grid1.J*(gam-1)*gam)
    t3 = lam3m*rho*(gam -1.)*(grid1.eta_y*(v**2 - u**2) + 2.*u*v*grid1.eta_x)/(2.*grid1.eta_y*grid1.J*gam)
    t4 = lam4m*rho*u*(gam-1.)*(u*grid1.eta_y - v*grid1.eta_x)/(grid1.eta_y*grid1.J*gam) 
    Gm[3,:,:] = t1 + t2 + t3 + t4


    lam1 = -a*sqrt(grid1.xi_y**2+grid1.xi_x**2)+v*grid1.xi_y+u*grid1.xi_x+grid1.xi_t
    lam2 = a*sqrt(grid1.xi_y**2+grid1.xi_x**2)+v*grid1.xi_y+u*grid1.xi_x+grid1.xi_t
    lam3 = v*grid1.xi_y+u*grid1.xi_x+grid1.xi_t
    lam4 = v*grid1.xi_y+u*grid1.xi_x+grid1.xi_t
    lam1p = 0.5*(lam1 + abs(lam1))
    lam2p = 0.5*(lam2 + abs(lam2))
    lam3p = 0.5*(lam3 + abs(lam3))
    lam4p = 0.5*(lam4 + abs(lam4))
    lam1m = 0.5*(lam1 - abs(lam1))
    lam2m = 0.5*(lam2 - abs(lam2))
    lam3m = 0.5*(lam3 - abs(lam3))
    lam4m = 0.5*(lam4 - abs(lam4))
    Uc = grid1.xi_t + grid1.xi_x*u + grid1.xi_y*v
    ## Positive  U Flux
    Fp[0,:,:] = (rho*(2*lam3p*(gam-1)+lam2p+lam1p))/(2*grid1.J*gam)
      
    Fp[1,:,:] = rho*(u*grid1.XI*(2.*lam4p*(gam-1.) + (lam2p + lam1p)) + (lam2p - lam1p)*a*grid1.xi_x)/(2.*grid1.J*grid1.XI*gam)

    Fp[2,:,:] = rho*(gam - 1.)/(grid1.xi_y*grid1.J*gam)*(lam3p*(Uc - grid1.xi_t) - lam4p*u*grid1.xi_x) + rho/(2.*grid1.XI**2*grid1.J*gam)*\
      (a*grid1.xi_y*grid1.XI*(lam2p - lam1p) + (v*grid1.xi_y**2 + v*grid1.xi_x**2)*(lam2p + lam1p) )
      
#    t1 = lam1p*(rho*(gam - 1.)*(2.*a*grid1.XI*(grid1.xi_y*v + grid1.xi_x*u ) - (v**2 + u**2)*grid1.XI) - 4.*a**2*grid1.XI**2)/(4.*grid1.XI**2*grid1.J*(gam - 1.)*gam)
#    t2 = lam2p*rho/(4*grid1.XI**2*grid1.J*(gam-1.)*gam)*((gam - 1.)*(v*(2.*a*grid1.xi_y*grid1.XI) + v**2*grid1.XI  + u*(2.*a*grid1.xi_x*grid1.XI) + u**2*grid1.XI ) + 2.*a**2*grid1.XI**2) 
    t1 = -(lam1p*rho*(2*a*v*grid1.xi_y*grid1.XI*gam+2*a*u*grid1.xi_x*grid1.XI*gam-v**2*grid1.xi_y**2*gam-u**2*grid1.xi_y**2*gam-v**2*grid1.xi_x**2*gam-u**2*grid1.xi_x**2*gam-2*a*v*grid1.xi_y*grid1.XI-2*a*u*grid1.xi_x*grid1.XI+v**2*grid1.xi_y**2+u**2*grid1.xi_y**2-2*a**2*grid1.xi_y**2+v**2*grid1.xi_x**2+u**2*grid1.xi_x**2-2*a**2*grid1.xi_x**2))/(4*(grid1.xi_y**2+grid1.xi_x**2)*grid1.J*(gam-1)*gam)
    t2 = (lam2p*rho*(2*a*v*grid1.xi_y*grid1.XI*gam+2*a*u*grid1.xi_x*grid1.XI*gam+v**2*grid1.xi_y**2*gam+u**2*grid1.xi_y**2*gam+v**2*grid1.xi_x**2*gam+u**2*grid1.xi_x**2*gam-2*a*v*grid1.xi_y*grid1.XI-2*a*u*grid1.xi_x*grid1.XI-v**2*grid1.xi_y**2-u**2*grid1.xi_y**2+2*a**2*grid1.xi_y**2-v**2*grid1.xi_x**2-u**2*grid1.xi_x**2+2*a**2*grid1.xi_x**2))/(4*(grid1.xi_y**2+grid1.xi_x**2)*grid1.J*(gam-1)*gam)
    t3 = lam3p*rho*(gam -1.)*(grid1.xi_y*(v**2 - u**2) + 2.*u*v*grid1.xi_x)/(2.*grid1.xi_y*grid1.J*gam)
    t4 = lam4p*rho*u*(gam-1.)*(u*grid1.xi_y - v*grid1.xi_x)/(grid1.xi_y*grid1.J*gam) 
    Fp[3,:,:] = t1 + t2 + t3 + t4

    ## Negative U flux
    Fm[0,:,:] = (rho*(2*lam3m*(gam-1.)+lam2m+lam1m))/(2*grid1.J*gam)
      
    Fm[1,:,:] = rho*(u*grid1.XI*(2.*lam4m*(gam-1.) + (lam2m + lam1m)) + (lam2m - lam1m)*a*grid1.xi_x)/(2.*grid1.J*grid1.XI*gam)

    Fm[2,:,:] = rho*(gam - 1.)/(grid1.xi_y*grid1.J*gam)*(lam3m*(Uc - grid1.xi_t) - lam4m*u*grid1.xi_x) + rho/(2.*grid1.XI**2*grid1.J*gam)*\
      (a*grid1.xi_y*grid1.XI*(lam2m - lam1m) + (v*grid1.xi_y**2 + v*grid1.xi_x**2)*(lam2m + lam1m) )

    t1 = -(lam1m*rho*(2*a*v*grid1.xi_y*grid1.XI*gam+2*a*u*grid1.xi_x*grid1.XI*gam-v**2*grid1.xi_y**2*gam-u**2*grid1.xi_y**2*gam-v**2*grid1.xi_x**2*gam-u**2*grid1.xi_x**2*gam-2*a*v*grid1.xi_y*grid1.XI-2*a*u*grid1.xi_x*grid1.XI+v**2*grid1.xi_y**2+u**2*grid1.xi_y**2-2*a**2*grid1.xi_y**2+v**2*grid1.xi_x**2+u**2*grid1.xi_x**2-2*a**2*grid1.xi_x**2))/(4*(grid1.xi_y**2+grid1.xi_x**2)*grid1.J*(gam-1)*gam)
    t2 = (lam2m*rho*(2*a*v*grid1.xi_y*grid1.XI*gam+2*a*u*grid1.xi_x*grid1.XI*gam+v**2*grid1.xi_y**2*gam+u**2*grid1.xi_y**2*gam+v**2*grid1.xi_x**2*gam+u**2*grid1.xi_x**2*gam-2*a*v*grid1.xi_y*grid1.XI-2*a*u*grid1.xi_x*grid1.XI-v**2*grid1.xi_y**2-u**2*grid1.xi_y**2+2*a**2*grid1.xi_y**2-v**2*grid1.xi_x**2-u**2*grid1.xi_x**2+2*a**2*grid1.xi_x**2))/(4*(grid1.xi_y**2+grid1.xi_x**2)*grid1.J*(gam-1)*gam)
    t3 = lam3m*rho*(gam -1.)*(grid1.xi_y*(v**2 - u**2) + 2.*u*v*grid1.xi_x)/(2.*grid1.xi_y*grid1.J*gam)
    t4 = lam4m*rho*u*(gam-1.)*(u*grid1.xi_y - v*grid1.xi_x)/(grid1.xi_y*grid1.J*gam) 
    Fm[3,:,:] = t1 + t2 + t3 + t4

    return Fp,Fm,Gp,Gm