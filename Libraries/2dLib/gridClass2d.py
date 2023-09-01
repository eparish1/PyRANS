from numpy import *
from numpy.linalg import*
from mesh import x,y
from importModule2d import *
class grid():
  nx,ny = shape(x)
  A = gridMetrics(x,y)
  from mesh import xperd as xperd
  from mesh import yperd as yperd
  Afx,Afy = getFiltMatrices(x,xperd,yperd) #get inverse matrices for filtering
  from mesh import yperd_grid
  from mesh import xperd_grid
  J,xi_x,eta_x,xi_y,eta_y,XI,ETA,xi_t,eta_t = gridMetrics(x,y)
  Cx,Cy = getCompactMat(nx,ny,xperd,yperd)
#  xi_y += 1e-9
#  eta_y += 1e-9
#  xi_x += 1e-9
#  eta_x += 1e-9
  
  #Cx,Cy = getCompactMat(x)

grid1 = grid()
