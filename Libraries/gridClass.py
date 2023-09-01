from numpy import *
from mesh import x,y,z
from mesh import xperd as xperd1
from mesh import yperd as yperd1
from importModule import *
class grid():
  nx,ny,nz = shape(x)
  A = gridMetrics(x,y,z)
  Afx,Afy,Afz = getFiltMatrices(x) #get inverse matrices for filtering
  xperd = xperd1
  yperd = yperd1
  Cx,Cy,Cz = getCompactMat(x)

grid1 = grid()
