from schemes import *

from libraries2d import getFiltMatrices
from libraries2d import bd_dX
from libraries2d import getCompactMat
from libraries2d import gradedSpace
from libraries2d import gridMetrics as gridMetrics

## Differincing schemes
if (differencing_scheme == 'KurgTadmor'):
    from KurgTadmor import kurgTadmorFlux as computeFlux
if (differencing_scheme == 'c2'):
  from libraries2d import d_dx as d_dx
  from libraries2d import d_dy as d_dy
  from central import centralFlux as computeFlux

if (differencing_scheme == 'StegerWarming1'):
  from libraries2d import d_dxF1 as d_dxF
  from libraries2d import d_dxB1 as d_dxB
  from libraries2d import d_dyF1 as d_dyF
  from libraries2d import d_dyB1 as d_dyB

if (differencing_scheme == 'StegerWarming2'):
  from libraries2d import d_dxF2 as d_dxF
  from libraries2d import d_dxB2 as d_dxB
  from libraries2d import d_dyF2 as d_dyF
  from libraries2d import d_dyB2 as d_dyB

if (differencing_scheme == 'StegerWarming3'):
  from libraries2d import d_dxF3 as d_dxF
  from libraries2d import d_dxB3 as d_dxB
  from libraries2d import d_dyF3 as d_dyF
  from libraries2d import d_dyB3 as d_dyB

#if (differencing_scheme == 'compact_c4'):
#  from libraries import getCompactMat
#  from libraries import gridMetricsC as gridMetrics
#  from libraries import d_dXC as d_dX

from KurgTadmor import *  
## Filtering Schemes
if (filter_scheme == '2'):
  from libraries2d import filt2 as filt
if (filter_scheme == '4'):
  from libraries2d import filt4 as filt
