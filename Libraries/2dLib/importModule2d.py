from schemes import *
from libraries2d import getFiltMatrices
from libraries2d import bd_dX
from libraries2d import getCompactMat
from libraries2d import gradedSpace
## Differincing schemes
if (differencing_scheme == 'KurgTadmor'):
    from KurgTadmor import kurgTadmorFlux as computeFlux
    from libraries2d import d_dx_central2 as d_dx
    from libraries2d import d_dy_central2 as d_dy
if (differencing_scheme == 'c2'):
  from libraries2d import d_dx_central2 as d_dx
  from libraries2d import d_dy_central2 as d_dy
  from central import centralFlux as computeFlux

if (differencing_scheme == 'c4'):
  from libraries2d import d_dx_compact4 as d_dx
  from libraries2d import d_dy_compact4 as d_dy
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


from Jacobians import inviscidJacobians

exists = 1
try:
  ddt_scheme
except NameError:
  exists = 0
if (exists != 0):
  if (ddt_scheme == 'rk4'):
    from timeIntegrationSchemes import rk4 as advanceTime
  if (ddt_scheme == 'implicitEuler'):
    from timeIntegrationSchemes import implicitEuler as advanceTime
else:
  print('ddt_scheme not found. Using rk4')
  from timeIntegrationSchemes import rk4 as advanceTime


exists = 1
try:
  gridMetrics
except NameError:
  exists = 0
if (exists != 0):
  if (gridMetrics == 'c4'):
    from libraries2d import gridMetrics_compact4 as gridMetrics
  if (gridMetrics == 'c2'):
    from libraries2d import gridMetrics_central2 as gridMetrics
else:
  from libraries2d import gridMetrics_central2 as gridMetrics
exists = 1
try:
  viscous
except NameError:
  exists = 0
if (exists != 0):
  if (viscous == 'on'):
    from viscousFluxes import viscFlux
    if (muEOS == 'sutherland'):
      from EOSs import sutherland as muEOS
    if (muEOS == 'constant_mu'):
      from EOSs import constant_mu as muEOS









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
