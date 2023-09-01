from schemes import *

from libraries import getFiltMatrices
from libraries import bd_dX
from libraries import getCompactMat

## Differincing schemes
if (differencing_scheme == 'c2'):
  from libraries import d_dX as d_dX
  from libraries import gridMetrics as gridMetrics
if (differencing_scheme == 'compact_c4'):
  from libraries import getCompactMat
  from libraries import gridMetricsC as gridMetrics
  from libraries import d_dXC as d_dX
  
## Filtering Schemes
if (filter_scheme == '2'):
  from libraries import filt as filt
if (filter_scheme == '4'):
  from libraries import filt4 as filt
