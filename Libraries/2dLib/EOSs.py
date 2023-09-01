from pylab import *
from gas import *
def constant_mu(mu,T):
  mu[:,:] = mu_0
  return mu
def sutherland(mu,T):
  C1 = 1.458e-6
  S = 110.4
  mu[:,:] = C1*T[:,:]**(1.5)/(T + S)
  return mu
