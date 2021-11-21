from __future__ import division, print_function, absolute_import
import numpy as np

from .Planck import Planck_n

# Here: choose integrator
from scipy.integrate import trapz as numint
##from scipy.integrate import simps as numint


####
# Transmission function between level i and j
def trans(i,j,tau):
    return np.exp( -np.abs(tau[i,:] - tau[j,:]) )

# Upward flux at level i (returns all spectral intervals)
# Implements eqn.(4.12) in PoPC.
#
# INPUT:
# 'data' object with following fields attached
# B,tau - 2D matrices, dimensions pressure x wavenumber
# B_surf - 1D matrix, 1 x wavenumber
#
# B,B_surf is emission, pi*Planck(T,nu)
# tau is the optical depth
#
# OUTPUT:
# upwards thermal flux, 1 x wavenumber
#
# Notes:
#  A[slice(i,None)] returns A[i:]
#  trans(i,slice(i,None)) is the transmission function between i-th level down to surface level
#  B[i:,:] is the emission between i-th level down to surface

def Fplus(i,data):
    return data.B[i,:] + (data.B_surf-data.B[-1,:])*trans(i,-1,data.tau) + \
      numint( trans(i,slice(i,None),data.tau ),x=data.B[i:,:],axis=0 )


### Implement eqn.(4.11) instead.
### Note the minus sign, which comes from Pierrehumbert's sign convention for tau.
def Fplus_alternative(i,data):
    return data.B_surf*trans(i,-1,data.tau) - numint( data.B[i:,:],x=trans(i,slice(i,None),data.tau), axis=0)


### Compute exp(-tau) weighted vertical integrals.
### For feedback calculations, I need to compute atmospheric integrals like this:
###     int_{0}^{surf}  f(tau) exp(-tau) dtau
### where f(tau) = pi B(T(tau)) for OLR
###       f(tau) = pi dB/dT(T(tau)) for the planck feedback etc.

def integrate_over_tau(i,data,integrand):
    return -1.* numint( integrand[i:,:],x=trans(i,slice(i,None),data.tau), axis=0)
