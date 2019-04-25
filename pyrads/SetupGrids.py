'''
***********************************************************
This script setups the vertical and spectral grids.
***********************************************************
'''
from __future__ import division, print_function
import numpy as np
import VerticalStructure
import Thermodynamics

# ---
## Helpers
class Dummy:
    pass

# ---
## Setup the grid (spectral, pressure) here.
## FOR A FIXED SURFACE TEMPERATURE.
## Put into a single object.

def make_grid( Ts,Tstrat,N_press,wavenr_min,wavenr_max,dwavenr,params, \
                   pTOA_decade=-4, RH=1., adiabat="full" ):

    grid = Dummy()
    grid.Ts = Ts
    grid.ps = params.ps_dry + RH * params.esat(Ts)

    # setup constants
    grid.Np = N_press
    grid.n0, grid.n1, grid.dn = wavenr_min, wavenr_max, dwavenr

    # setup arrays
    grid.n = np.arange( grid.n0, grid.n1, grid.dn )  # 1d
    grid.wave = np.tile( grid.n,(grid.Np,1) )         # 2d
    grid.Nn = len(grid.n)

    # start filling arrays: compute p,T,q profiles
    #  -> careful about starting pressure!
    #  -> note: stratosphere only used for a full moist adiabat!

    grid.p = np.logspace( pTOA_decade,0,endpoint=True,num=grid.Np ) * grid.ps  # 1d [Pa]
    if adiabat=="full":
        grid.T,grid.q = VerticalStructure.get_Tq_moist(grid.p,Ts,grid.ps,params,Tstrat=Tstrat,RH=RH)
    elif adiabat=="dry":
        grid.T = VerticalStructure.get_TofP_DryAdiabat(grid.p,Ts,grid.ps,params)
        grid.q = Thermodynamics.get_q(grid.T,grid.p,params,RH=RH)
    elif adiabat=="steam":
        grid.T = VerticalStructure.get_TofP_SingleCompAdiabat(grid.p,params)
        grid.q = grid.T * 0. + 1.
    else:
        print( "(setupgrids) error: adiabat choice=",adiabat, " -> not recognized!  *****")

    return grid
