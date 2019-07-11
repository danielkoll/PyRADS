'''
***********************************************************
This script setups the vertical and spectral grids.
***********************************************************
'''
from __future__ import division, print_function, absolute_import
import numpy as np
from . import VerticalStructure
from . import Thermodynamics

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



## AGAIN, but now allow vertical resolution to vary with surface pressure!
##     -> keep adding more points as pressure keeps rising to make sure resolution near TOA is fixed
def make_grid_fixedTop( Ts,Tstrat,wavenr_min,wavenr_max,dwavenr,params, \
                   ptop=0.1, Npress_per_decade=6., RH=1., adiabat="full" ):

    grid = Dummy()
    grid.Ts = Ts
    grid.ps = params.ps_dry + RH * params.esat(Ts)

    # setup constants
    ptop = min( ptop,grid.ps*0.0002 )       # make sure ptop < ps
    grid.Np = int(round( (np.log10(grid.ps)-np.log10(ptop))*Npress_per_decade))   # 
    grid.n0, grid.n1, grid.dn = wavenr_min, wavenr_max, dwavenr

    # setup arrays
    grid.n = np.arange( grid.n0, grid.n1, grid.dn )  # 1d
    grid.wave = np.tile( grid.n,(grid.Np,1) )         # 2d
    grid.Nn = len(grid.n)
    
    # start filling arrays: compute p,T,q profiles
    #  -> careful about starting pressure!
    #
    #  -> DISORT also needs temps at mid-points
    #     say pmid[i] = p[i] + (p[i+1]-p[i])/2.
    #
    # beware:
    #  roundoff can cause max(p) > ps. Make sure that doesn't happen.

    # ---
    # OLD: log spaced
#     grid.p = np.logspace( np.log10(ptop),np.log10(grid.ps),endpoint=True,num=grid.Np ) # 1d [Pa]

#     if max(grid.p) > grid.ps:
#         #grid.p = grid.p * grid.ps/max(grid.p)   # fixes roundoff? (OLD)
#         grid.p = grid.p * 0.99999 * grid.ps/max(grid.p)   # fixes roundoff? (NEW)

    # ---
    # NEW v2: log spaced; 2x resolution in last decade; 1.5x in the 2nd last
    Np1 = int(Npress_per_decade * 1.5)
    Np2 = int(Npress_per_decade * 2.)
    Np0 = grid.Np - (Np1+Np2)

    p0 = np.logspace( np.log10(ptop),np.log10(grid.ps)-2,endpoint=False,num=Np0 )
    p1 = np.logspace( np.log10(grid.ps)-2,np.log10(grid.ps)-1,endpoint=False,num=Np1 )
    p2 = np.logspace( np.log10(grid.ps)-1,np.log10(grid.ps),endpoint=True,num=Np2 )
    grid.p = np.concatenate( (p0,p1,p2) ) * 0.99999    # 1d [Pa]
    # ---

    #grid.pmid = np.concatenate( (np.array([0.]),grid.p[0:-1] + np.diff(grid.p)/2.) )
    grid.pmid = grid.p[0:-1] + np.diff(grid.p)/2.      # NOTE: len(pmid) = len(p)-1!

    if adiabat=="full":
        grid.T,grid.q = VerticalStructure.get_Tq_moist(grid.p,Ts,grid.ps,params,Tstrat=Tstrat,RH=RH)
        grid.Tmid,tmp = VerticalStructure.get_Tq_moist(grid.pmid,Ts,grid.ps,params,Tstrat=Tstrat,RH=RH)
    elif adiabat=="dry":
        grid.T = VerticalStructure.get_TofP_DryAdiabat(grid.p,Ts,grid.ps,params)
        grid.q = Thermodynamics.get_q(grid.T,grid.p,params,RH=RH)
        grid.Tmid = VerticalStructure.get_TofP_DryAdiabat(grid.pmid,Ts,grid.ps,params)
    elif adiabat=="steam":
        grid.T = VerticalStructure.get_TofP_SingleCompAdiabat(grid.p,params)
        grid.q = grid.T * 0. + 1.
        grid.Tmid = VerticalStructure.get_TofP_SingleCompAdiabat(grid.pmid,params)
    else:
        print( "(setupgrids) error: adiabat choice=",adiabat, " -> not recognized!  *****" )

    return grid
