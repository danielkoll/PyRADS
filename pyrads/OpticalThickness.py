'''
***********************************************************
This script computes absorption coefficients and
optical thicknesses.

Notes:
- A Voigt line shape is necessary at low pressures, where the Lorentz line shape produces kappa->infinity.
  This can matter e.g. if some GHGs exist in the stratosphere -- e.g., CO2 radiative forcing with non-isothermal stratospheres.
  Voigt lines are ~3-4x slower to compute than Lorentz approx.
***********************************************************
'''
from __future__ import division, print_function, absolute_import
import numpy as np
from .Absorption_Crosssections_HITRAN2016_lorentz import getKappa_HITRAN as getKappa_HITRAN_lorentz
from .Absorption_Crosssections_HITRAN2016_voigt import getKappa_HITRAN as getKappa_HITRAN_voigt
from . import Absorption_Continuum_MTCKD
from .Absorption_Continuum_MTCKD import get_H2OContinuum
from scipy.integrate import cumtrapz
from numba import jit

from .Thermodynamics import convert_molar_to_mass_ratio


# ---
##
def compute_tau_H2ON2(p,T,q,grid,params,RH=1.,use_numba=False,lineShape="voigt"):
    # setup fn wrapper:   all these fns need to accept exact same args
    if use_numba:
        from .Absorption_Crosssections_HITRAN2016_numba import getKappa_HITRAN_numba
        getKappa = getKappa_HITRAN_numba
    else:
        if lineShape=="lorentz":
            getKappa = getKappa_HITRAN_lorentz
        elif lineShape=="voigt":
            getKappa = getKappa_HITRAN_voigt
        else:
            print( "(OpticalThickness) error -- choice of lineshape not recognized!" )
    # ..
        
    kappa = np.zeros( (grid.Np,grid.Nn) )
    for pres,temp,q_H2O in zip(p,T,q):
        p_H2O = RH * params.esat(temp)  # ...

        print( "compute kappa at p,T = ",pres,temp)
        kappaH2O = getKappa(grid.n,grid.n0,grid.n1,grid.dn, \
                                "H2O",press=pres,press_self=p_H2O, \
                                temp=temp,broadening="mixed", lineWid=25., \
                                cutoff_option="fixed",remove_plinth=True)

        # add continuum:
        #    here I'm only using kappa from mtckd crosssection file,
        #    which doesn't include N2-N2 and similar continua.
        kappaH2O_cont = get_H2OContinuum(grid.n,temp,pres,p_H2O, \
                            exe_file=Absorption_Continuum_MTCKD.mtckd_exe_H2O_N2)

        kappa[ p==pres,: ] = kappaH2O*q_H2O + kappaH2O_cont*q_H2O  # save
    print( "done! \n")

    # Integrate to get optical thickness:
    p2d = np.tile( p,(grid.Nn,1) ).T
    tau = 1./(params.g*params.cosThetaBar) * cumtrapz( kappa,x=p2d,initial=0.,axis=0 )

    return tau


# ---
## Here: assume CO2 is a minor trace gas!
##     (I'm using params.R to compute R_mean, so ignoring mass contribution of CO2)
def compute_tau_H2ON2_CO2dilute(p,T,q,ppv_CO2,grid,params,RH=1.,use_numba=False,lineShape="voigt"):
    # setup fn wrapper:   all these fns need to accept exact same args
    if use_numba:
        from .Absorption_Crosssections_HITRAN2016_numba import getKappa_HITRAN_numba
        getKappa = getKappa_HITRAN_numba
    else:
        if lineShape=="lorentz":
            getKappa = getKappa_HITRAN_lorentz
        elif lineShape=="voigt":
            getKappa = getKappa_HITRAN_voigt
        else:
            print( "(OpticalThickness) error -- choice of lineshape not recognized!" )
    # ..        
    
    kappa = np.zeros( (grid.Np,grid.Nn) )
    kappa_h2o = np.zeros( (grid.Np,grid.Nn) )
    kappa_co2 = np.zeros( (grid.Np,grid.Nn) )
    for pres,temp,q_H2O in zip(p,T,q):
        p_H2O = RH * params.esat(temp)  # ...
        R_mean = q_H2O*params.Rv + (1.-q_H2O)*params.R
        q_CO2 = convert_molar_to_mass_ratio(ppv_CO2,params.R_CO2,R_mean)

        print( "compute kappa at p,T = ",pres,temp)
        kappaH2O = getKappa(grid.n,grid.n0,grid.n1,grid.dn, \
                                "H2O",press=pres,press_self=p_H2O, \
                                temp=temp,broadening="mixed", lineWid=25., \
                                cutoff_option="fixed",remove_plinth=True)

        kappaCO2 = getKappa(grid.n,grid.n0,grid.n1,grid.dn, \
                                "CO2",press=pres,press_self=0., \
                                temp=temp,broadening="air", lineWid=25., \
                                cutoff_option="fixed",remove_plinth=False)  # use of "air" broadening consistent with trace gas assumption

        # add continuum:
        #    here I'm only using kappa from mtckd crosssection file,
        #    which doesn't include N2-N2 and similar continua.
        kappaH2O_cont = get_H2OContinuum(grid.n,temp,pres,p_H2O, \
                            exe_file=Absorption_Continuum_MTCKD.mtckd_exe_H2O_N2)

        kappa[ p==pres,: ] = kappaH2O*q_H2O + kappaH2O_cont*q_H2O + kappaCO2*q_CO2  # save
        kappa_h2o[ p==pres,: ] = kappaH2O*q_H2O + kappaH2O_cont*q_H2O  # save
        kappa_co2[ p==pres,: ] = kappaCO2*q_CO2  # save
    print( "done! \n")

    # Integrate to get optical thickness:
    p2d = np.tile( p,(grid.Nn,1) ).T
    tau     = 1./(params.g*params.cosThetaBar) * cumtrapz( kappa,x=p2d,initial=0.,axis=0 )
    tau_h2o = 1./(params.g*params.cosThetaBar) * cumtrapz( kappa_h2o,x=p2d,initial=0.,axis=0 )
    tau_co2 = 1./(params.g*params.cosThetaBar) * cumtrapz( kappa_co2,x=p2d,initial=0.,axis=0 )
    return tau, tau_h2o, tau_co2


