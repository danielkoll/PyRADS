'''
***********************************************************
This script provides a number of functions to compute
vertical temperature,humidity structure.
***********************************************************
'''
from __future__ import division, print_function, absolute_import

import numpy as np
import scipy
from scipy.optimize import fsolve  # [!]
from .Thermodynamics import get_qsat,get_q,get_satvps


# Returns T-p, q-p profiles, assuming a moist adiabat.
# Allow for a isothermal stratosphere with q equal to
# its value at the tropopause.
#
# HERE: T from moist adiabat, q scales with qsat & RH.
# (~T is set by a convecting region, but mean q is lower due to additional processes)
#
# NOTE: q = RH*qsat is only valid if water vapor is dilute!
#
# USAGE:
# T_vector,q_vector = get_Tq_moist( p_vector,Ts,ps,Tstrat=200. )

def get_Tq_moist(p,Ts,ps,params,RH=1.,Tstrat=None):

    # careful: I might need to flip p-array for my solver to work
    # --> need p to increase monotonically!
    T_adiabat = get_TofP_MoistAdiabat(p[::-1],Ts,ps,params)[::-1]

    q_adiabat = get_q(T_adiabat,p,params,RH=RH)

    if Tstrat is not None and np.any(T_adiabat < Tstrat):
        # caution: only implements stratosphere if T<Tstrat anywhere
        # this can fail at crappy vertical resolution!

        mask = T_adiabat < Tstrat
        T = T_adiabat
        T[mask] = Tstrat

        # pick highest pressure level in strat, assume q is uniform at that value
        q = q_adiabat
        q[mask] = q_adiabat[p==(p[mask].max())]
    else:
        T = T_adiabat
        q = q_adiabat

    return T,q


# Numerical moist adiabat, T=T(p),
# following Ding & Pierrehumbert (2016).
# which allows for fully non-dilute atmospheres.
# p is the total pressure, pa is the partial pressure of dry component,
# with p = pa + pc.
#
# CAREFUL: is 'ps_in' the total surface pressure, or partial dry pressure?
# --> for now I'm assuming it's the total pressure!

def get_TofP_MoistAdiabat(p,Ts,ps_in,params):
    esat = lambda T: get_satvps(T,params.satvap_T0,params.satvap_e0,params.Rv,params.Lvap)
    rsat = lambda T,pa: params.R/params.Rv * esat(T)/pa

    #ps = ps_in + esat(Ts)  # here: assume ps_in is dry pressure only
    ps = ps_in

    # (integrate dlogT/dlog(p/ps)): based on eqn. 12 in Ding & Pierrehumbert (2016)
    def slope(logT,logpps):
        T = np.exp(logT)
        p = np.exp(logpps) * ps

        pc = esat(T)
        pa = p - pc

        alpha = params.Lvap/(params.R*T)
        beta = params.Lvap/(params.Rv*T)
        gamma = params.Lvap/(params.cp*T)

        return 1./( pc/p*beta + \
                        pa/p*params.cp/params.R * \
                        (1.+ (params.cpv/params.cp + (beta-1.)*gamma)*rsat(T,pa)) / (1.+alpha*rsat(T,pa)) )

    # First element of log(p/ps) vector needs to be the initial condition.
    #      log(Ts)<->log(p/ps=1)=0.
    # Rest contains p levels that we want to know T for.
    if hasattr(p,"__len__"):
        pps_vector = np.append([1.],p/ps)
    else:
        pps_vector = [1.,p/ps]

    logInt = scipy.integrate.odeint( slope,np.log(Ts),np.log(pps_vector) )
    logInt = np.squeeze( logInt )

    if hasattr(p,"__len__"):
        return np.exp(logInt[1:])  # don't return initial value...
    else:
        return np.exp(logInt[-1])  # don't return initial value...


# analytical dry adiabat, T=T(p)
def get_TofP_DryAdiabat(p,Ts,ps,params):
    return Ts*( p/ps )**(params.R/params.cp)

# analytical single-component moist adiabat, T=T(p)
def get_TofP_SingleCompAdiabat(p,params):
    return params.satvap_T0 / (1. - params.Rv*params.satvap_T0/params.Lvap * np.log(p/params.satvap_e0) )
