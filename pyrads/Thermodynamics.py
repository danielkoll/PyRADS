'''
Thermodynamic helper functions.
'''
from __future__ import division, print_function
import numpy as np


# Saturation vapor pressure from the Clausius-Clapeyron relation.
# --> assumes L is constant with temperature!
def get_satvps(T,T0,e0,Rv,Lv):
    return e0*np.exp(-(Lv/Rv)*(1./T - 1./T0))


# ---
# Input: T and *total* pressure. (Partial dry pressure is inferred.)
# Output: saturation specific humidity
def get_qsat(T,p_tot,params):
    eps = params.R/params.Rv
    p_dry = p_tot - params.esat(T)
    qsat = eps* params.esat(T)/( eps* params.esat(T) + p_dry )
    return qsat


# ---
# Input: T and *total* pressure. (Partial dry pressure is inferred.)
# Output: specific humidity.
#
# --> Valid for any value of RH! If RH=1, reduces to get_qsat.
def get_q(T,p_tot,params,RH=1.):
    eps = params.R/params.Rv
    p_h2o = RH*params.esat(T)
    p_dry = p_tot - p_h2o
    q = eps*p_h2o/(p_dry + eps*p_h2o)
    return q


# ---
# Input: T and *total* pressure.
# Output: mass mixing ratio r := rho_H2O / rho_dry = m_H2O*p_H2O / (m_dry * p_dry)
#
# NOTE: this is NOT the volume/number mixing ratio !
#       for that simply use n_H2O/n_dry = p_H2O/p_dry
def get_rsat(T,p_tot,params):
    eps = params.R/params.Rv
    e_sat = params.esat(T)
    p_dry = p_tot - e_sat
    r_sat = eps * e_sat / p_dry
    return r_sat





# ---
# Input: abundance of gas, by volume and/or number of molecules
#         molar concentration   eta := n_i / (n_i + n_rest)
# Output: abundance of gas, by mass
#         mass mixing ratio   q := rho_i / (rho_i + rho_rest)
#
# where rho_i is (mass) density, and n_i is number density.
#
# Formula is (see, e.g., PoPC p.86-87):
#    q_i = eta_i * <R>/R_i,    where <R>  is the mass-weighted average gas constant
#           <R> = q_a*R_a + q_b*R_b + ...
#  If I'm only given molar concentrations, use instead
#           <R> = 1./( eta_a/R_a + eta_b/R_b + ...)
#
# USES, e.g.:
#   - convert 300 ppmv of CO2 into a mass mixing ratio.

def convert_molar_to_mass_ratio(molar_i,R_i,R_air):
    molar_air = 1. - molar_i
    R_mean = 1./(molar_i/R_i + molar_air/R_air)
    q_i = molar_i * R_mean/R_i
    return q_i
