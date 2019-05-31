from __future__ import division, print_function, absolute_import
import numpy as np
from . import phys

###
# Planck function (of frequency)
# From Pierrehumbert - PoPC
def Planck_nu(nu,T):
    u = phys.h*nu/(phys.k*T)
    u[u>500.] = 500. #To prevent overflow
    #return (2.*phys.h*nu**3/phys.c**2)/(np.exp(u)-1.) # [?]
    return 2.*(phys.k*T)**3/((phys.h*phys.c)**2) * u**3/(np.exp(u)-1.)

def Planck_lambda(lam,T):
    u = phys.h*phys.c/(lam*phys.k*T)
    u[u>500.] = 500. #To prevent overflow
    return (2.*(phys.k*T)**5/(phys.h**4*phys.c**3)) * u**5/(np.exp(u)-1.)


## Based on Ray's book.
## NOTE: Returns the radiance in W/m2/stearadian / [length]
## where [length] is the unit in which I'm measuring the spectral axis.
##
## Conversion factor 'k' is necessary because Planck_lambda (above) returns
## W/m2/stearadian/m. 'k' allows me to return W/m2/stearadian/cm.
def Planck_n(n,T,unit="cm^-1"):
    if unit=="m^-1":
        k = 1.
    elif unit=="cm^-1":
        n = n*100.
        k = 100.
    else:
        print( "(Planck_n) Error: unit not recognized!")
    return k*phys.c*Planck_nu(n*phys.c,T)



## -----------------------------------------------------------
## BELOW: compute Brightness temperatures

## Get Brightness temperature from inverting Planck function:
## [lam] = m                       # wavelength
## [I] = W/m^2 / stearadian / m    # intensity
def Tbrightness_lambda(lam,I):
    return phys.h*phys.c / (phys.k*lam* np.log(2.*phys.h*phys.c**2/(I * lam**5) + 1.) )

## [nu] = s^-1 = Hz                # frequency
## [I] = W/m^2 / stearadian / Hz   # intensity
def Tbrightness_nu(nu,I):
    return phys.h * nu / (phys.k * np.log(2.*phys.h*nu**3/(phys.c**2 * I) + 1.) )

## [n] = cm^-1
## [intensity] = W/m^2 / stearadian / cm^-1
def Tbrightness_n(n,I,unit="cm^-1"):
    if unit=="m^-1":
        k = 1.
    elif unit=="cm^-1":
        n = n*100.
        k = 100.
    else:
        print( "(Planck_n) Error: unit not recognized!")

    return Tbrightness_nu(n*phys.c, I/(k*phys.c))



## -----------------------------------------------------------
## BELOW: compute derivative with respect to temperature

def dPlanckdT_nu(nu,T):
    u = phys.h*nu/(phys.k*T)
    u[u>500.] = 500. #To prevent overflow
    return 2.*phys.h**2/(phys.k*phys.c**2) * (nu**4/T**2) * np.exp(u)/( (np.exp(u)-1.)**2 )

# note: assume 'n' is in cm^-1!
def dPlanckdT_n(n,T,unit="cm^-1"):
    if unit=="m^-1":
        k = 1.
    elif unit=="cm^-1":
        k = 100.
    else:
        print( "(Planck_n) Error: unit not recognized!")

    return (k*phys.c) * dPlanckdT_nu( phys.c*n*k ,T)
