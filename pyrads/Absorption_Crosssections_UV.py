from __future__ import division, print_function, absolute_import
import numpy as np
from . import phys
import os

'''
Implement UV scattering cross-sections.
Either data or fits.
'''

### -----------------------------------
### Global definitions here


### -----------------------------------
### Absorption crosssection for CO2
### based on eqn 6 in Venot+ (2013).

def get_kappaAbs_CO2(wavenr,T):
    lam = 1e7 / wavenr # cm-1 -> nm
    
    # fitting formula:
    a = lambda T: -42.26 + (9593.*1.44/T)
    b = lambda T: 4.82e-3 - 61.5*1.44/T
    Q = lambda T: (1.-np.exp(-667.4*1.44/T))**(-2) * (1.-np.exp(-1388.2*1.44/T))**(-1) * (1.-np.exp(-2449.1*1.44/T))**(-1)

    sigma = Q(T) * np.exp(a(T)+b(T)*lam)
    kappaAbs = sigma * phys.N_avogadro/phys.CO2.MolecularWeight * 1e-1  # cm^2/molec -> m^2/kg
    return kappaAbs
