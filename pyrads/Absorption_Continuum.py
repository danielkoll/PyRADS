from __future__ import division, print_function
import numpy as np
import phys

'''
Implement functions that return absorption cross-sections
for various modes of continuum absorption.
'''


###
# Implement fits from Pierrehumbert, PoPC.

# Here: Water vapor self-continuum; Eqn. (4.91)
# NOTE: 'wave' is a vector/matrix, T is a scalar
# OTHER NOTE: Ray's values are at 100mbar of H2O.
#           I need to rescale them to whatever the actual partial pressure of H2O is.
def get_H2OSelfContinuum(wave,T,pH2O,extrapolate=True):

    kappaC = np.zeros( wave.shape )  # kappa = 0 outside of region...

    mask = np.logical_and( wave>=500., wave<=1400. )
    kappaC[mask] = np.exp(12.167-0.050898*wave[mask] + 8.3207e-05*wave[mask]**2
                        -7.0748e-08*wave[mask]**3 + 2.3261e-11*wave[mask]**4)

    x = wave - 2500.
    mask = np.logical_and( wave>=2100., wave<=3000. )
    kappaC[mask] = np.exp(-6.0055 + -0.0021363*x[mask] + 6.4723e-07*x[mask]**2 +
                              -1.493e-08*x[mask]**3 + 2.5621e-11*x[mask]**4 +
                              7.328e-14*x[mask]**5)

    # just extrapolate the edges to ever higher/lower wavenumbers?
    if extrapolate:
        # extend left
        mask = wave < 500.
        kappa0 = np.exp(12.167-0.050898* 500. + 8.3207e-05* 500.**2
                        -7.0748e-08* 500.**3 + 2.3261e-11* 500.**4)
        kappaC[mask] = kappa0

        # extend right
        mask = wave>3000.
        kappa1 = np.exp(-6.0055 + -0.0021363*(3000.-2500.) + 6.4723e-07*(3000.-2500.)**2 +
                               -1.493e-08*(3000.-2500.)**3 + 2.5621e-11*(3000.-2500.)**4 +
                               7.328e-14*(3000.-2500.)**5)
        kappaC[mask] = kappa1

        # interpolate in between:
        mask = np.logical_and( wave>=1400., wave<=2100. )
        kappa2 = np.exp(12.167-0.050898* 1400. + 8.3207e-05* 1400.**2
                        -7.0748e-08* 1400.**3 + 2.3261e-11* 1400.**4)
        kappa3 = np.exp(-6.0055 + -0.0021363*(2100.-2500.) + 6.4723e-07*(2100.-2500.)**2 +
                               -1.493e-08*(2100.-2500.)**3 + 2.5621e-11*(2100.-2500.)**4 +
                               7.328e-14*(2100.-2500.)**5)
        #kappaC[mask] = kappa2 + (kappa3 - kappa2) * (wave[mask] - 1400.)/(2100. - 1400.)   # a) linear

        kappaC[mask] = 10**(np.log10(kappa2) + (np.log10(kappa3) - np.log10(kappa2)) * (wave[mask] - 1400.)/(2100. - 1400.))   # b) logarithmic


    # Temperature dependence & rescale for partial pressure:
    Tfac = (296./T)**4.25
    return kappaC *Tfac *(pH2O/1e4)
