#!/usr/bin/env python
"""
Test of the Python wrapper to the DISORT library

Module '_disort' is auto-generated with f2py (version:2).
"""

import numpy as np
from pylab import plt
import disort

##########################################################################################################

if __name__ == '__main__':

    # read Rayleigh optical thickness for 325.8 nm [Ozone fitting window]
    xy     = np.loadtxt('rayleigh_layer_opd.txt')
    dTau   = xy[::-1,1]  # FROM TOP TO BOTTOM
    z_atm  = xy[::-1,0]  # last altitude value missing, find in header
    z_atm  = np.insert(z_atm, 0, 120.)

    N_tau  = len(dTau)
    w0     = np.ones(N_tau)*1.
    iphas  = np.ones(N_tau,dtype='int')*2
    gg     = np.zeros(N_tau)

    cumTau = np.hstack([0.,dTau.cumsum()])
    uTau   = cumTau
    phi    = np.array([0.,60.,120.])
    umu0   = 1./np.sqrt(2.)
    fbeam  = 1./umu0  # Ensures fluxes to be normalized to one
    phi0   = 0.0
    albedo = 0.1
    umu    = np.array([-1.,-0.5,0.5,1.])
    prnt   = np.array([True, True, True, False, True])

    [rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed] =\
                                      disort.run(dTau, w0=w0, iphas=iphas, gg=gg,
                                                 umu0=umu0, phi0=phi0, albedo=albedo, fbeam=fbeam,
                                                 utau=uTau, umu=umu, phi=phi, prnt=prnt)

    rfltot = rfldir + rfldn
    print '\n# Energy conservation, R(TOA)+T(BOA)*(1-albedo) ~ 1:  %.3f' % (flup[0] + rfltot[-1]*(1.-albedo))

    plt.figure()
    plt.plot(rfltot, z_atm)

    plt.figure()
    plt.plot(flup, z_atm)
