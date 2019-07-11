#!/usr/bin/env python
"""
Python wrapper to the DISORT library

Module '_disort' is auto-generated with f2py (version:2).
"""
from __version__ import __version__

import _disort

##########################################################################################################
#
# functions
#
##########################################################################################################

def run(dTau, w0=1., iphas=2, gg=0.85,
        umu0=1., phi0=0., albedo=0.1, fbeam=1.0,
        UsrTau=True, utau=0., UsrAng=True, umu=1.,
        phi=0., Nstr=32, maxmom=299, lamber=True,
        onlyFl=False, accur=0., plank=False,
        temp=300., wvnmlo=999., wvnmhi=1000.,
        ibcnd=0, fisot=0., btemp=300., ttemp=300.,
        temis=1., prnt=[False]*5, verbose=1):

    """

    performs radiative transfer simulations by means of the DISORT RT solver

    Parameters
    ----------
    dTau : array
        optical thickness in atmospheric layers
    w0 : float, array
        single scattering albedo (Default: 1.)
    iphas : int, array
        scattering phase function type (Default 2).
        1 : Isotropic
        2 : Rayleigh
        3 : Henyey-Greenstein with asymmetry factor GG
        4 : Haze L as specified by Garcia/Siewert
        5 : Cloud C.1 as specified by Garcia/Siewert
    gg : float, array
        scattering asymmetry parameter (Default: 0.85)
    umu0 : float
        cosine of solar zenith angle (Default: 1.)
    phi0 : float
        solar azimuth angle (Default: 0.)
    albedo : float
        surface albedo (Default: 0.1)
    fbeam : float
        solar irradiance (Default: 1.)
    utau : float, array
        optical thickness where to output the RT fields (Default: 0.)
    umu : float, array
        cosine of viewing zenith angle where to output the RT fields (Default: 1.)
    phi : float, array
        viewing azimuth angle where to output the RT fields (Default: 0.)
    maxmom : int
        Max. number of Legendre coefficients. (Default: 299)
    Nstr : int
        Number of computational polar angles to be used
        (= number of 'streams')  ( should be even and .GE. 2 ).
        (Default: 32)
    temp : float, array
        LEV = 0 to NLYR, Temperatures (K) of levels.
        (Note that temperature is specified at LEVELS
        rather than for layers.)  Be sure to put top level
        temperature in TEMPER(0), not TEMPER(1).  Top and
        bottom level values do not need to agree with top and
        bottom boundary temperatures (i.e. temperature
        discontinuities are allowed).  (Default: 300.)
    wvnmlo, wvnmhi : float
        Wavenumbers (inv cm) of spectral interval of interest
        ( used only for calculating Planck function ).
        Needed only if PLANK is TRUE, or in multiple runs, if
        LAMBER is FALSE and BDREF depends on spectral interval.
        (Default: wvnmlo=999., wvnmhi=1000.)
    UsrTau : logical
        = FALSE, Radiant quantities are to be returned
                     at boundary of every computational layer.
        = TRUE,  Radiant quantities are to be returned
                     at user-specified optical depths
        (Default: True)
    UsrAng : logical
        = FALSE, Radiant quantities are to be returned
                     at computational polar angles.
        = TRUE,  Radiant quantities are to be returned
                     at user-specified polar angles.
        (Default: True)
    ibcnd : int
        = 0, General case: boundary conditions any combination of:
             * beam illumination from the top ( see FBEAM )
             * isotropic illumination from the top ( see FISOT )
             * thermal emission from the top ( see TEMIS, TTEMP )
             * internal thermal emission sources ( see TEMPER )
             * reflection at the bottom ( see LAMBER, ALBEDO, BDREF )
             * thermal emission from the bottom ( see BTEMP )
        = 1, Return only albedo and transmissivity of the entire
             medium vs. incident beam angle; see S2 for details.
        (Default: 0)
    fisot : float
        Intensity of top-boundary isotropic illumination.
        [same units as PLKAVG (default W/sq m) if thermal
        sources active, otherwise arbitrary units].
        Corresponding incident flux is pi (3.14159...)
        times FISOT.
        (Default: 0.)
    lamber : bool
        = TRUE, isotropically reflecting bottom boundary.
        = FALSE, bidirectionally reflecting bottom boundary.
        (Default: True)
    btemp : float
        Temperature of bottom boundary (K)  (bottom emissivity
        is calculated from ALBEDO or function BDREF, so it need
        not be specified). Needed only if PLANK is TRUE.
        (Default: 300.)
    ttemp : float
        Temperature of top boundary (K).
        Needed only if PLANK is TRUE.
        (Default: 300.)
    temis : float
        Emissivity of top boundary.
        Needed only if PLANK is TRUE.
        (Default: 1.)
    plank : bool
        = TRUE,  include thermal emission
        = FALSE, ignore all thermal emission (saves computer time)
        (Default: False)
    onlyFl : bool
        = TRUE, return fluxes, flux divergences, and mean
                intensities.
        = FALSE, return fluxes, flux divergences, mean
                 intensities, AND intensities.
        (Default: False)
    accur : float
        Convergence criterion for azimuthal (Fourier cosine)
        series.  Will stop when the following occurs twice:
        largest term being added is less than ACCUR times
        total series sum.  (Twice because there are cases where
        terms are anomalously small but azimuthal series has
        not converged.)  Should be between 0 and 0.01 to avoid
        risk of serious non-convergence.  Has no effect on
        problems lacking a beam source, since azimuthal series
        has only one term in that case.
        (Default: 0.)
    PRNT : array(dtype=bool)
        Array of LOGICAL print flags causing the following prints
           ===      ===============================================
           L        quantities printed
           ===      ===============================================
           1        input variables (except PMOM)
           2        fluxes
           3        intensities at user levels and angles
           4        planar transmissivity and planar albedo
                    as a function solar zenith angle ( IBCND = 1 )
           5        phase function moments PMOM for each layer
                    ( only if PRNT(1) = TRUE, and only for layers
                    with scattering )
           ===      ===============================================
        (Default: array([False False False False False]))

    Returns
    -------
    ds_fields : list of arrays
        [rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed]

        rfldir : Downward Direct
        rfldn  : Downward Diffuse
        flup   : Upward Diffuse
        dfdt   : d(Net Flux) / d(Op Dep)
        uu     : Intensity
        uavg   : Mean intensity (including the direct beam)
                 (Not corrected for delta-M-scaling effects)
        albmed : Albedo of the medium as a function of incident
                 beam angle cosine UMU(IU)  (IBCND = 1 case only)
        trnmed : Transmissivity of the medium as a function of incident
                 beam angle cosine UMU(IU)  (IBCND = 1 case only)

    Examples
    --------
    >>> import disort
    >>> D_dir, D_diff, U_up, dFdt, I = disort.run(dTau, ssalb, iphas='Rayleigh')
    """

    import numpy as np

    if not hasattr(w0,'__iter__'):
        w0 = w0 * np.ones_like(dTau)
    if not hasattr(iphas,'__iter__'):
        iphas = iphas * np.ones_like(dTau)
    if not hasattr(gg,'__iter__'):
        gg = gg * np.ones_like(dTau)
    if not hasattr(temp,'__iter__'):
        temp = temp * np.ones(len(dTau)+1)
    if not hasattr(utau,'__iter__'):
        utau = np.array([utau])
    if not hasattr(umu,'__iter__'):
        umu = np.array([umu])
    if not hasattr(phi,'__iter__'):
        phi = np.array([phi])
    if not hasattr(prnt,'__iter__'):
        prnt = prnt * np.ones(5, dtype='bool')

    if verbose > 0:
        print ' **************************************************************'+\
            '**************************************'
        print ' DISORT: Python wrapper to the DISORT radiative transfer solver'
        print ' **************************************************************'+\
            '**************************************'

    rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed = \
           _disort.run(dTau, w0, maxmom, temp, iphas, gg,
                       wvnmlo, wvnmhi, UsrTau, utau, Nstr,
                       UsrAng, umu, phi, ibcnd, fbeam,
                       umu0, phi0, fisot, lamber, albedo, btemp,
                       ttemp, temis, plank, onlyFl, accur, prnt)

    return rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed
