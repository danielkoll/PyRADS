from __future__ import division, print_function, absolute_import
import numpy as np
from . import phys
import disort

### -----------------------------------
### HELPERS:

# ...

### -----------------------------------
###
# Use pyDISORT to compute scattering.
#
# INPUT:
#   p,     pressure at which we want to know fluxes, [0d,1d]
#   pres,  pressure grid, 1d
#   temp,  temperature grid, 1d
#   wave,  wavenr grid, 1d
#   Tau,   optical thickness, 2d [pressure x wavenr]
#   Omega, single scattering albedo, 2d [pressure x wavenr]
#
# OUTPUT:
#   
#
# USAGE:
#
#
# ---
# Notes:
# Pydisort throws errors when tau becomes very large.
# I think that's because fortran has fewer significant digits than numpy,
#   so num.truncation can lead to max(utau)>max(tau)...
#   --> truncate Tau,Omega input

def get_fluxes(pres_output,pres,temp,wave,Tau,Omega,cosz,Lstar,alpha_surf, \
                   output_angle = None,
                   N_streams=4,doThermal=False,Ts=None,sourceFn=None, \
                   UsrTau=None,UsrAng=None ):

    # ---
    # setup pyDISORT:
    params_disort = {}
    params_disort['verbose'] = 0          # 0 = suppress header

    # cosine of viewing zenith angle where to output the RT fields (Default: 1.)
    if output_angle is None:
        params_disort['umu'] = 1.
    else:
        params_disort['umu'] = output_angle    # (CAREFUL: might have to be monotonically increasing?)

    # Output returned at computational or user-defined tau,angles?
    #     (default: returns ouput at user-defined angles)
    #  -> note: I haven't gotten the computational option to work..
    if UsrTau is not None:
        params_disort['UsrTau'] = UsrTau
    if UsrAng is not None:
        params_disort['UsrAng'] = UsrAng

    # ---
    # if output is at user-defined levels:
    #    define those levels
    if pres_output=="toa":
        # (if we only want TOA output)
        p = np.array([0.])
        def get_utau(tau=None,p=None,pres=None):
            return np.array([0.])
    elif pres_output=="input":
        # (output on input grid)
        p = pres
        def get_utau(tau,p=None,pres=None):
            utau = tau
            return utau
    else:
        # (output at user-defined pressures)
        #     interpolate on log tau
        #     at toa: tau=0-> -inf -> NaN  with this method. manually fix
        #     also can get utau>max(tau) because of round-off?? manually fix -> still crashes!
        p = pres_output
        def get_utau(tau,p,pres):
            utau = 10**( np.interp( p,pres,np.log10(tau) ))
            utau[np.isnan(utau)] = 0.
            return utau

    # ---
    params_disort['iphas'] = 2          # phase function: 1=isotropic, 2=rayleigh
    params_disort['gg'] = 0.            # asymmetry parameter. =0 for Rayleigh
    params_disort['umu0'] = cosz
    params_disort['albedo'] = alpha_surf
    params_disort['Nstr'] = N_streams   # number of streams  -> ??

    params_disort['ibcnd'] = 0

    # using spectrally resolved stellar source function?
    # need same dims as 'wave'!
    if sourceFn is not None:
        source = sourceFn(wave)
    else:
        source = wave*0. + Lstar


    # show more/less output?
    params_disort['prnt'] = np.array([False, False, False, False, False])
    
    # include thermal emission?
    params_disort['plank'] = doThermal

    if doThermal:        
        # careful about mid-level vs interfaces..
        params_disort['temp'] = temp
        params_disort['btemp'] = Ts
        params_disort['temis'] = 0.  # emissivity of top boundary
        params_disort['ttemp'] = 0.  #Tstrat   # doesn't matter if temis=0?

        # HERE: assume fixed dnu in wave grid!!
        dwave = np.diff(wave)[0]

    # ---
    # iterate over spectral grid:
    #    cast p into array form, to deal with float input..

    Idirect_grid = np.zeros( (len(np.array(p)),len(wave)) )
    Iplus_grid = np.zeros( (len(np.array(p)),len(wave)) )
    Iminus_grid = np.zeros( (len(np.array(p)),len(wave)) )

    # --
    for nu in wave:        
        # ...
        if doThermal:
            params_disort['wvnmlo'] = nu - dwave/2.
            params_disort['wvnmhi'] = nu + dwave/2.
        
        # direct beam source:
        params_disort['fbeam'] = source[nu==wave] # (don't include cosz here!)

        # optical thickness:
        tau = np.squeeze(Tau[:,nu==wave])     # 1d

        # single scattering albedo:
        #    careful: tau is computed from kappa via trapz rule -> apply same to omega!
        #    i.e., int[f(x)dx] = sum [f(i+1)+f(i)]*0.5 * dx(i)
        omega = np.squeeze(Omega[:,nu==wave]) # 1d
        params_disort['w0'] = 0.5*(omega[0:-1] + omega[1:])
        
        #
        params_disort['utau'] = get_utau(tau,p,pres)

        #
        dtau = np.diff(tau)

        #
        Idirect, Iminus, Iplus, dfdt, uavg, uu, albmed, trnmed = disort.run(dtau, **params_disort)

        Idirect_grid[:,nu==wave] = Idirect[:,np.newaxis]  # cast into right shape..
        Iplus_grid[:,nu==wave] = Iplus[:,np.newaxis]
        Iminus_grid[:,nu==wave] = Iminus[:,np.newaxis]

    if doThermal:
        # THERMAL EMISSION has units=W/m2. transform to W/m2/cm-1:
        Idirect_grid = Idirect_grid/dwave
        Iplus_grid = Iplus_grid/dwave
        Iminus_grid = Iminus_grid/dwave

    return Idirect_grid,Iplus_grid,Iminus_grid
