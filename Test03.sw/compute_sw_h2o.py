from __future__ import division, print_function
import numpy as np
import sys,os

sys.path.append("..")
import pyrads

from scipy.integrate import trapz,simps,cumtrapz
import time

### -----------------------------------
### Helpers
class Dummy:
    pass


### -----------------------------------

# ---
## setup thermodynamic parameters
params = Dummy()

params.Rv = pyrads.phys.H2O.R # moist component
params.cpv = pyrads.phys.H2O.cp
params.Lvap = pyrads.phys.H2O.L_vaporization_TriplePoint
params.satvap_T0 = pyrads.phys.H2O.TriplePointT
params.satvap_e0 = pyrads.phys.H2O.TriplePointP
params.esat = lambda T: pyrads.Thermodynamics.get_satvps(T,params.satvap_T0,params.satvap_e0,params.Rv,params.Lvap)

params.R = pyrads.phys.air.R  # dry component
params.cp = pyrads.phys.air.cp
params.ps_dry = 1e5           # surface pressure of dry component

params.g = 9.8             # surface gravity
params.cosThetaBar = 3./5. # average zenith angle used in 2stream eqns
params.RH = 1.             # relative humidity


## scattering parameters:
params.Lstar = 1365.         # ...
params.alpha_surf = 0.12     # surface albedo

params.Tstar = 5800.
params.sourceFn = lambda wave: np.pi * pyrads.Planck.Planck_n(wave,params.Tstar) * \
    (params.Lstar)/(pyrads.phys.sigma* params.Tstar**4)   # units = W/m2/cm^-1

params.zenith_list = np.array([7.5,22.5,37.5,52.5,67.5,82.5])   # angles for scattering?

## other parameters:
params.Tstrat = 150.

## numerical parameters:
params.N_press = 15   # MAKE SURE THIS IS AN 'INT'! -- chosen low res for testing only!
params.RH = 1.



### -----------------------------------
### Main

def main(wavenr_min,wavenr_max,dwavenr,Ts,OUTDIR='./'):
    print( "Ts, wavenr_min, wavenr_max, dwave = ", Ts,wavenr_min,wavenr_max,dwavenr )


    # ----
    # setup grid:
    grid = pyrads.SetupGrids.make_grid( Ts,params.Tstrat,params.N_press, \
                          wavenr_min,wavenr_max,dwavenr, \
                          params, pTOA_decade=-4, RH=params.RH )


    # ----
    # compute optical thickness:
    print( "Compute tau, omega... " )
    grid.tau,grid.omega = pyrads.OpticalThickness.compute_tau_omega_H2ON2(grid.p,grid.T,grid.q,grid,params, RH=params.RH )
    print( "Done " )

    # ---
    # compute SW fluxes
    # run pyDISORT over all zenith angles  -> this can get long!
    
    grid.output_angles = []
    
    start01 = time.time()  # TESTING
    for zenith in params.zenith_list:
        print( "\n---" )
        print( "\nZenith angle = ", zenith )

        cosz = np.cos( zenith * np.pi/180. )

        #
        out = Dummy()
        out.zenith = zenith
        out.cosz = cosz
        
        # run pyDISORT
        SWdirect_spec,SWplus_spec,SWminus_spec = pyrads.Get_Fluxes_pyDISORT.get_fluxes( \
                    grid.p*0.9999, \
                    grid.p,grid.T, \
                    grid.n, grid.tau, grid.omega, \
                    cosz, params.Lstar,params.alpha_surf, \
                    N_streams=4,doThermal=False, \
                    sourceFn=params.sourceFn )

        # convert to total fluxes
        #   output has same units as source fn
        #   integrate to get W/m2
        out.SWdirect = np.trapz( SWdirect_spec,x=grid.n,axis=-1 )
        out.SWplus = np.trapz( SWplus_spec,x=grid.n,axis=-1 )
        out.SWminus = np.trapz( SWminus_spec,x=grid.n,axis=-1 )

        #
        grid.output_angles.append( out )

    print( "Done! \n" )
    end01 = time.time()  # TESTING
    print( "***\nTiming for all solar angles: ", end01-start01, "\n***" )  # TESTING

    # ---
    # save to file
    OUTDIR = OUTDIR + "/"    # (make sure to add a '/' at end!!)
    print( "Saving output to ",OUTDIR )
    if not os.path.isdir( OUTDIR ):
        os.makedirs( OUTDIR )

    #
    zenith_list = [out.zenith for out in grid.output_angles]
    SWdirect_list = [out.SWdirect for out in grid.output_angles]
    SWminus_list = [out.SWminus for out in grid.output_angles]
    SWplus_list = [out.SWplus for out in grid.output_angles]

    pyrads.Write_Data.save_profile(OUTDIR,grid.n0,grid.n1,grid.dn, \
                                grid.p,grid.T,grid.q, \
                                zenith_list,SWdirect_list,SWminus_list,SWplus_list )



## Finally, allow main() to be called from command line:
##    USAGE:
##    python compute_sw_h2o.py 1000. 1100. 0.01 300. dir_output

if __name__ == "__main__":
    n0 = float(sys.argv[1])
    n1 = float(sys.argv[2])
    dn = float(sys.argv[3])
    Ts = float(sys.argv[4])
    outdir = sys.argv[5]
    main(n0,n1,dn,Ts,OUTDIR=outdir)


