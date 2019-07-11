from __future__ import division, print_function
import numpy as np
import sys,os

sys.path.append("..")
import pyrads

from scipy.integrate import trapz,simps,cumtrapz

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

# ---
## setup resolution (vertical,spectral)

#N_press = 60       #
N_press = 15       #
wavenr_min = 0.1   # [cm^-1]
wavenr_max = 3500. #
#dwavenr = 0.01     #
dwavenr = 0.1     #

Tstrat = 150.      # stratospheric temperature

# ---
## setup range of temperatures, and if/where output is saved to:

#Ts_grid = np.arange(170.,370.1,10.)
Ts_grid = np.arange(170.,370.1,20.)
filename = 'output.compute_olr_h2o.01.100RH.txt'

saveOutput = True  # Save the output/plots? [Yes/No]
if saveOutput:
    OUTDIR = "./"
    print( "Saving output to ",OUTDIR)
    if not os.path.isdir( OUTDIR ):
        os.makedirs( OUTDIR )


### -----------------------------------
## MAIN LOOP

# save resolution etc for a given loop to file:
if saveOutput:
    f = open(OUTDIR+filename,'w')
    f.write("wavenr_min,wavenr_max,dwave [cm^-1] = %.4f,%.4f,%.4f" % (wavenr_min,wavenr_max,dwavenr) )
    f.write("\n")
    f.write("N_press = %.1f" % N_press )
    f.write("\n")
    f.write("\n")
    f.write("Ts [K],\tps [bar],\tolr [W/m2],\tsurface contribution to olr [W/m2],\tTransmission int[T dBdT(Ts) dn]/int[ dBdT(Ts) dn],\tSimple feedback model int[dB/dT*T]")
    f.write("\n")

    f.close()

    ## main loop here
    for Ts in Ts_grid:
        f = open(OUTDIR+filename,'a')

        # setup grid:
        g = pyrads.SetupGrids.make_grid( Ts,Tstrat,N_press,wavenr_min,wavenr_max,dwavenr,params, RH=params.RH )

        # compute optical thickness:
        #   -> this is the computationally most intensive step
        g.tau,g.omega = pyrads.OpticalThickness.compute_tau_omega_H2ON2(g.p,g.T,g.q,g,params, RH=params.RH )
        
        # compute Planck functions etc:
        #   -> here: fully spectrally resolved!
        T_2D = np.tile( g.T, (g.Nn,1) ).T               # [press x wave]
        g.B_surf = np.pi* pyrads.Planck.Planck_n( g.n,Ts )     # [wave]
        g.B = np.pi* pyrads.Planck.Planck_n( g.wave, T_2D )    # [press x wave]

        # compute OLR etc:
        olr_spec = pyrads.Get_Fluxes.Fplus_alternative(0,g) # (spectrally resolved=irradiance)
        olr = simps(olr_spec,g.n)

        # compute fraction of surface flux that makes it to space
        surf_spec = g.B_surf * np.exp(-g.tau[-1,:])
        surf = simps(surf_spec,x=g.n)

        # compute spectrally averaged transmission function...
        weight = np.pi* pyrads.Planck.dPlanckdT_n( g.n,Ts )
        trans = trapz( np.exp(-g.tau[-1,:]) * weight,x=g.n ) / trapz( weight,x=g.n )

        # Simple feedback model (like above, without normalization)
        weight = np.pi* pyrads.Planck.dPlanckdT_n( g.n,Ts )
        lam = trapz( np.exp(-g.tau[-1,:]) * weight,x=g.n )

        print( "\n",Ts,g.ps/1e5,olr,surf, "\n")

        f.write("%.2f,\t%.4f,\t%.8f,\t%.8f,\t%.8f,\t%.8f" % (Ts,g.ps/1e5,olr,surf,trans,lam) )
        f.write("\n")

        f.close()
