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

N_press = 15       # for testing only!
dwavenr = 0.1     #  for testing only!

#N_press = 60       #
wavenr_min = 0.1   # [cm^-1]
wavenr_max = 3500. #
#dwavenr = 0.01     #

Tstrat = 150.      # stratospheric temperature

# ---
## setup range of temperatures, and if/where output is saved to:

Ts = 300.

### -----------------------------------
## MAIN LOOP

print( "wavenr_min,wavenr_max,dwave [cm^-1] = %.4f,%.4f,%.4f" % (wavenr_min,wavenr_max,dwavenr))
print( "\n")
print( "N_press = %.1f" % N_press)
print( "\n")
print( "Surface temperature = %.1f K" % Ts)


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

print( "OLR = ",olr)
