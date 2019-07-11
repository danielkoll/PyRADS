'''
***********************************************************
This script computes absorption coefficients and
optical thicknesses.
***********************************************************
'''
from __future__ import division, print_function, absolute_import
import numpy as np
from .Absorption_Crosssections_HITRAN2016 import getKappa_HITRAN
from . import Absorption_Continuum_MTCKD
from .Absorption_Continuum_MTCKD import get_H2OContinuum
from . import hitran_cia_fast as hitran_cia
from .Scattering_Crosssections_Rayleigh import get_KappaSca_H2O,get_KappaSca_Air,get_KappaSca_H2

from scipy.integrate import cumtrapz
from .Thermodynamics import convert_molar_to_mass_ratio


# ---
#
# Compute Kappa, return grid with Tau and Omega:
#   {{ INPUT:   (p,T,q) 1d vectors; (n0,n1,dn) floats
#   {{ OUTPUT:  (tau,omega)=2d

def compute_tau_omega_H2ON2(p,T,q,grid,params,RH=1.):

    kappa = np.zeros( (grid.Np,grid.Nn) )
    omega = np.zeros( (grid.Np,grid.Nn) )

    for pres,temp,q_H2O in zip(p,T,q):
        print( "compute kappa at p,T = ",pres,temp )  ## FOR TESTING
        p_H2O = RH * params.esat(temp)  # ...
        q_N2 = 1.-(q_H2O)   # no other species!
        R_mean = q_H2O*params.Rv + q_N2*params.R
        rho = pres/(R_mean*temp)
        
        # --
        # add lines
        # a) for wavenumbers < 10000cm, use the MTCKD continuum:
        #    to be consistent, also truncate lines at 25cm-1 from center.
        # b) for wavenumbers > 10000cm, no continuum data:
        #    instead truncate lines x decay widths away from center.
        # c) for wavenumbers > 40000cm, H2O dissociates anyway, so just skip ahead.
        #    (note, there is no line data in HITRAN2016 above ~25000 cm, so use 30000 cm here) 

        linewid_sw = 1000.
        
        if grid.n.max() < 10000.:
            kappaH2O = getKappa_HITRAN(grid.n,grid.n0,grid.n1,grid.dn,"H2O",press=pres, \
                                       temp=temp,broadening="mixed", lineWid=25., \
                                       cutoff_option="fixed",remove_plinth=True, \
                                       press_self=p_H2O)
            
        elif grid.n.min() > 30000.:
            kappaH2O = grid.n * 0.

        elif grid.n.min() > 10000.:
            kappaH2O = getKappa_HITRAN(grid.n,grid.n0,grid.n1,grid.dn,"H2O",press=pres, \
                                       temp=temp,broadening="mixed", lineWid=linewid_sw, \
                                       cutoff_option="relative",remove_plinth=False, \
                                       press_self=p_H2O)
            
        else:
            kappa_fixed = getKappa_HITRAN(grid.n,grid.n0,grid.n1,grid.dn,"H2O",press=pres, \
                                   temp=temp,broadening="mixed", lineWid=25., \
                                   cutoff_option="fixed",remove_plinth=True, \
                                   press_self=p_H2O)
            kappa_relative = getKappa_HITRAN(grid.n,grid.n0,grid.n1,grid.dn,"H2O",press=pres, \
                                   temp=temp,broadening="mixed", lineWid=linewid_sw, \
                                   cutoff_option="relative",remove_plinth=False, \
                                   press_self=p_H2O)

            kappaH2O = grid.n * 0.
            kappaH2O[grid.n<10000.] = kappa_fixed[grid.n<10000.]
            kappaH2O[grid.n>=10000.] = kappa_relative[grid.n>=10000.]
        
        # --
        # add H2O-H2O, H2O-N2 continuum:
        kappaH2O_cont = get_H2OContinuum(grid.n,temp,pres,p_H2O, \
                               exe_file=Absorption_Continuum_MTCKD.mtckd_exe_H2O_N2)  # careful!?
        kappaH2O_cont[grid.n>=10000.] = 0.   # ...MTCKD returns NaN's past 10000cm-1
        
        # add N2-N2 cont from HITRAN CIA.
        # note: this returns qN2 * kappa!
        kappaN2_self = hitran_cia.get_kappa_cia("N2","N2",grid.n,temp,q_N2,q_N2,rho)
        
        # --
        # add scattering:
        kappaH2O_sca = get_KappaSca_H2O(grid.n)
        kappaAir_sca = get_KappaSca_Air(grid.n)    # FOR AIR BACKGROUND

        # --
        # add a really simple UV parametrization:
        #   based on MPI-Mainz Atlas, say:
        #   kappa_H2O >~ 1e-18 cm/molec = 3e3 m2/kg beyond ~200nm=50000cm-1
        #   kappa_{N2,H2}>~3e3 m2/kg beyond ~100nm=100000cm-1 -> ignore
        kappaUV = 3e3
        nUV = 50000.

        kappaH2O_UV = grid.n * 0.
        kappaH2O_UV[grid.n>=nUV] = kappaUV
        
        # --
        # here: careful to use correct mass ratios!
        kappa[ p==pres,: ] = q_H2O*(kappaH2O + kappaH2O_cont + kappaH2O_sca + kappaH2O_UV) + \
                             (1.-q_H2O)*kappaAir_sca + kappaN2_self   # save
        omega[ p==pres,: ] = ( q_H2O*kappaH2O_sca + (1.-q_H2O)*kappaAir_sca ) / kappa[ p==pres,: ]  # save
    print( "Done! \n" )

    # Integrate to get optical thickness:
    p2d = np.tile( p,(grid.Nn,1) ).T
    tau = 1./(params.g) * cumtrapz( kappa,x=p2d,initial=0.,axis=0 )  # DONT INCLUDE cosTHETA

    return tau,omega



# ---
#
# Compute Kappa, return grid with Tau and Omega:
#   {{ INPUT:   (p,T,q) 1d vectors; (n0,n1,dn) floats
#   {{ OUTPUT:  (tau,omega)=2d

def compute_tau_omega_H2OH2(p,T,q,grid,params,RH=1.):

    kappa = np.zeros( (grid.Np,grid.Nn) )
    omega = np.zeros( (grid.Np,grid.Nn) )

    for pres,temp,q_H2O in zip(p,T,q):
        print( "compute kappa at p,T = ",pres,temp )  ## FOR TESTING
        p_H2O = RH * params.esat(temp)  # ...
        q_H2 = 1.-(q_H2O)   # no other species!
        R_mean = q_H2O*params.Rv + q_H2*params.R
        rho = pres/(R_mean*temp)

        # --
        # add lines
        # a) for wavenumbers < 10000cm, use the MTCKD continuum:
        #    to be consistent, also truncate lines at 25cm-1 from center.
        # b) for wavenumbers > 10000cm, no continuum data:
        #    instead truncate lines x decay widths away from center.

        linewid_sw = 1000.
        
        if grid.n.max() < 10000.:
            kappaH2O = getKappa_HITRAN(grid.n,grid.n0,grid.n1,grid.dn,"H2O",press=pres, \
                                       temp=temp,broadening="mixed", lineWid=25., \
                                       cutoff_option="fixed",remove_plinth=True, \
                                       press_self=p_H2O)
            
        elif grid.n.min() > 10000.:
            kappaH2O = getKappa_HITRAN(grid.n,grid.n0,grid.n1,grid.dn,"H2O",press=pres, \
                                       temp=temp,broadening="mixed", lineWid=linewid_sw, \
                                       cutoff_option="relative",remove_plinth=False, \
                                       press_self=p_H2O)
            
        else:
            kappa_fixed = getKappa_HITRAN(grid.n,grid.n0,grid.n1,grid.dn,"H2O",press=pres, \
                                   temp=temp,broadening="mixed", lineWid=25., \
                                   cutoff_option="fixed",remove_plinth=True, \
                                   press_self=p_H2O)
            kappa_relative = getKappa_HITRAN(grid.n,grid.n0,grid.n1,grid.dn,"H2O",press=pres, \
                                   temp=temp,broadening="mixed", lineWid=linewid_sw, \
                                   cutoff_option="relative",remove_plinth=False, \
                                   press_self=p_H2O)

            kappaH2O = grid.n * 0.
            kappaH2O[grid.n<10000.] = kappa_fixed[grid.n<10000.]
            kappaH2O[grid.n>=10000.] = kappa_relative[grid.n>=10000.]
        
        # --
        # add continuum (only where it is defined):
        kappaH2O_cont = get_H2OContinuum(grid.n,temp,pres,p_H2O, \
                               exe_file=Absorption_Continuum_MTCKD.mtckd_exe_H2O_N2)  # careful!?
        kappaH2O_cont[grid.n>=10000.] = 0.

        # add H2-H2 cont from HITRAN CIA.
        # note: this returns qH2 * kappa!
        kappaH2_self = hitran_cia.get_kappa_cia("H2","H2",grid.n,temp,q_H2,q_H2,rho)
        
        # --
        # add scattering:
        kappaH2O_sca = get_KappaSca_H2O(grid.n)
        kappaH2_sca = get_KappaSca_H2(grid.n)    # FOR H2 BACKGROUND

        # --
        # add a really simple UV parametrization:
        #   based on MPI-Mainz Atlas, say:
        #   kappa_H2O >~ 1e-18 cm/molec = 3e3 m2/kg beyond ~200nm=50000cm-1
        #   kappa_{N2,H2}>~3e3 m2/kg beyond ~100nm=100000cm-1 -> ignore
        kappaUV = 3e3
        nUV = 50000.

        kappaH2O_UV = grid.n * 0.
        kappaH2O_UV[grid.n>=nUV] = kappaUV
        
        # --
        # here: careful to use correct mass ratios!
        kappa[ p==pres,: ] = q_H2O*(kappaH2O + kappaH2O_cont + kappaH2O_sca + kappaH2O_UV) + \
                             (1.-q_H2O)*kappaH2_sca + kappaH2_self   # save
        omega[ p==pres,: ] = ( q_H2O*kappaH2O_sca + (1.-q_H2O)*kappaH2_sca ) / kappa[ p==pres,: ]  # save
    print( "Done! \n" )

    # Integrate to get optical thickness:
    p2d = np.tile( p,(grid.Nn,1) ).T
    tau = 1./(params.g) * cumtrapz( kappa,x=p2d,initial=0.,axis=0 )  # DONT INCLUDE cosTHETA

    return tau,omega







# ---
#
# Compute Kappa, return grid with Tau and Omega:
#   {{ INPUT:   (p,T,q) 1d vectors; (n0,n1,dn) floats
#   {{ OUTPUT:  (tau,omega)=2d

def compute_tau_omega_H2OCO2(p,T,q,grid,params,RH=1.):

    kappa = np.zeros( (grid.Np,grid.Nn) )
    omega = np.zeros( (grid.Np,grid.Nn) )

    for pres,temp,q_H2O in zip(p,T,q):
        print( "compute kappa at p,T = ",pres,temp )  ## FOR TESTING
        p_H2O = RH * params.esat(temp)  # ...
        p_CO2 = pres - p_H2O
        q_CO2 = 1.-(q_H2O)   # no other species!
        R_mean = q_H2O*params.Rv + q_CO2*params.R
        rho = pres/(R_mean*temp)
        
        # --
        # add lines
        # a) for wavenumbers < 10000cm, use the MTCKD continuum:
        #    to be consistent, also truncate lines at 25cm-1 from center.
        # b) for wavenumbers > 10000cm, no continuum data:
        #    instead truncate lines x decay widths away from center.

        linewid_sw = 1000.
        
        if grid.n.max() < 10000.:
            kappaH2O = getKappa_HITRAN(grid.n,grid.n0,grid.n1,grid.dn,"H2O",press=pres, \
                                       temp=temp,broadening="mixed", lineWid=25., \
                                       cutoff_option="fixed",remove_plinth=True, \
                                       press_self=p_H2O)

        elif grid.n.min() > 10000.:
            kappaH2O = getKappa_HITRAN(grid.n,grid.n0,grid.n1,grid.dn,"H2O",press=pres, \
                                       temp=temp,broadening="mixed", lineWid=linewid_sw, \
                                       cutoff_option="relative",remove_plinth=False, \
                                       press_self=p_H2O)
            
        else:
            kappa_fixed = getKappa_HITRAN(grid.n,grid.n0,grid.n1,grid.dn,"H2O",press=pres, \
                                   temp=temp,broadening="mixed", lineWid=25., \
                                   cutoff_option="fixed",remove_plinth=True, \
                                   press_self=p_H2O)
            kappa_relative = getKappa_HITRAN(grid.n,grid.n0,grid.n1,grid.dn,"H2O",press=pres, \
                                   temp=temp,broadening="mixed", lineWid=linewid_sw, \
                                   cutoff_option="relative",remove_plinth=False, \
                                   press_self=p_H2O)

            kappaH2O = grid.n * 0.
            kappaH2O[grid.n<10000.] = kappa_fixed[grid.n<10000.]
            kappaH2O[grid.n>=10000.] = kappa_relative[grid.n>=10000.]


        # for CO2 lines: always use some relative cutoff!
        kappaCO2 = getKappa_HITRAN(grid.n,grid.n0,grid.n1,grid.dn,"CO2",press=pres, \
                                       temp=temp,broadening="mixed", lineWid=1000., \
                                       cutoff_option="relative",remove_plinth=False, \
                                       press_self=p_CO2)

        
        # --
        # add H2O-H2O, H2O-CO2 continuum:
        #    ** here: assume H2O-CO2 is same as H2O-N2 **
        kappaH2O_cont = get_H2OContinuum(grid.n,temp,pres,p_H2O, \
                               exe_file=Absorption_Continuum_MTCKD.mtckd_exe_H2O_N2)  # careful!?
        kappaH2O_cont[grid.n>=10000.] = 0.   # ...MTCKD returns NaN's past 10000cm-1
        
        # add CO2-CO2 cont from Ray's fits.
        kappaCO2_self = Absorption_Continuum.get_CO2SelfContinuum(grid.n,temp,p_CO2)
        
        # --
        # add scattering:
        kappaH2O_sca = get_KappaSca_H2O(grid.n)
        kappaAir_sca = get_KappaSca_Air(grid.n)    # HERE: I should used CO2 instead of air!

        # --
        # add a really simple UV parametrization:
        #   based on MPI-Mainz Atlas, say:
        #   kappa_H2O >~ 1e-18 cm/molec = 3e3 m2/kg beyond ~200nm=50000cm-1
        #   kappa_{N2,H2}>~3e3 m2/kg beyond ~100nm=100000cm-1 -> ignore
        kappaUV = 3e3
        nUV = 50000.

        kappaH2O_UV = grid.n * 0.
        kappaH2O_UV[grid.n>=nUV] = kappaUV
        
        # --
        # here: careful to use correct mass ratios!
        kappa[ p==pres,: ] = q_H2O*(kappaH2O + kappaH2O_cont + kappaH2O_sca + kappaH2O_UV) + \
                             (1.-q_H2O)*(kappaAir_sca + kappaCO2 + kappaCO2_self)   # save
        omega[ p==pres,: ] = ( q_H2O*kappaH2O_sca + (1.-q_H2O)*kappaAir_sca ) / kappa[ p==pres,: ]  # save
    print( "Done! \n" )

    # Integrate to get optical thickness:
    p2d = np.tile( p,(grid.Nn,1) ).T
    tau = 1./(params.g) * cumtrapz( kappa,x=p2d,initial=0.,axis=0 )  # DONT INCLUDE cosTHETA

    return tau,omega

