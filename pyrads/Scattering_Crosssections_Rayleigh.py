from __future__ import division, print_function, absolute_import
import numpy as np
from . import phys
import os

'''
Implement Rayleigh scattering cross-sections.
'''

### -----------------------------------
### Global definitions here

#Path to the datasets
datapath = '/'.join( os.path.abspath(__file__).split('/')[:-2] ) + '/DATA/'  # !
goldblatt_tableS1 = datapath + "Goldblatt2013_archive_revised/Goldblatt_TableS1.csv"

# # hard-code data?
# wavelen = np.array([1.00E-02,1.50E-01,1.60E-01,1.70E-01,1.80E-01,1.90E-01,2.00E-01,2.10E-01,2.20E-01,
#                     2.30E-01,2.40E-01,2.50E-01,2.60E-01,2.70E-01,2.80E-01,2.90E-01,3.00E-01,3.10E-01,
#                     3.20E-01,3.30E-01,3.40E-01,3.50E-01,3.60E-01,3.70E-01,3.80E-01,3.90E-01,4.00E-01, 
#                     4.10E-01,4.20E-01,4.30E-01,4.40E-01,4.50E-01,4.60E-01,4.70E-01,4.80E-01,4.90E-01,
#                     5.00E-01,5.10E-01,5.20E-01,5.30E-01,5.40E-01,5.50E-01,5.60E-01,5.70E-01,5.80E-01,
#                     5.90E-01,6.00E-01,6.10E-01,6.20E-01,6.30E-01,6.40E-01,6.50E-01,6.60E-01,6.70E-01,
#                     6.80E-01,6.90E-01,7.00E-01,7.10E-01,7.20E-01,7.30E-01,7.40E-01,7.50E-01,7.60E-01,
#                     7.70E-01,7.80E-01,7.90E-01,8.00E-01,8.10E-01,8.20E-01,8.30E-01,8.40E-01,8.50E-01,
#                     8.60E-01,8.70E-01,8.80E-01,8.90E-01,9.00E-01,9.10E-01,9.20E-01,9.30E-01,9.40E-01,
#                     9.50E-01,9.60E-01,9.70E-01,9.80E-01,9.90E-01,1.00E+00,1.10E+00,1.20E+00,1.30E+00,
#                     1.40E+00,1.50E+00,1.60E+00,1.70E+00,1.80E+00,1.90E+00,2.00E+00,3.00E+00,4.00E+00,
#                     5.00E+00,1.00E+01])
# cAir = ...
# cH2O = ...


### -----------------------------------
### Implement fits from Goldblatt et al (2013), Table S1
# 
# Data are given as m^2/molecule.
#     careful: the description given in table caption doesn't match these units.
# Linearly interpolate to get empirical constant 'c', then use analytical ~1/lambda^4.
#
# INPUT:
#     wave [NEED cm^-1!]
# OUTPUT:
#     kappaSca_H2O    [m^2/kg]
#
# NOTE: here, load table each time kappa is called.
# NOTE 2: the eqn in Table S1 has an error '128*pi^{5/3}' -> '128*pi^{5}/3'!

def get_KappaSca_H2O(wavenr):
    wavelen = 1e-2/wavenr            # [cm^-1] -> [m]

    # read in data:
    data = np.genfromtxt(goldblatt_tableS1,skip_header=1,delimiter=',')

    wavelen_data = data[:,0] * 1e-6  # [micron] -> [m]
    cH2O_data = data[:,2] * phys.N_avogadro/phys.H2O.MolecularWeight*1e3  # [m^2/molecule] * (molec/mole*mole/g*g/kg) = [m^2/kg]

    cH2O = np.interp( wavelen, wavelen_data, cH2O_data )
    kappaH2O = 128./3. *np.pi**5 *cH2O *(1./wavelen)**4.

    return kappaH2O


def get_KappaSca_Air(wavenr):
    wavelen = 1e-2/wavenr            # [cm^-1] -> [m]

    # read in data:
    data = np.genfromtxt(goldblatt_tableS1,skip_header=1,delimiter=',')

    wavelen_data = data[:,0] * 1e-6  # [micron] -> [m]
    cAir_data = data[:,1] * phys.N_avogadro/phys.air.MolecularWeight*1e3  # [m^2/molecule] * (molec/mole*mole/g*g/kg) = [m^2/kg]

    cAir = np.interp( wavelen, wavelen_data, cAir_data )
    kappaAir = 128./3. *np.pi**5 *cAir *(1./wavelen)**4.

    return kappaAir


### -----------------------------------
### Implements simple analytical formula from Dalgarno & Williams (1962), Eq.3:
# Note: in their formulation, [lambda] = Angstrom!
# Note: these values reproduce the crosssections in Table 5.2 in PoPC pretty well.

def get_KappaSca_H2(wavenr):
    wavelen = 1e-2/wavenr            # [cm^-1] -> [m]
    
    # for fit: [input]=A, [output]=cm^2
    get_sigma = lambda x: (8.14e-13/(x**4) + 1.28e-6/(x**6) + 1.61/(x**8)  )

    # [cm^2/molecule] * (molec/mole*mole/g) * (m^2/cm^2*g/kg) = [m^2/kg]
    kappaH2 = get_sigma(wavelen*1e10) * phys.N_avogadro/phys.H2.MolecularWeight *1e-1

    return kappaH2


### -----------------------------------
### Implements fits to Rayleigh scattering for CO2 from PoPC!
###  -> very simple fits

def get_KappaSca_PoPC(wavenr,gasname):
    wavelength = 1e-2/wavenr            # [cm^-1] -> [m]
    

    # / ------ /
    # scale data to given wavelength from
    # nearest reference point that is closest to the wavelength
    # / ------ /
    # Actually, the above takes much too long:
    # just pick the reference value for 1 micron

    # Reference value: H2 xsec at 1 micron
    wavelength0 = 1e-6        # [m]
    chi_per_mass0 = 2.49e-6

    x = wavelength / wavelength0
    if gasname == "H2":
        relative_chi = 1.
        relative_chi_per_mass = 1.
    elif gasname == "H2O":
        relative_chi = 3.3690
        relative_chi_per_mass = 0.3743
    elif gasname == "He":
        relative_chi = 0.0641
        relative_chi_per_mass = 0.0321
    elif gasname == "air":
        relative_chi = 4.4459
        relative_chi_per_mass = 0.3066
    elif gasname == "N2":
        relative_chi = 4.6035
        relative_chi_per_mass = 0.3288
    elif gasname == "O2":
        relative_chi = 3.8634
        relative_chi_per_mass = 0.2415
    elif gasname == "CO2":
        relative_chi = 10.5611
        relative_chi_per_mass = 0.4800
    elif gasname == "NH3":
        relative_chi = 7.3427
        relative_chi_per_mass = 0.8638
    elif gasname == "CH4":
        relative_chi = 10.1509
        relative_chi_per_mass = 1.2689
    else:
        print( "Gas not recognized!" )

    kappa = chi_per_mass0 * relative_chi_per_mass / (x**4)

    return kappa
