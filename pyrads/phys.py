from __future__ import division, print_function, absolute_import
### THIS SCRIPT IS TAKEN FROM THE COURSEWARE OF
###  Pierrehumbert, 2010, Principles of Planetary Climate
###
###  Downloaded from:
###  https://geosci.uchicago.edu/~rtp1/PrinciplesPlanetaryClimate/Courseware/coursewarePortal.html

import math
#from .ClimateUtilities import * #To get the math methods routines
#
#All units are mks units
#
#ToDo:
#       * Add conversion factors (calories to joule,Megatons, etc.)
#       * Unit conversion calculating object
#       * Database of interesting constants (e.g energy and C
#          content of coal)
#       * Find a better way to organize database of contants,
#         allowing indexing,unit information and description
#       *Add the rest of the gases to the database, and find a
#        better way to index them and provide help. Note that
#        a printed table is more "transparent" in terms of what
#        data is available. Find a way of providing a similar
#        tabular summary of available data, and perhaps values
#        This applies to the planet data table as well, and perhaps
#        even to the table of physical constants.
#
#       *Finish putting data into the gas database,
#        perhaps including Van der Waals, Shomate and Antoine
#        coefficients, at least in selected cases
#
#

#-------------Basic physical constants-------------------
#
#The following five are the most accurate 1986 values
#
h = 6.626075540e-34    #Planck's constant
c = 2.99792458e8       #Speed of light
k =1.38065812e-23      #Boltzman thermodynamic constant
sigma = 5.67051196e-8  #Stefan-Boltzman constant
G = 6.67428e-11        #Gravitational constant (2006 measurements)
#

#-----------Thermodynamic constants------------
#Following will come out in J/(deg kmol), so
#that dividing Rstar by molecular weight gives
#gas constant appropriate for mks units
N_avogadro = 6.022136736e23  #Avogadro's number
Rstar = 1000.*k*N_avogadro   #Universal gas constant
#


#----------Properties of gases-----------------
#The following are approximate mean values
#for "normal" temperatures and pressures, suitable only
#for rough calculations.
#

#This class allows convenient access
#to the basic thermodynamic properties of
#a gas, and selected properties of its
#solid and liquid phases. The properties
#do not need to be specified in the __init__
#method, since they can be created dynamically.
#We specify them anyway, and set them to None,
#as a guide to the naming conventions for those
#creating their own gas objects. A utility is also
#available which will turn a LaTeX formatted thermodynamic table
#into a Python script defining new gas objects.
#
#
#
#ToDo: *Add more documentation and help features
#
#      *provide dictionary of properties and units.
#      *Add some methods or lists that make it easier
#       for the user to create new gas objects and insert
#       the data
class gas:
    '''
    A gas object stores thermodynamic data for
    a gas, and selected properties of its condensed
    phases.  You can create a gas object for gas
    G by executing:
             G = gas()
    and then setting the attributes individually (see
    below for explanation of names).  A collection
    of gas objects for common gases is provided as part
    of the phys module.  These gas objects were not actually
    created "by hand" but rather using a utility script that
    automatically translates a LaTeX formatted table into a
    Python script defining the objects.


    Attributes of a gas object:
        CriticalPointT:  Critical point temperature (K)
        CriticalPointP:  Critical point pressure (Pa)
        TriplePointT: Triple point temperature (K)
        TriplePointP: Triple point pressure (Pa)
        L_vaporization_BoilingPoint: Latent heat of vaporization (J/kg)
            at boiling point.  The so-called "boiling point" is
            the temperature at which the saturation vapor pressure
            equals 1 atmosphere (1.013 bar). For CO2, the "boiling point"
            occurs below the triple point temperature, so the condensed
            phase would not be a liquid. Hence, for CO2 the
            latent heat is given at the arbitrary reference point
            of 253K and 29Pa.
        L_vaporization_TriplePoint: Latent heat of vaporization (J/kg)
            at the triple point
        L_fusion: Latent heat of fusion (J/kg) at the triple point
        L_sublimation: Latent heat of sublimation (J/kg) at triple point
        rho_liquid_BoilingPoint Liquid phase density (kg/m**3)
            at the boiling point
        rho_liquid_TriplePoint: Liquid phase density (kg/m**3)
            at the triple point
        rho_solid: Solid phase density (kg/m**3) at (or sometimes near)
            the triple point
        cp: Gas phase specific heat (J/(kg K)), at 298K and 1 bar
        gamma: ratio of specific heat at constant pressure
            to specific heat at constant volume. (Generally
            stated at 298K and 1bar)
        MolecularWeight: Molecular weight of the dominant isotope
        name: Name of the gas
        formula: Chemical formula (e.g. 'CH4')

        L_vaporization: Default value to use for latent heat of
            vaporization.  Set to triple point value, if available,
            else to boiling point value
        rho_liquid: Default value to use for liquid phase density.
            Set to triple point value if available otherwise
            set to boiling point value

        R: Gas constant for the individual gas. Computed from
            other data as Rstar/MolecularWeight, when the update()
            method is called
        Rcp: The adiabatic exponent R/cp. Computed from other
            data when the update() method is called.

    '''

    #The __repr__ method allows us to print out
    #a help string when the user types the name
    #of a gas object
    def __repr__(self):
        firstline =\
        'This gas object contains thermodynamic data on %s\n'%self.formula
        secondline = \
        'Type \"help(gas)\" for more information\n'
        return firstline+secondline
    def __init__(self):
        self.CriticalPointT = None
        self.CriticalPointP = None
        self.TriplePointT = None
        self.TriplePointP = None
        self.L_vaporization_BoilingPoint = None
        self.L_vaporization_TriplePoint = None
        self.L_fusion = None
        self.L_sublimation = None
        self.rho_liquid_BoilingPoint = None
        self.rho_liquid_TriplePoint = None
        self.rho_solid = None
        self.cp = None
        self.gamma = None
        self.MolecularWeight = None
        self.name = None
        self.formula = None
        #
        #Default values for latent heat and liquid
        #density.  The triple point value is set
        #as the default if it is available, otherwise
        #the boiling point values are used.
        self.L_vaporization= None
        self.rho_liquid= None
        #
        #Computed quantities
        #
        self.R = None
        self.Rcp = None
    #Function to compute derived properties
    def update(self):
        self.R = Rstar/self.MolecularWeight
        self.Rcp = self.R/self.cp

#Set properties of individual gases
#
#Thermodynamic properties of dry Earth air
air = gas()
air.name = 'Earth Air'
air.cp = 1004.
air.MolecularWeight = 28.97
air.gamma = 1.4003

#------------------------
H2O = gas()
H2O.CriticalPointT = 6.471000e+02
H2O.CriticalPointP = 2.210000e+07
H2O.TriplePointT = 2.731500e+02
H2O.TriplePointP = 6.110000e+02
H2O.L_vaporization_BoilingPoint = 2.255000e+06
H2O.L_vaporization_TriplePoint = 2.493000e+06
H2O.L_fusion = 3.340000e+05
H2O.L_sublimation = 2.840000e+06
H2O.rho_liquid_BoilingPoint = 9.584000e+02
H2O.rho_liquid_TriplePoint = 9.998700e+02
H2O.rho_solid = 9.170000e+02
H2O.cp = 1.847000e+03
H2O.gamma = 1.331000e+00
H2O.MolecularWeight = 1.800000e+01
H2O.name = 'Water'
H2O.formula = 'H2O'
H2O.L_vaporization=2.493000e+06
H2O.rho_liquid=9.998700e+02
#------------------------
CH4 = gas()
CH4.CriticalPointT = 1.904400e+02
CH4.CriticalPointP = 4.596000e+06
CH4.TriplePointT = 9.067000e+01
CH4.TriplePointP = 1.170000e+04
CH4.L_vaporization_BoilingPoint = 5.100000e+05
CH4.L_vaporization_TriplePoint = 5.360000e+05
CH4.L_fusion = 5.868000e+04
CH4.L_sublimation = 5.950000e+05
CH4.rho_liquid_BoilingPoint = 4.502000e+02
CH4.rho_liquid_TriplePoint = None
CH4.rho_solid = 5.093000e+02
CH4.cp = 2.195000e+03
CH4.gamma = 1.305000e+00
CH4.MolecularWeight = 1.600000e+01
CH4.name = 'Methane'
CH4.formula = 'CH4'
CH4.L_vaporization=5.360000e+05
CH4.rho_liquid=4.502000e+02
#------------------------
CO2 = gas()
CO2.CriticalPointT = 3.042000e+02
CO2.CriticalPointP = 7.382500e+06
CO2.TriplePointT = 2.165400e+02
CO2.TriplePointP = 5.185000e+05
CO2.L_vaporization_BoilingPoint = None
CO2.L_vaporization_TriplePoint = 3.970000e+05
CO2.L_fusion = 1.960000e+05
CO2.L_sublimation = 5.930000e+05
CO2.rho_liquid_BoilingPoint = 1.032000e+03
CO2.rho_liquid_TriplePoint = 1.110000e+03
CO2.rho_solid = 1.562000e+03
CO2.cp = 8.200000e+02
CO2.gamma = 1.294000e+00
CO2.MolecularWeight = 4.400000e+01
CO2.name = 'Carbon Dioxide'
CO2.formula = 'CO2'
CO2.L_vaporization=3.970000e+05
CO2.rho_liquid=1.110000e+03
#------------------------
N2 = gas()
N2.CriticalPointT = 1.262000e+02
N2.CriticalPointP = 3.400000e+06
N2.TriplePointT = 6.314000e+01
N2.TriplePointP = 1.253000e+04
N2.L_vaporization_BoilingPoint = 1.980000e+05
N2.L_vaporization_TriplePoint = 2.180000e+05
N2.L_fusion = 2.573000e+04
N2.L_sublimation = 2.437000e+05
N2.rho_liquid_BoilingPoint = 8.086000e+02
N2.rho_liquid_TriplePoint = None
N2.rho_solid = 1.026000e+03
N2.cp = 1.037000e+03
N2.gamma = 1.403000e+00
N2.MolecularWeight = 2.800000e+01
N2.name = 'Nitrogen'
N2.formula = 'N2'
N2.L_vaporization=2.180000e+05
N2.rho_liquid=8.086000e+02
#------------------------
O2 = gas()
O2.CriticalPointT = 1.545400e+02
O2.CriticalPointP = 5.043000e+06
O2.TriplePointT = 5.430000e+01
O2.TriplePointP = 1.500000e+02
O2.L_vaporization_BoilingPoint = 2.130000e+05
O2.L_vaporization_TriplePoint = 2.420000e+05
O2.L_fusion = 1.390000e+04
O2.L_sublimation = 2.560000e+05
O2.rho_liquid_BoilingPoint = 1.141000e+03
O2.rho_liquid_TriplePoint = 1.307000e+03
O2.rho_solid = 1.351000e+03
O2.cp = 9.160000e+02
O2.gamma = 1.393000e+00
O2.MolecularWeight = 3.200000e+01
O2.name = 'Oxygen'
O2.formula = 'O2'
O2.L_vaporization=2.420000e+05
O2.rho_liquid=1.307000e+03
#------------------------
H2 = gas()
H2.CriticalPointT = 3.320000e+01
H2.CriticalPointP = 1.298000e+06
H2.TriplePointT = 1.395000e+01
H2.TriplePointP = 7.200000e+03
H2.L_vaporization_BoilingPoint = 4.540000e+05
H2.L_vaporization_TriplePoint = None
H2.L_fusion = 5.820000e+04
H2.L_sublimation = None
H2.rho_liquid_BoilingPoint = 7.097000e+01
H2.rho_liquid_TriplePoint = None
H2.rho_solid = 8.800000e+01
H2.cp = 1.423000e+04
H2.gamma = 1.384000e+00
H2.MolecularWeight = 2.000000e+00
H2.name = 'Hydrogen'
H2.formula = 'H2'
H2.L_vaporization=4.540000e+05
H2.rho_liquid=7.097000e+01
#------------------------
He = gas()
He.CriticalPointT = 5.100000e+00
He.CriticalPointP = 2.280000e+05
He.TriplePointT = 2.170000e+00
He.TriplePointP = 5.070000e+03
He.L_vaporization_BoilingPoint = 2.030000e+04
He.L_vaporization_TriplePoint = None
He.L_fusion = None
He.L_sublimation = None
He.rho_liquid_BoilingPoint = 1.249600e+02
He.rho_liquid_TriplePoint = None
He.rho_solid = 2.000000e+02
He.cp = 5.196000e+03
He.gamma = 1.664000e+00
He.MolecularWeight = 4.000000e+00
He.name = 'Helium'
He.formula = 'He'
He.L_vaporization=2.030000e+04
He.rho_liquid=1.249600e+02
#------------------------
NH3 = gas()
NH3.CriticalPointT = 4.055000e+02
NH3.CriticalPointP = 1.128000e+02
NH3.TriplePointT = 1.954000e+02
NH3.TriplePointP = 6.100000e+03
NH3.L_vaporization_BoilingPoint = 1.371000e+06
NH3.L_vaporization_TriplePoint = 1.658000e+06
NH3.L_fusion = 3.314000e+05
NH3.L_sublimation = 1.989000e+06
NH3.rho_liquid_BoilingPoint = 6.820000e+02
NH3.rho_liquid_TriplePoint = 7.342000e+02
NH3.rho_solid = 8.226000e+02
NH3.cp = 2.060000e+03
NH3.gamma = 1.309000e+00
NH3.MolecularWeight = 1.700000e+01
NH3.name = 'Ammonia'
NH3.formula = 'NH3'
NH3.L_vaporization=1.658000e+06
NH3.rho_liquid=7.342000e+02
#------------------------

#------------------------
#Synonym for H2O
water = H2O
#Make a list of all the gases
#
#This clever little fragment uses the fact that
#Python can execute any string as a Python statement,
#in order to find all the gases and build a list of them.
#I don't know if there is a more straightforward way to
#get a list of all the objects of a certain type, but this
#works.  Some of the trickery below is needed because
#the dir() command returns a list of strings, which
#give the names of the objects. It doesn't give the objects
#themselves.
gases = []
for ob in dir():
    exec('isGas=isinstance('+ob+',gas)')
    if isGas:
        exec('gases.append('+ob+')')

#Update all the gases
for gas1 in gases:
    gas1.update()
#
#
#----------------Radiation related functions-------------

#Planck function (of frequency)
def B(nu,T):
    u = min(h*nu/(k*T),500.) #To prevent overflow
    return (2.*h*nu**3/c**2)/(math.exp(u)-1.)
#
#

