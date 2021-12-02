from __future__ import division, print_function, absolute_import
### THIS SCRIPT IS BASED ON PyTran, WHICH IS PART OF THE COURSEWARE IN
###  Pierrehumbert, 2010, Principles of Planetary Climate
###
### MODIFIED BY dkoll

#--------------------------------------------------------
#Description:
# [...]
#
#Note that a limitation of PyTran is that it uses a cutoff
#Lorentz line shape to synthesize absorption from line data.
#This is mainly to keep things simple and easy to understand,
#and it covers most cases of interest adequately well. However,
#a "professional strength" code should use the Voigt line shape
#instead, to allow for the dominance of Doppler broadening near
#line centers when pressure is low. In addition, to get the line
#overlap right in high pressure situations (including Early Mars)
#one ought to consider a more careful treatment of the non-Lorentz
#far tails of collisionally broadened lines.  The student who wishes
#to explore these effects (the latter of which is leading-edge research)
#will find it easy to modify the code to accomodate different assumptions
#on line shape.
#
#As currently written, Pytran loads in the lines for the dominant
#isotopologue for each molecule (based on abundance in Earth's
#atmosphere). If you want to modify the code to look at minor
#isotopologues, it is important to note that the HITRAN database
#downweights the line strengths for each isotopologue according
#to relative abundance in Earth's atmosphere.
#--------------------------------------------------------
#
#Change Log
#   3/13/2012: Corrected algebraic prefactor in temperature
#              scaling of line strength, and put in a more general
#              line-dependent scaling formula for the exponential factor
#              In this release, a generic power-law dependence for the
#              partition function Q(T) is used, but in the next release
#              I will implement the exact Q(T) for selected molecules,
#              based on routines provided as part of the HITRAN distribution.
#
#   2/22/2018: DKOLL - adapt PyTran to a newer database, like HITRAN2016
#
#   Dec 2021: DKOLL - clean up old functions;
#             replace math with numpy;
#             include Voigt line shape
#             a) based on canned routine/scipy Faddeeva fn: https://scipython.com/book/chapter-8-scipy/examples/the-voigt-profile/
#                -> this implementation is 3-4x slower than the Lorentz approx!
#
#             NOTES:
#             - potential alternatives for voigt implemtation, https://atmos.eoc.dlr.de/tools/lbl4IR/
#             - The 'relative' line cutoff option causes issues at low pressures!!
#               First is physical: kappa in line center blows up as p->0. Not a bug, in that it's consistent with Lorentz line approx.
#               Second is numerical:  see "nsum = int(numWidths*gam/dWave)"
#                      the line gets narrower & narrower, until it falls below the numerical wave spacing, numWidths*gam < dWave,
#                      at which point numpy.arange(i1-iw,i2-iw)=numpy.arange(iw,iw) produces an empty array.
#               So first kappa gets large, but once p-broadened linewidths drop below grid spacing kappa=0.
#             - The 'absolute' line cutoff just blows up as p->0, for lorentz lines...
#---------------------------------------------------------

#import string,math
import numpy as np
from .ClimateUtilities import *
from . import phys
import os

from scipy.special import wofz  # DKOLL -- for voigt profile: accurate but slower than lorentz


#Path to the datasets
datapath = '/'.join( os.path.abspath(__file__).split('/')[:-2] ) + '/DATA/HITRAN_DATA/'

#Path to the hitran by-molecule database
hitranPath = datapath+'HITRAN2016/ThermalOnly_0-5000cm.MainIsotopesOnly/'

#------------Constants and file data------------
#
#Hitran field locations
fieldLengths = [2,1,12,10,10,5,5,10,4,8,15,15,15,15,6,12,1,7,7]
Sum = 0
fieldStart = [0]
for length in fieldLengths:
    Sum += length
    fieldStart.append(Sum)
iso = 1
waveNum = 2
lineStrength = 3
airWidth =  5
selfWidth = 6
Elow = 7
TExp = 8
#
#
#Total internal partition functions (or generic approximations).
#These are used in the temperature scaling of line strength.
#The generic partition functions are OK up to about 500K
#(somewhat less for CO2)
def QGenericLin(T): #Approx. for linear molecules like CO2
    return T
def QGenericNonLin(T): #Approx for nonlinear molecules like H2O
    return T**1.5
#**ToDo: Provide actual partition functions for CO2, H2O and CH4

#Molecule numbers and molecular weights
#Add more entries here if you want to do other
#molecules in the HITRAN database. These entries are
#for the major isotopologues, but by using different
#molecule numbers you can do other isotopologues.
#The dictionary below maps the molecule name to the HITRAN
#molecule number (see documentation) and corresponding molecular
#weight.
#
#**ToDo:*Correct this to allow for minor isotopomers.
#       *Improve structure of the molecule dictionary,
#           e.g. use objects instead of arrays, allow for isotopologues
#       *Add in entries for the rest of the molecules
molecules = {} #Start an empty dictionary
molecules['H2O'] = [1,18.,QGenericNonLin] #Entry is [molecule number,mol wt,partition fn]
molecules['CO2'] = [2,44.,QGenericLin]
molecules['O3']  = [3,48.,QGenericNonLin]
molecules['N2O'] = [4,44.,QGenericLin]

molecules['CH4'] = [6,16.,QGenericNonLin]
molecules['NH3'] = [11,17.,QGenericNonLin] # linear structure?

molecules['HCN'] = [23,27.,QGenericNonLin]
molecules['C2H6'] = [27,30.,QGenericNonLin]

molecules['SF6'] = [30,146.,QGenericNonLin]  # careful: old file!

#-----------------------------------------------
# DKOLL: line shape functions
#        baed on https://scipython.com/book/chapter-8-scipy/examples/the-voigt-profile/

""" Return Gaussian line shape at x with HWHM alpha """
def lineshape_G(x, alpha):
    return np.sqrt(np.log(2) / np.pi)/alpha * np.exp(-(x/alpha)**2 * np.log(2))

""" Return Lorentzian line shape at x with HWHM gamma """
def lineshape_L(x, gamma):
    return gamma / (np.pi* (x**2 + gamma**2))

"""
Return the Voigt line shape at x with Lorentzian component HWHM gamma
and Gaussian component HWHM alpha. """
def lineshape_V(x, alpha, gamma):
    sigma = alpha / np.sqrt(2 * np.log(2))
    return np.real(wofz((x + 1j*gamma)/(sigma*np.sqrt(2)))) / (sigma*np.sqrt(2*np.pi))


#-----------------------------------------------
#Gets the fieldNum'th data item from a Hitran2004 record
def get(line,fieldNum):
    return line[fieldStart[fieldNum]:fieldStart[fieldNum]+fieldLengths[fieldNum]]


#Computes the absorption spectrum on a wave grid, by summing up
#contributions from each line.  numWidths is the number
#of line widths after which the line contribution is cut off.
#Typically we use 100-1000 for Earth troposphere, but in low pressure
#(like Mars or upper strat) values as high as 10000 might be needed.
#The validity of the Lorentz shape at such large cutoffs is dubious.
#At The cutoff can affect the absorption significantly in the
#water vapor or CO2 window, or "continuum" regions
def computeAbsorption(waveGrid,getGamma,p,T,dWave,numWidths = 1000.):
    N_grid = len(waveGrid)
    absGrid = numpy.zeros(N_grid,numpy.Float)

    # DKOLL ..
    alpha_factor = 1./phys.c * np.sqrt(phys.N_avogadro*phys.k*T*np.log(2.)/(molecules[molName][1]*1e-3)) # [unitless]; 1e-3 from g->kg
    
    for i in range(len(waveList)):
        n = waveList[i] # Wavenumber of the line
        #gam = gamList[i]*(p/1.013e5)*(296./T)**TExpList[i]  # DKOLL: old
        gam = getGamma(i)*(296./T)**TExpList[i]              # DKOLL: new. getGamma includes p-scaling
        #Temperature scaling of line strength
        Tfact = np.exp(-100.*(phys.h*phys.c/phys.k)*ElowList[i]*(1/T-1/296.))
        #The following factor is usually pretty close to unity
        #for lines that aren't far from the peak of the Planck spectrum
        #for temperature T, but it can become important on the low frequency
        #side, and is easy to incorporate.
        Tfact1 = (1.- np.exp(-100.*(phys.h*phys.c/phys.k)*n/T))/ \
              (1.- np.exp(-100.*(phys.h*phys.c/phys.k)*n/296.))
        #The following has the incorrect algebraic prefactor used in
        #  the original version of PyTran (see Errata/Updates document)
        #S = sList[i]*(T/296.)**TExpList[i]*Tfact
        #The following is the corrected version, including also the
        #  low frequency factor Tfact1
        #S = sList[i]*(Q(296.)/Q(T))*TExpList[i]*Tfact*Tfact1
        #Preceding line didn't delete "*TExpList" factor.  Results now
        #checked against LMD kspectrum code, for Lorentz line case
        #-->Corrected on 6/10/2013
        S = sList[i]*(Q(296.)/Q(T))*Tfact*Tfact1
        #
        iw = int(N_grid*(n-waveStart)/(waveEnd-waveStart))
        nsum = int(numWidths*gam/dWave)
        i1 = max(0,iw-nsum)
        i2 = min(N_grid-1,iw+nsum)
        if i2>0:
            dn = numpy.arange(i1-iw,i2-iw)*dWave
            #abs = S*gam/(np.pi*( dn**2 + gam**2))   # old

            ## New - lorentz only
            #abs = S*lineshape_L(dn,gam)
            ## New - doppler onlu
            #alpha = n*alpha_factor
            #abs = S*lineshape_G(dn,alpha)   # units of alpha=[n], so cm-1            
            ## New - voigt line
            alpha = n*alpha_factor
            abs = S*lineshape_V(dn,alpha,gam)   # units of alpha=[n], so cm-1            
            absGrid[i1:i2] += abs
    return absGrid


### DKOLL: add option to have a fixed cutoff.
###        i.e., truncate line at N cm^-1 away from center instead of N halfwidths
###        For example, MT_CKD continuum is defined as everything beyond 25cm^-1.
###
### DKOLL: also allow for option to remove the Lorenz line 'plinth',
##         cf. MTCKD continuum references
def computeAbsorption_fixedCutoff(waveGrid,getGamma,p,T,dWave,numWidths=25.,remove_plinth=False):
    N_grid = len(waveGrid)
    absGrid = numpy.zeros(N_grid,numpy.Float)

    # DKOLL ..
    alpha_factor = 1./phys.c * np.sqrt(phys.N_avogadro*phys.k*T*np.log(2.)/(molecules[molName][1]*1e-3)) # [unitless]; 1e-3 from g->kg

    for i in range(len(waveList)):
        n = waveList[i] # Wavenumber of the line
        gam = getGamma(i)*(296./T)**TExpList[i]              # DKOLL: new. getGamma includes p-scaling
        #Temperature scaling of line strength
        Tfact = np.exp(-100.*(phys.h*phys.c/phys.k)*ElowList[i]*(1/T-1/296.))
        #The following factor is usually pretty close to unity
        #for lines that aren't far from the peak of the Planck spectrum
        #for temperature T, but it can become important on the low frequency
        #side, and is easy to incorporate.
        Tfact1 = (1.- np.exp(-100.*(phys.h*phys.c/phys.k)*n/T))/ \
              (1.- np.exp(-100.*(phys.h*phys.c/phys.k)*n/296.))
        #The following has the incorrect algebraic prefactor used in
        #  the original version of PyTran (see Errata/Updates document)
        #S = sList[i]*(T/296.)**TExpList[i]*Tfact
        #The following is the corrected version, including also the
        #  low frequency factor Tfact1
        #S = sList[i]*(Q(296.)/Q(T))*TExpList[i]*Tfact*Tfact1
        #Preceding line didn't delete "*TExpList" factor.  Results now
        #checked against LMD kspectrum code, for Lorentz line case
        #-->Corrected on 6/10/2013
        S = sList[i]*(Q(296.)/Q(T))*Tfact*Tfact1
        #
        iw = int(N_grid*(n-waveStart)/(waveEnd-waveStart))
        #nsum = int(numWidths*gam/dWave)   # DKOLL: old
        nsum = int( numWidths/dWave )  # DKOLL: new
        i1 = max(0,iw-nsum)
        i2 = min(N_grid-1,iw+nsum)
        # DKOLL:
        if (i2>0) and (remove_plinth==False):
            dn = numpy.arange(i1-iw,i2-iw)*dWave
            #abs = S*gam/(np.pi*( dn**2 + gam**2))  # old

            ## New - lorentz
            #abs = S*lineshape_L(dn,gam)
            ## New - doppler
            #alpha = n*alpha_factor
            #abs = S*lineshape_G(dn,alpha)   # units of alpha=[n], so cm-1            
            ## New - voigt
            alpha = n*alpha_factor
            abs = S*lineshape_V(dn,alpha,gam)   # units of alpha=[n], so cm-1

            absGrid[i1:i2] += abs
        # DKOLL: option to remove plinth if using with MTCKD
        elif (i2>0) and (remove_plinth==True):
            dn = numpy.arange(i1-iw,i2-iw)*dWave
            #abs = S*gam/np.pi * (1./(dn**2 + gam**2) - 1./(numWidths**2 + gam**2) )   # old

            ## New - lorentz
            #abs = S*(lineshape_L(dn,gam) - lineshape_L(numWidths,gam))
            ## New - doppler
            #alpha = n*alpha_factor
            #abs = S*(lineshape_G(dn,alpha) - lineshape_G(numWidths,alpha))   # units of alpha=[n], so cm-1            
            ## New - voigt
            alpha = n*alpha_factor
            abs = S*(lineshape_V(dn,alpha,gam) - lineshape_V(numWidths,alpha,gam))   # units of alpha=[n], so cm-1
            
            absGrid[i1:i2] += abs
    return absGrid



#Open the Hitran file for the molecule in question,
#and read in the line data
#**ToDo: add isotope index to argument
#        Note that HITRAN weights the line strengths
#        according to relative abundance of each isotopologue
#        in Earth's atmosphere. When loading minor isotopologues,
#        it is important to undo this weighting, so that the
#        correct abundance weighting for the actual mixture of
#        interest can be computed.
#**DKOLL:
#        modified to speed up long line lists.
#        skip lines until you're close to your spectral region of interest.
#        then lead lines
#        if you're outside region, break, skip rest of file.
#        for this to work, make use of fact that line lists are ordered
def loadSpectralLines(molName,minWave=None,maxWave=None):
    global waveList,sList,gamAirList,gamSelfList,ElowList,TExpList,Q
    molNum = molecules[molName][0]
    Q = molecules[molName][2] #Partition function for this molecule
    if molNum < 10:
        file = hitranPath+'0%d_hit16.par'%molNum
    elif molNum == 30:
        # need an exception for SF6 because old available data is old!
        ##file = hitranPath+'%d_hit08.par'%molNum
        file = hitranPath+'SF6_data_reducedByRay.par'  # only using a subset of lines!
    else:
        file = hitranPath+'%d_hit16.par'%molNum
    f = open(file)
    waveList = []
    sList = []
    gamAirList = []
    gamSelfList = []
    TExpList = []
    ElowList = [] #Lower state energy
    #
    # DKOLL: only consider lines close to
    #        frequencies for which we actually want kappa.
    #        to be conservative, say 1000 cm-1.
    #        ->that might imply cutting off some far tails!
    if minWave is not None:
        minWave = minWave - 1000.
    else:
        minWave = 0.
    if maxWave is not None:
        maxWave = maxWave + 1000.
    else:
        maxWave = 1e10 # something big
    #
    count = 0
    line = f.readline()
    while len(line)>0:
        count += 1
        isoIndex = int(get(line,iso))
        n = float(get(line,waveNum))
        #DKOLL:
        if n<minWave:
            line = f.readline()   # skip to next line
        elif n>maxWave:
            break  # ignore rest of file
        else:
            S = float(get(line,lineStrength))
            El = float(get(line,Elow))
            #Convert line strength to (m**2/kg)(cm**-1) units
            #The cm**-1 unit is there because we still use cm**-1
            #as the unit of wavenumber, in accord with standard
            #practice for IR
            #
            #**ToDo: Put in correct molecular weight for the
            #        isotope in question.
            S = .1*(phys.N_avogadro/molecules[molName][1])*S
            gamAir = float(get(line,airWidth))
            gamSelf = float(get(line,selfWidth))
            TemperatureExponent = float(get(line,TExp))
            if  isoIndex == 1:
                waveList.append(n)
                sList.append(S)
                gamAirList.append(gamAir)
                gamSelfList.append(gamSelf)
                ElowList.append(El)
                TExpList.append(TemperatureExponent)
            #Read the next line, if there is one
            line = f.readline()
    #Convert the lists to numpy/Numeric arrays
    waveList = numpy.array(waveList)
    sList = numpy.array(sList)
    gamAirList = numpy.array(gamAirList)
    gamSelfList = numpy.array(gamSelfList)
    ElowList =  numpy.array(ElowList)
    TExpList = numpy.array(TExpList)
    f.close()


#====================================================================
#
#                End of function definitions
#                   Main script starts here
#====================================================================


#Program data and initialization


#Standard wavenumbers used for spectral survey

def getKappa_HITRAN(waveGrid,wave0,wave1,delta_wave,molecule_name,\
                        press=1e4,temp=300.,lineWid=1000.,broadening="mixed", \
                        press_self=None, \
                        cutoff_option="relative",remove_plinth=False):

    # ! not a great solution !
    global waveStart, waveEnd, dWave, p, T, molName
    waveStart = wave0
    waveEnd = wave1
    dWave = delta_wave
    molName = molecule_name

    p = float(press) # DKOLL: make sure (p,T) are floats!
    T = float(temp)

    loadSpectralLines(molName,minWave=wave0,maxWave=wave1)

#-->Decide whether you want to compute air-broadened
#or self-broadened absorption. The plots of self/air
#ratio in the book were done by running this script for
#each choice and computing the ratio of the absorption coefficients
    p_tot  = press/1.013e5              # need [atm]!
    p_self = press_self/1.013e5
    p_air  = (press-press_self)/1.013e5
    if broadening=="air":
        getGamma = lambda i: gamAirList[i]*p_tot
    elif broadening=="self":
        getGamma = lambda i: gamSelfList[i]*p_tot
    elif broadening=="mixed":
        # see HITRAN documentation.
        #   e.g., http://hitran.org/docs/definitions-and-units/
        if press_self==None:
            press_self = press
        getGamma = lambda i: gamAirList[i]*p_air + gamSelfList[i]*p_self

#-->Compute the absorption on the wavenumber grid
#Set nWidths to the number of line widths to go out to in
#superposing spectral lines to compute the absorption coefficient.
#This is a very crude implementation of the far-tail cutoff.
#I have been using 1000 widths as my standard value.
    nWidths = lineWid

    if cutoff_option=="relative":
        absGrid = computeAbsorption(waveGrid,getGamma,p,T,dWave,nWidths)
    elif cutoff_option=="fixed":
        absGrid = computeAbsorption_fixedCutoff(waveGrid,getGamma,p,T,dWave,nWidths,remove_plinth=remove_plinth)

    return absGrid
