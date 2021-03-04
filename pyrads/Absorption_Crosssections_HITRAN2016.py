from __future__ import division, print_function, absolute_import
### THIS SCRIPT IS BASED ON PyTran, WHICH IS PART OF THE COURSEWARE IN
###  Pierrehumbert, 2010, Principles of Planetary Climate
###
### MODIFIED BY Daniel Koll and Andrew Williams

#--------------------------------------------------------
#Description:
#Utilities for reading the HITRAN line database, synthesizing
#absorption spectrum on a regular wavenumber grid.  
# It also can be used in simple
#line-by-line transmission function calculations
#
#The functions of interest to most users are:
#      loadSpectralLines(molName)
#             Loads spectral line data for molecule molName.
#
#      computeAbsorption(p,T,dWave,LineTailCutoff)
#             Synthesizes the absorption coefficient vs. wavenumber
#             from the line data, using an assumed pressure-broadened
#             Lorentz line shape with tails cut off LineTailCutoff
#             cm**-1 from the line centers.
#
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
#
#---------------------------------------------------------

import string,math
from .ClimateUtilities import *
from . import phys
import os

#Path to the datasets
datapath = '/'.join( os.path.abspath(__file__).split('/')[:-2] ) + '/DATA/HITRAN_DATA/'

#Path to the hitran by-molecule database
hitranPath = datapath+'HITRAN2016/ThermalOnly_0-5000cm.MainIsotopesOnly/'

#------------Constants and file data------------
#
#Hitran field locations (for reading par files)
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

"""
Molecule numbers and molecular weights
Add more entries here if you want to do other
molecules in the HITRAN database. These entries are
for the major isotopologues, but by using different
molecule numbers you can do other isotopologues.
The dictionary below maps the molecule name to the HITRAN
molecule number (see documentation) and 
corresponding molecular weight.
"""
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


def Planck(wavenum,T):
    """ Planck Function """
    return 100.*math.pi*phys.c*phys.B(100.*wavenum*phys.c,T)

def computeAbsorption(waveGrid,waveStart,waveEnd,
                      hitran_data,getGamma,
                      p,T,dWave,
                      numWidths = 1000.):
    """
    Computes the absorption spectrum on a wave grid, by summing up
    contributions from each line.  numWidths is the number
    of line widths after which the line contribution is cut off.
    Typically we use 100-1000 for Earth troposphere, but in low pressure
    (like Mars or upper strat) values as high as 10000 might be needed.
    The validity of the Lorentz shape at such large cutoffs is dubious.
    At The cutoff can affect the absorption significantly in the
    water vapor or CO2 window, or "continuum" regions
    """
    absGrid = numpy.zeros(len(waveGrid))
    waveList,sList,gamAirList,gamSelfList,ElowList,TExpList,Q = hitran_data
    
    for i in range(len(waveList)):
        n = waveList[i] # Wavenumber of the line
        gam = getGamma(i)*(296./T)**TExpList[i]              
        
        #Temperature scaling of line strength
        Tfact = math.exp(-100.*(phys.h*phys.c/phys.k)*ElowList[i]*(1/T-1/296.))

        Tfact1 = (1.- math.exp(-100.*(phys.h*phys.c/phys.k)*n/T))/ \
              (1.- math.exp(-100.*(phys.h*phys.c/phys.k)*n/296.))

        S = sList[i]*(Q(296.)/Q(T))*Tfact*Tfact1
        
        iw = int(len(waveGrid)*(n-waveStart)/(waveEnd-waveStart))
        nsum = int(numWidths*gam/dWave)
        i1 = max(0,iw-nsum)
        i2 = min(len(waveGrid)-1,iw+nsum)
        if i2>0:
            dn = numpy.arange(i1-iw,i2-iw)*dWave
            abs = S*gam/(math.pi*( dn**2 + gam**2))
            absGrid[i1:i2] += abs
    return absGrid


def computeAbsorption_fixedCutoff(waveGrid,waveStart,waveEnd,
                                  hitran_data,getGamma,
                                  p,T,dWave,
                                  numWidths=25.,remove_plinth=False):
    """
    DKOLL: add option to have a fixed cutoff.
        i.e., truncate line at N cm^-1 away from center instead of N halfwidths
        For example, MT_CKD continuum is defined as everything beyond 25cm^-1.

     DKOLL: also allow for option to remove the Lorenz line 'plinth',
        cf. MTCKD continuum references
    """
    absGrid = numpy.zeros(len(waveGrid))
    waveList,sList,gamAirList,gamSelfList,ElowList,TExpList,Q = hitran_data

    for i in range(len(waveList)):
        n = waveList[i] # Wavenumber of the line
        gam = getGamma(i)*(296./T)**TExpList[i]              
        
        #Temperature scaling of line strength
        Tfact = math.exp(-100.*(phys.h*phys.c/phys.k)*ElowList[i]*(1/T-1/296.))
        Tfact1 = (1.- math.exp(-100.*(phys.h*phys.c/phys.k)*n/T))/ \
              (1.- math.exp(-100.*(phys.h*phys.c/phys.k)*n/296.))
        
        S = sList[i]*(Q(296.)/Q(T))*Tfact*Tfact1
        
        iw = int(len(waveGrid)*(n-waveStart)/(waveEnd-waveStart))
        nsum = int( numWidths/dWave )  
        i1 = max(0,iw-nsum)
        i2 = min(len(waveGrid)-1,iw+nsum)
        # DKOLL: old
        if (i2>0) and (remove_plinth==False):
            dn = numpy.arange(i1-iw,i2-iw)*dWave
            abs = S*gam/(math.pi*( dn**2 + gam**2))
            absGrid[i1:i2] += abs
        # DKOLL: new
        elif (i2>0) and (remove_plinth==True):
            dn = numpy.arange(i1-iw,i2-iw)*dWave
            abs = S*gam/math.pi * (1./(dn**2 + gam**2) - 1./(numWidths**2 + gam**2) )
            absGrid[i1:i2] += abs
    return absGrid


def get(line,fieldNum):
    """
    Gets the fieldNum'th data item from a Hitran2016 record
    """
    return line[fieldStart[fieldNum]:fieldStart[fieldNum]+fieldLengths[fieldNum]]


def loadSpectralLines(molName,minWave=None,maxWave=None):
    """
    Open the Hitran file for the molecule in question,
    and read in the line data
    
    ToDo: 
    Add isotope index to argument
    Note that HITRAN weights the line strengths
    according to relative abundance of each isotopologue
    in Earth's atmosphere. When loading minor isotopologues,
    it is important to undo this weighting, so that the
    correct abundance weighting for the actual mixture of
    interest can be computed.
    
    DKOLL:
    modified to speed up long line lists.
    skip lines until you're close to your spectral region of interest.
    then lead lines
    if you're outside region, break, skip rest of file.
    for this to work, make use of fact that line lists are ordered
    """
    
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
    
    hitran_lists = [waveList,sList,gamAirList,gamSelfList,ElowList,TExpList,Q]
    return hitran_lists


#====================================================================
#                   Main script starts here
#====================================================================


def getKappa_HITRAN(waveGrid,wave0,wave1,delta_wave, \
                    press=1e4,temp=300.,lineWid=1000.,broadening="mixed", \
                    press_self=None, molecule_name=None, hitran_data=None, \
                    cutoff_option="relative",remove_plinth=False):

    waveStart = wave0
    waveEnd = wave1
    dWave = delta_wave
    molName = molecule_name

    p = float(press) # DKOLL: make sure (p,T) are floats!
    T = float(temp)

    # Only run loadSpectralLines if necessary
    if (molName is not None) and (hitran_data is None):
        print("Running loadSpectralLines...")
        hitran_data = loadSpectralLines(molName,minWave=wave0,maxWave=wave1)
    
    #-->Decide whether you want to compute air-broadened
    #or self-broadened absorption. The plots of self/air
    #ratio in the book were done by running this script for
    #each choice and computing the ratio of the absorption coefficients
    gamAirList,gamSelfList = hitran_data[2], hitran_data[3]
    
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
        absGrid = computeAbsorption(waveGrid,waveStart,waveEnd,hitran_data,getGamma,p,T,dWave,nWidths)
    elif cutoff_option=="fixed":
        absGrid = computeAbsorption_fixedCutoff(waveGrid,waveStart,waveEnd,hitran_data,getGamma,p,T,dWave,nWidths,remove_plinth=remove_plinth)

    return absGrid
