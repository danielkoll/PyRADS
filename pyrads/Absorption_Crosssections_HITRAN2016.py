from __future__ import division, print_function, absolute_import
### THIS SCRIPT IS BASED ON PyTran, WHICH IS PART OF THE COURSEWARE IN
###  Pierrehumbert, 2010, Principles of Planetary Climate
###
### MODIFIED BY dkoll

#--------------------------------------------------------
#Description:
#Utilities for reading the HITRAN line database, synthesizing
#absorption spectrum on a regular wavenumber grid, and computing
#exponential sum coefficients.  It also can be used in simple
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
#      plotLogQuartiles(bandWidth)
#             Makes the plots of absorption statistics in
#             bands of width bandWidth, as discussed in
#             the greenhouse gas spectral survey in Chapter 4.
#
#      makeEsumTable(waveStart,waveEnd,Delta=50.,N=21)
#             Constructs the exponential sum tables used in
#             the homebrew radiation code.
#
#
#Scroll down to "Main script starts here" to see a basic example
#of use of the routines to plot absorption coefficients and make
#the plot of the quartile summaries in the spectral survey in the book.
#
#
#There are a few other things in PyTran for computing line-by-line
#band averaged transmission functions, and transmission along
#a path computed exactly at a fixed wavenumber without making
#pressure scaling assumptions. These are used in Chapter 4 to
#do some of the LBL transmission function calculations, but are
#probably of less general interest.
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
# #hitranPath = datapath+'HITRAN2016/ThermalOnly_0-5000cm.MainIsotopesOnly/'
hitranPath = datapath+'HITRAN2016/Full_0-60000cm.MainIsotopesOnly/'    # FOR SW calculations!

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


#Planck function
def Planck(wavenum,T):
    return 100.*math.pi*phys.c*phys.B(100.*wavenum*phys.c,T)
#
#Function to make an exponential sum table
#
def makeEsumTable(waveStart,waveEnd,Delta=50.,N=21):
	wave = waveStart
	c = Curve()
	while (wave+Delta ) <= waveEnd:
		header1 = 'logKappa.%d.%d'%(wave,wave+Delta)
		header2 = 'dH.%d.%d'%(wave,wave+Delta)
		bins,hist,cum = logHisto([wave,wave+Delta],N)
		c.addCurve(bins,header1)
		c.addCurve(hist,header2)
		wave += Delta
	return c

#Note: Log quartile values are exponentiated for plotting
def plotLogQuartiles(bandWidth):
    waveList = []
    statsList = []
    wave = waveStart
    while wave + bandWidth <= waveEnd:
        waveList.append(wave+bandWidth/2.)
        statsList.append(logQuartiles([wave,wave+bandWidth])[0])
        wave = wave+bandWidth
    c = Curve()
    c.addCurve(waveList,'Wavenumber')
    c.addCurve([math.exp(stat[0]) for stat in statsList],'Min')
    c.addCurve([math.exp(stat[1]) for stat in statsList],'q25')
    c.addCurve([math.exp(stat[2]) for stat in statsList],'Median')
    c.addCurve([math.exp(stat[3]) for stat in statsList],'q75')
    c.addCurve([math.exp(stat[4]) for stat in statsList],'Max')
    return c



def logQuartiles(waveRange):
    i1 = numpy.searchsorted(waveGrid,waveRange[0])
    i2 = numpy.searchsorted(waveGrid,waveRange[1])
    logAbs = [math.log(kappa) for kappa in absGrid[i1:i2] if kappa>0.]
    numZeros = (i2-i1) - len(logAbs)
    return quartiles(logAbs),numZeros

def quartiles(data):
    vs = numpy.sort(data)
    N = len(data)
    if N > 0:
        return [vs[0],vs[N/4],vs[N/2],vs[3*N/4],vs[-1]]
    else:
        return [-1.e20,-1.e20,-1.e20,-1.e20,-1.e20]


#Utility to compute histogram of absorption data
#Returns bins, histogram and cumulative PDF
def histo(data,nBins=10,binMin=None,binMax=None):
    #Flag bands with no lines in them
    if len(data) == 0:
        binMid = numpy.zeros(nBins-1,numpy.Float)
        hist = numpy.zeros(nBins-1,numpy.Float)
        hist[0] = 1.
        binMid[0] = -1.e30
        cumHist = numpy.zeros(nBins-1,numpy.Float)+1.
        return binMid,hist,cumHist
    #
    if binMin == None:
        binMin = min(data)
    if binMax == None:
        binMax = max(data)
    d = (binMax-binMin)/(nBins-1)
    bins = [binMin + d*i for i in range(nBins)]
    histo = [0. for i in range(nBins)]
    for value in data:
        i = int((value-binMin)/d - 1.e-5)
        if (i>=0) and (i<nBins):
            histo[i] += 1.
    n = sum(histo)
    histo = [h/n for h in histo]
    cum = 0.
    cumHisto = []
    for h in histo:
        cum += h
        cumHisto.append(cum)
    #
    bins = numpy.array(bins)
    histo = numpy.array(histo)
    cumHisto = numpy.array(cumHisto)
    binMid = .5*(bins[:-1]+bins[1:])
    return binMid,histo[0:-1],cumHisto[:-1]

#ToDo: Should count and report instances of zero kappa
def logHisto(waveRange,nBins = 10):
    i1 = numpy.searchsorted(waveGrid,waveRange[0])
    i2 = numpy.searchsorted(waveGrid,waveRange[1])
    logAbs = [math.log(kappa) for kappa in absGrid[i1:i2] if kappa>0.]
    return histo(logAbs,nBins)

def plotLogHisto(waveRange,data,nBins=10):
    bins,hist,cum = logHisto(waveRange,data,nBins)
    c = Curve()
    c.addCurve(bins)
    c.addCurve(hist)
    return c

#Computes the absorption spectrum on a wave grid, by summing up
#contributions from each line.  numWidths is the number
#of line widths after which the line contribution is cut off.
#Typically we use 100-1000 for Earth troposphere, but in low pressure
#(like Mars or upper strat) values as high as 10000 might be needed.
#The validity of the Lorentz shape at such large cutoffs is dubious.
#At The cutoff can affect the absorption significantly in the
#water vapor or CO2 window, or "continuum" regions
def computeAbsorption(waveGrid,getGamma,p,T,dWave,numWidths = 1000.):
    absGrid = numpy.zeros(len(waveGrid),numpy.Float)

    for i in range(len(waveList)):
        n = waveList[i] # Wavenumber of the line
        #gam = gamList[i]*(p/1.013e5)*(296./T)**TExpList[i]  # DKOLL: old
        gam = getGamma(i)*(296./T)**TExpList[i]              # DKOLL: new. getGamma includes p-scaling
        #Temperature scaling of line strength
        Tfact = math.exp(-100.*(phys.h*phys.c/phys.k)*ElowList[i]*(1/T-1/296.))
        #The following factor is usually pretty close to unity
        #for lines that aren't far from the peak of the Planck spectrum
        #for temperature T, but it can become important on the low frequency
        #side, and is easy to incorporate.
        Tfact1 = (1.- math.exp(-100.*(phys.h*phys.c/phys.k)*n/T))/ \
              (1.- math.exp(-100.*(phys.h*phys.c/phys.k)*n/296.))
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
        iw = int(len(waveGrid)*(n-waveStart)/(waveEnd-waveStart))
        nsum = int(numWidths*gam/dWave)
        i1 = max(0,iw-nsum)
        i2 = min(len(waveGrid)-1,iw+nsum)
        if i2>0:
            dn = numpy.arange(i1-iw,i2-iw)*dWave
            abs = S*gam/(math.pi*( dn**2 + gam**2))
            absGrid[i1:i2] += abs
    return absGrid


### DKOLL: add option to have a fixed cutoff.
###        i.e., truncate line at N cm^-1 away from center instead of N halfwidths
###        For example, MT_CKD continuum is defined as everything beyond 25cm^-1.
###
### DKOLL: also allow for option to remove the Lorenz line 'plinth',
##         cf. MTCKD continuum references
def computeAbsorption_fixedCutoff(waveGrid,getGamma,p,T,dWave,numWidths=25.,remove_plinth=False):
    absGrid = numpy.zeros(len(waveGrid),numpy.Float)

    for i in range(len(waveList)):
        n = waveList[i] # Wavenumber of the line
        #gam = gamList[i]*(p/1.013e5)*(296./T)**TExpList[i]  # DKOLL: old
        gam = getGamma(i)*(296./T)**TExpList[i]              # DKOLL: new. getGamma includes p-scaling
        #Temperature scaling of line strength
        Tfact = math.exp(-100.*(phys.h*phys.c/phys.k)*ElowList[i]*(1/T-1/296.))
        #The following factor is usually pretty close to unity
        #for lines that aren't far from the peak of the Planck spectrum
        #for temperature T, but it can become important on the low frequency
        #side, and is easy to incorporate.
        Tfact1 = (1.- math.exp(-100.*(phys.h*phys.c/phys.k)*n/T))/ \
              (1.- math.exp(-100.*(phys.h*phys.c/phys.k)*n/296.))
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
        iw = int(len(waveGrid)*(n-waveStart)/(waveEnd-waveStart))
        #nsum = int(numWidths*gam/dWave)   # DKOLL: old
        nsum = int( numWidths/dWave )  # DKOLL: new
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


# #Function to  to compute absorption coefficient
# #vs. pressure at a given wavenumber wavenum. and temperature T. Line parameter data
# #is global and must be read in first.
# def absPath(p,wavenum,T=296.,numWidths = 1000.):
#         width = .1*p/1.e5
#         i1 = numpy.searchsorted(waveList,wavenum-numWidths*width)-1
#         i2 = numpy.searchsorted(waveList,wavenum+numWidths*width)+1
#         i1 = max(i1,0)
#         i2 = min(i2,len(waveList))
#         lineCenters = waveList[i1:i2]
#         lineStrengths = sList[i1:i2]
#         lineWidths = gamAirList[i1:i2]  #gamList[i1:i2]
#         TExpList1 = TExpList[i1:i2]
#         ElowList1 = ElowList[i1:i2]
#         gam = lineWidths*(p/1.013e5)*(296./T)**TExpList1
#         #Temperature scaling of line strength
#         Tfact = numpy.exp(-100.*(phys.h*phys.c/phys.k)*ElowList1*(1/T-1/296.))
#         S = lineStrengths*(T/296.)**TExpList1*Tfact
#         terms = S*gam/ \
#           (math.pi*( (lineCenters-wavenum)**2 + gam**2))
#         return sum(terms)

#Function to compute a list of values of the transmission
#at a given frequency, for each of a set of pressures.
#This can be used to compute the band averaged transmission
#in a given band between wave1 and wave2 by calling transPath
#for a regular list of wavenumbers covering the interval, and summing
#the results. For the transmission computed in the book, we actually
#did this integration in a Monte-Carlo fashion, by picking random
#wavenumbers in the interval and accumulating the results, This
#has the advantage that you can watch the progress of the computation
#as it converges, and you get some information at once from the
#entire band.  For a typical band of width 50 cm**-1, even 1000
#samples can give quite good results
#
#**ToDo:  Modify this to take a list of pressure and temperature
#instead. Also, fix the generation of the pressure
#list, so it ends at p2 instead of running over. Perhaps also
#allow for changing T and numWidths in the absorption computation
def transPath(wave,p1,p2,q,ndiv):
	def f(p,q):
		return absPath(p,wave)*q/g
	dtau1 = romberg(f)
	dtau = 0.
	dp = (p2-p1)/ndiv
	p = p1
	pList = [p1]
	transList = [1.]
	trans = 1.
	while p <= p2:
		dtau += dtau1([p,p+dp],q)
		p += dp
		pList.append(p)
		transList.append(math.exp(-abs(dtau)))
	return numpy.array(pList),numpy.array(transList)

#Returns a curve for making a plot of absorption vs. pressure
#at a given frequency
def plotP(n,T=296.,numWidths = 1000.):
	press = [10.**(.25*i) for i in range(35)]
	absp = [absPath(p,n,T,numWidths) for p in press]
	c = Curve()
	c.addCurve(press)
	c.addCurve(absp)
	c.YlogAxis = c.XlogAxis = True
	return c

#Gets the fieldNum'th data item from a Hitran2004 record
def get(line,fieldNum):
    return line[fieldStart[fieldNum]:fieldStart[fieldNum]+fieldLengths[fieldNum]]

#Computes the mean transmission between two specified
#wavenumbers, for a given absorber path. This
#version does not take into account pressure broadening.
#It should be modified to do so.
def TransBar(wave1,wave2,massPath):
    i1 = numpy.searchsorted(waveGrid,wave1)
    i2 = numpy.searchsorted(waveGrid,wave2)
    trans = numpy.exp(-massPath*absGrid[i1:i2])
    return numpy.average(trans)
def plotTransBar(wave1,wave2,massPathStart,massPathEnd):
    nplot = 100
    r = math.log(massPathEnd/massPathStart)/nplot
    mPath = [massPathStart*math.exp(r*i) for i in range(nplot)]
    trans = [TransBar(wave1,wave2,path) for path in mPath]
    c = Curve()
    c.addCurve(mPath)
    c.addCurve(trans)
    return c

#Function to compute OLR spectrum line-by-line This
#version is to demonstrate the concept. It doesn't take
#into account the pressure dependence of the absorption
def Tprof(pps):
    return max(280.*pps**(2./7.),200.)
def OLRspec(pathMax):
    ps = 1.e5
    dp = .05*ps
    transLast = 1.+ numpy.zeros(len(absGrid),numpy.Float)
    OLR = numpy.zeros(len(absGrid),numpy.Float)
    pp = 0.
    while pp < ps:
        pp += dp
        path = pathMax*pp/ps
        trans = numpy.exp(-path*absGrid)
        TT = (Tprof(pp/ps)+ Tprof((pp-dp)/ps))/2.
        B = numpy.array([Planck(w,TT) for w in waveGrid])
        OLR += B*(transLast-trans)
        transLast = trans
    TT = Tprof(1.)
    B = numpy.array([Planck(w,TT) for w in waveGrid])
    OLR += B*trans
    return OLR

#Computes a function smoothed over specified wavenumber band
def smooth(wAvg,data):
    navg = int(wAvg/dWave)
    nmax = len(waveGrid)-navg
    return [numpy.average(data[i:i+navg]) for i in range(0,nmax,navg)]



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
        #isoIndex = string.atoi(get(line,iso))
        isoIndex = int(get(line,iso))
        #n = string.atof(get(line,waveNum))
        n = float(get(line,waveNum))
        #DKOLL:
        if n<minWave:
            line = f.readline()   # skip to next line
        elif n>maxWave:
            break  # ignore rest of file
        else:
            #S = string.atof(get(line,lineStrength))
            S = float(get(line,lineStrength))
            #El = string.atof(get(line,Elow))
            El = float(get(line,Elow))
            #Convert line strength to (m**2/kg)(cm**-1) units
            #The cm**-1 unit is there because we still use cm**-1
            #as the unit of wavenumber, in accord with standard
            #practice for IR
            #
            #**ToDo: Put in correct molecular weight for the
            #        isotope in question.
            S = .1*(phys.N_avogadro/molecules[molName][1])*S
            #gamAir = string.atof(get(line,airWidth))
            gamAir = float(get(line,airWidth))
            #gamSelf = string.atof(get(line,selfWidth))
            gamSelf = float(get(line,selfWidth))
            #TemperatureExponent = string.atof(get(line,TExp))
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
                        press=1e4,temp=300.,lineWid=1000.,broadening="air", \
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
