### -------------------------
### Utilities for reading in HITRAN CIA data.
###
### note that the files are arranged with
### "chunks"//multiple T sets and wave bands, stacked
### within a single file!
### 
### DKOLL 11.10.2018
###
### COMMENTS:
### - careful about how to extrapolate data beyond the given bounds!
### - missing data = zeros! deal with zeros appropriately.
### - it's very inefficient to read in a whole file.
###    solve by reading through only once and extracting header lines
###    then can go back to specific headers.

from __future__ import division, print_function, absolute_import
from . import phys
import numpy as np
import os

from scipy.interpolate import interp1d

### -------------------------
### Global variables and helpers

datapath = '/'.join( os.path.abspath(__file__).split('/')[:-2] ) + '/DATA/'  # !

cia_path = datapath + "HITRAN_DATA/CIA/"

# Below: dict with available continua. Add as needed.
cia_dict = {"N2-H2" : "N2-H2_2011.cia", \
            #"H2-H2" : "H2-H2_2011.cia", \
            "H2-H2" : "H2-H2_2011.T_less_800K_only.cia", \
            # "H2-He" : "H2-He_2011.cia", \
            "N2-N2" : "N2-N2_2011.0-554cm.cia", \
            "N2-N2_2000cm_cold" : "N2-N2_2011.2000-2700cm.cia", \
            "N2-N2_2000cm_hot" : "N2-N2_2011.1850-3000cm.cia", \
            }  # map molec to filename

class Dummy():
    pass


### -------------------------
### Parse cia header
### see: http://hitran.org/data/CIA/CIA_Readme.pdf
### 
### note: I'm splitting files at "CO2-CO2", so only use rest of header past that point...
def parse_cia_header(header):

    nu0 = float( header[0:10].lstrip() )
    nu1 = float( header[10:20].lstrip() )
    N_spec = float( header[20:27].lstrip() )
    T = float( header[27:34].lstrip() )

    return nu0,nu1,T


### -------------------------
### Read through a cia file to isolate header lines.
### Return lists of (nu0,nu1,temp,line nr in datafile)
def get_header_lines(file,cia):
    wave0,wave1,temps,lines = [],[],[],[]
    with open(file) as f:
        N = 0
        for line in f:
            if cia in line:
                # (need to remove leading whitespace and cia name)
                nu0,nu1,T = parse_cia_header( line.lstrip()[len(cia):] )
                
                wave0.append(nu0)
                wave1.append(nu1)
                temps.append(T)
                lines.append(N)
            N+=1

    wave0 = np.array(wave0)  # make sure to return np arrays
    wave1 = np.array(wave1)
    temps = np.array(temps)
    lines = np.array(lines)

    return wave0,wave1,temps,lines

### -------------------------
### Input: range of lines between (n0,n1)
###        including line nr n0 but not line nr n1
### Output: data for those lines
def get_data(file,n0,n1):
    with open(file) as f:
        N = 0
        data = []
        for line in f:
            if N<n0:
                pass
            elif N>=n1:
                break
            else:
                data.append(line)
            N+=1

    return data

### -------------------------
### Returns one or two hitran text chunks, depending on temp
def read_cia(file,cia,T):

    # first, read in the header lines
    # CAREFUL - list should be sorted by temp = cia file should be sorted!
    wave0,wave1,temps,lines = get_header_lines(file,cia)

    # second, find the closest 1 or 2 headers
    if T<=min(temps):
        n0 = lines[0]
        n1 = lines[1]
        chunk0 = get_data(file,n0,n1)
        chunk1 = None

    elif T>=max(temps):
        n0 = lines[-1]
        n1 = np.inf
        chunk0 = get_data(file,n0,n1)
        chunk1 = None

    else:
        n0 = lines[temps<T][-1]
        n1 = lines[temps>=T][0]
        try:
            n2 = lines[temps>=T][1]
        except:
            n2 = np.inf

        chunk0 = get_data(file,n0,n1)
        chunk1 = get_data(file,n1,n2)

    return chunk0,chunk1

### -------------------------
### Input: chunk of textfile
### Output: object with (wavenr,cia coef)
def parse_chunk(chunk,cia,zero_bad_data=True):
    c = Dummy()

    nu0,nu1,T = parse_cia_header( chunk[0].lstrip()[len(cia):] ) # ?
    c.T = T
    
    data = np.genfromtxt( chunk, skip_header=1 )
    c.nu = data[:,0]
    c.b = data[:,1]

    # only return non-zero,non-NaN data?
    if zero_bad_data:
        mask = np.logical_or( c.b < 0., np.isnan(c.b) )  # all *bad* values
            
        nu_subset = c.nu[~mask].copy()
        b_subset = c.b[~mask].copy()

        c.nu = nu_subset
        c.b = b_subset

    return c


### -------------------------------------------------------
###       MAIN FUNCTION FOR OUTSIDE ACCESS:
###
### Get kappa at desired (nu,T) by interpolating HITRAN CIA data.
###   The binary absorption cross-sec is given in cm^5/molec^2, so need to multiply by number
###   density of both colliding molecs. e.g:
###     tau_H2-He   = b*n_H2*n_He*dz
###
###   To make this compatible with mass-specific cross-sec in rest of my code, 
###   in terms of partial pressures, use
###     tau_A-B   = b*(p_A/p_tot)*(p_B/p_tot)*(R_A/R_mean)*(R_B/R_mean)*1/(M_A*M_B)*rho_tot * dp/g
###
###   in terms of specific masses, use
###     tau_A-B   = b*(q_A/M_A)*(q_B/M_B)*rho_tot * dp/g
###
### WATCH OUT:
###    This function returns 'kappa*q', To get tau, just multiply *dp/g,
###    not *q dp/g!
###
### USAGE:
###    get_kappa_cia("CO2","CO2",[1000.],300.,qA,qB,rho_tot)
###    get_kappa_cia("CO2","CO2",[10.,15.,20.,...],300.,qA,qB,rho_tot)
###
### CAREFUL:
###    - interpolate linearly, or log-linearly?
###       --> makes a difference bc it shapes the tails of the CIA features!
###    - my current interpolation does *NOT* work for datasets with overlapping values,
###      e.g., the full N2-N2 or O2-O2 datasets!
###      --> simple fix, manually take datasets apart
###    - also, the fn spends a lot of time trying to interpolate where I set things to zero anyway.

def get_kappa_cia(molecA,molecB,nu,T,qA,qB,rho,interpolate_linearly=True):
    gasA = getattr(phys,molecA)
    gasB = getattr(phys,molecB)

    # make sure 'wavenr' is a np array
    nu = np.asarray( nu )

    # get data:
    continuum = molecA + "-" + molecB   # e.g, "H2" + "He"

    c0,c1 = read_cia( cia_path + cia_dict[continuum],continuum,T )

    dat0 = parse_chunk(c0,continuum)
    if c1 is not None:
        dat1 = parse_chunk(c1,continuum)

        
    # interpolate first to desired frequencies at each temp, then interpolate in temp
    # FOR TEMP: do interpolation always linearly!
    # FOR FREQUECY:  'interpolate_linearly' is probably better justified?
    #                -> assume zero beyond bounds!
    if interpolate_linearly:
        interp_fn = lambda x,x_data,y_data: np.interp(x,x_data,y_data,left=0.,right=0.)
        #interp_fn = lambda x,x_data,y_data: interp1d(x_data,y_data,fill_value=0.,bounds_error=False,assume_sorted=True)(x)
    else:
        interp_fn = lambda x,x_data,y_data: np.exp(np.interp(x,x_data,np.log(y_data)))
        #interp_fn = lambda x,x_data,y_data: np.exp( interp1d(x_data, np.log(y_data),fill_value='extrapolate',bounds_error=False,assume_sorted=True)(x) )

    # ---
    # interpolate in frequency
    # here: don't try to interpolate beyond data bounds
    b0 = np.zeros_like( nu )
    mask = np.logical_and( dat0.nu.min()<=nu, nu<=dat0.nu.max() )
    b0[mask] = interp_fn( nu[mask], dat0.nu, dat0.b )
    if c1 is not None:
        b1 = np.zeros_like( nu )
        mask = np.logical_and( dat1.nu.min()<=nu, nu<=dat1.nu.max() )
        b1[mask] = interp_fn( nu[mask], dat1.nu, dat1.b )

    # ---
    # interpolate in temp, if needed
    # linear interp: b(T)=b0 + db/dT*(T-T0)
    if c1 is not None:
        b = b0 + (b1-b0)*(T-dat0.T)/(dat1.T-dat0.T)
    else:
        b = b0

    # g/mole -> kg/mole
    M_A = gasA.MolecularWeight*1e-3
    M_B = gasB.MolecularWeight*1e-3
        
    # convert cm^5/molec^2 -> m^5/kg^2
    b_per_mass = b *(phys.N_avogadro**2)/(M_A*M_B) *(1e-2)**5

    # multiply by total density of gas [kg/m^3]
    kappa = b_per_mass * qA * qB * rho

    return kappa
