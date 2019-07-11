from __future__ import division, print_function, absolute_import
import numpy as np
from . import phys
import os

####
# DEFINITIONS:


#Path to the datasets
datapath = '/'.join( os.path.abspath(__file__).split('/')[:-2] ) + '/DATA/VPL_solar_spectrum/'  # !

#
solar_file = datapath + 'vpl_sunum_extended.txt'


### 
# Return a normalized solar spectrum by interpolating the VPL spectrum.
#   returns units = W/m^2/cm-1
#
#   paper: Segura et al, 2003, Astrobiology
#   link: vpl.astro.washington.edu/spectra/stellar/other_stars.htm
#
#     CAREFUL: how to deal with out of bounds interpolation?
#              for now, set to zero!

def get_solar_fn(Lstar):

    # 1) read vpl file
    data = np.genfromtxt(solar_file,skip_header=1)
    wvl,sun = data[:,0],data[:,1]
    sun = sun*1e4   # convert W/cm2/micron to W/m^2/micron

    # 2) now convert everything to W/m^2/cm-1
    #    note: B_n = B_lambda * (-1)*lambda^2
    n = 1./(wvl*1e-4)
    n = n[::-1]   # flip so monotonically increasing...

    sun_n = sun*wvl*(wvl*1e-4)   # W/m2/micron *(micron) *(cm)
    sun_n = sun_n[::-1]

    sun_n = sun_n/np.trapz(sun_n,x=n) * Lstar     # normalize
    print( "(get_solar_fn) integrated solar source fn: %.4f W/m2" % (np.trapz(sun_n,x=n)) ) # TESTING


    # 3) define interpolation fn
    oob = 0.    # out of bounds: zero or nan?
    get_solar = lambda nu: np.interp(nu,n,sun_n,left=oob,right=oob)

    return get_solar
