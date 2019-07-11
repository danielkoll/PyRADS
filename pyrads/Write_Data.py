from __future__ import division, print_function, absolute_import
import numpy as np
from . import phys

### -----------------------------------
### Helpers
class Dummy:
    pass


### -----------------------------------
# save (T,p,q) structure and corresponding SW fluxes.
#    wave0,wave1,dwave - floats
#    p,T,q             - 1d arrays
#    zenith_list       - list
#    SW flux lists     - lists of 1d arrays:
#                          [ shape(pressure),shape(pressure),...] -> overall len(zenith_list) long

def save_profile( OUTDIR,wave0,wave1,dwave,p,T,q, \
                      zenith_list,SWdirect_list,SWminus_list,SWplus_list ):
    
    filename = "output" + "_nA%.1f" % wave0 + "_nB%.1f" % wave1 + "_dn%.2e" % dwave + ".txt"

    ## create new/blank file:
    f = open(OUTDIR+filename,'w')
    f.close()
    
    ## save data:
    f = open(OUTDIR+filename,'a')
    for zenith,SWdirect,SWminus,SWplus in zip(zenith_list,SWdirect_list,SWminus_list,SWplus_list):
        f.write("wavenr_min,wavenr_max,dwave [cm^-1] = %.4f,%.4f,%.4f" % (wave0,wave1,dwave) )
        f.write("\n")
        f.write("zenith angle = %.4f" % (zenith) )
        f.write("\n")
        f.write("pressure [bar],\tTemperature [K],\tSpecific Humidity[-],\tdir SW flux [W/m2],\tdn SW flux [W/m2],\tup SW flux [W/m2]")
        f.write("\n")
        for p_i,T_i,q_i,SWdirect_i,SWminus_i,SWplus_i in zip(p,T,q,SWdirect,SWminus,SWplus):
            f.write("%.6e,\t%.6e,\t%.6e,\t%.6e,\t%.6e,\t%.6e" % (p_i/1e5,T_i,q_i,SWdirect_i,SWminus_i,SWplus_i) )
            f.write("\n")
        f.write("\n")
        f.write("\n")

    f.close()
    
