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
### MAIN

output_angles = pyrads.Merge_Spectral_Output.merge_output(".",prefix='output')

# COMPUTE ANGLE-AVERAGED QUANTITIES
Nangles = len(output_angles)

# if Nangles != 6:
#     print( "*** ERROR: number of angles doesn't match the expected value!! ***" )

## ---
SWdir_angular = np.zeros( (len(output_angles[0].SWdir),Nangles) )
SWdn_angular = np.zeros( (len(output_angles[0].SWdir),Nangles) )
SWup_angular = np.zeros( (len(output_angles[0].SWdir),Nangles) )
LWdn_angular = np.zeros( (len(output_angles[0].SWdir),Nangles) )
LWup_angular = np.zeros( (len(output_angles[0].SWdir),Nangles) )
albedo_angular = np.zeros( Nangles )
    
zenith_angular = [angle.zenith for angle in output_angles]
cosz_angular = [angle.cosz for angle in output_angles]

for angle,i in zip(output_angles,range(Nangles)):
    SWdir_angular[:,i] = angle.SWdir
    SWdn_angular[:,i] = angle.SWdn
    SWup_angular[:,i] = angle.SWup

    # FOR NOW: skip LW, not implemented yet...
    # LWdn_angular[:,i] = angle.LWdn
    # LWup_angular[:,i] = angle.LWup
        
    albedo_angular[i] = angle.SWup[0]/(angle.SWdir[0]+angle.SWdn[0])

# average over angles, using correct weights:
weight = np.sin( np.array(zenith_angular) * np.pi/180. )

SWdir = np.average( SWdir_angular,axis=-1,weights=weight )
SWdn = np.average( SWdn_angular,axis=-1,weights=weight )
SWup = np.average( SWup_angular,axis=-1,weights=weight )
LWdn = np.average( LWdn_angular,axis=-1,weights=weight )
LWup = np.average( LWup_angular,axis=-1,weights=weight )
albedo = np.average( albedo_angular,axis=-1,weights=weight )


print( "zenith angles", zenith_angular )
print( "albedo at each angle", albedo_angular )
print( "average albedo", albedo )
