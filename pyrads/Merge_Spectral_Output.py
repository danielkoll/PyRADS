from __future__ import division, print_function, absolute_import
import numpy as np
import glob


### -----------------------------------
### Helpers
class Dummy:
    pass

# use this class to sum up over all spectral intervals:
class Output_Angle:
    def __init__(self,zenith,p,T,q):
        self.zenith = zenith
        self.cosz = np.cos( zenith * np.pi/180. )
        
        self.p = p
        self.T = T
        self.q = q
        self.SWdir = np.zeros_like(p)
        self.SWdn = np.zeros_like(p)
        self.SWup = np.zeros_like(p)

    def add_to_fluxes(self,SWdir,SWdn,SWup):
        self.SWdir += SWdir
        self.SWdn += SWdn
        self.SWup += SWup


# ---
def read_output(file):
    f = open(file,'r')
    txt = f.read()
    
    # files are arranged by zenith angle chunks:
    chunks = txt.split('wavenr_min,wavenr_max,dwave [cm^-1] =')

    angle_list = []

    for chunk in chunks[1:]:
        # note: genfromtxt needs to be fed a list of strings..
        data = np.genfromtxt( chunk.split('\n'), delimiter=',\t', skip_header=3 )

        # get zenith angle:
        #   is in 2nd line of each chunk
        zenith = float( chunk.split('\n')[1].split('zenith angle =')[-1] )
        cosz = np.cos( zenith*np.pi/180. )

        # save into a out object:
        out = Dummy()
        out.zenith = zenith
        out.cosz = cosz
        
        out.p = data[:,0]
        out.T = data[:,1]
        out.q = data[:,2]
        out.SWdir = data[:,3]
        out.SWdn = data[:,4]
        out.SWup = data[:,5]

        angle_list.append(out)

    return angle_list


### -----------------------------------
### MAIN UTILITY FUNCTION:

def merge_output(path,prefix='output',remove_out=False,remove_err=False):
    # if path doesn't end in '/', make sure if does:
    if path[-1] is not '/':
        path = path + '/'
    
    files = glob.glob( path + prefix + '*.txt' )

    # get data from each spectral interval:
    spectral_list = []
    for ff in files:
        x = read_output(ff)
        spectral_list.append( x )

    # prepare output:
    #    (p,T,q) should be same for all spectral intervals, so just pick any one.
    spec0 = spectral_list[-1]  # pick any one
    N_angles = len(spec0)
    output_list = [Output_Angle(angle.zenith,angle.p,angle.T,angle.q) for angle in spec0]

    # to get net fluxes: add up all fluxes
    for i in range(N_angles):
        output = output_list[i]

        for spectral_interval in spectral_list:
            SWdir = spectral_interval[i].SWdir
            SWdn = spectral_interval[i].SWdn
            SWup = spectral_interval[i].SWup

            output.add_to_fluxes(SWdir,SWdn,SWup)

    # ** UNDER CONSTRUCTION ** 
    # (remove std out/err files?)
    if remove_out:
        pass
    if remove_err:
        pass

    return output_list
