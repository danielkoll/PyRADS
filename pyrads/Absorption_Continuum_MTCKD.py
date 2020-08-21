from __future__ import division, print_function, absolute_import
import os
import numpy
from . import phys
import time

'''
A number of utilities to run the stand-alone MT_CDK continuum using python.
dkoll (07.11.2017)

UPDATE: (27.01.2018)
Fixed the scaling in 'get_H2OContinuum'.
Tau computed with these coefficients now matches tau natively computed by MTCKD.

UPDATE: (29.01.2018)
To allow these scripts to run parallel within the same dir,
create all input/output files within a dummy directory.
For a unique directory label, use time.time() stamp.
'''


### -----------------------------------
### Global definitions here

datapath = '/'.join( os.path.abspath(__file__).split('/')[:-2] ) + '/DATA/MT_CKD_continuum/'

mtckd_exe_earthlike = datapath + "cntnm/cntnm_v3.2_OS_X_gnu_dbl"
mtckd_exe_N2N2 = datapath + "cntnm.N2N2only/cntnm_v3.2_OS_X_gnu_dbl"
mtckd_exe_H2O_N2 = datapath + "cntnm.H2O_N2/cntnm_v3.2_OS_X_gnu_dbl"
mtckd_exe_H2O_N2_foreign = datapath + "cntnm.H2O_N2.foreign_only/cntnm_v3.2_OS_X_gnu_dbl"
mtckd_exe_H2O_N2_self = datapath + "cntnm.H2O_N2.self_only/cntnm_v3.2_OS_X_gnu_dbl"

mtckd_exe_default = mtckd_exe_H2O_N2      # !


### -----------------------------------
### Given wavenumber, (T,p,q), return total specific absorption crosssection for
### continuum. Include both H2O self-continuum plus N2-H2O continuum.
###
### INPUT: wave -> VECTOR;      p,T,q -> SCALARS!
###
### OUTPUT: kappa, units = m2/kg
###         Note: Assume input is in MKS!
###
# From MTCKD:
# "the self-continuum scales as ( Rself/Ro )
#  the foreign continuum scales as ( (Rtot-Rself)/Ro )
#  where R is the density rho [ R = (P/Po)*(To/T) ]."
#
# Also see Shine et al (2016), Eqn.1 (?)
#
### CAREFUL:
### even though the above suggest that I should rescale using rho_i, I'm rescaling using p_i.
### this is the only way I can get tau from this kappa to match up with the tau computed by MTCKD.
### --> check to see how source code is doing this!

def get_H2OContinuum(wave,T,p_tot,p_H2O,exe_file=None):
    # get crosssections MTCKD's 'pre-defined' (I'd have to change MTCKD code) spectral grid
    # Note: convert pressure in [Pa] to [mbar=hPa]
    # NOTE: the coefs are actually independent of pressure,pathlength,vmr!
    wave0,cross_self,cross_foreign = get_coefficients(p_tot/1e2,T,pathlength=0.,vmr_h2o=0.,exe_file=exe_file)

    kappa_self0 = cross_self * 0.1 * phys.N_avogadro/phys.H2O.MolecularWeight
    kappa_foreign0 = cross_foreign * 0.1 * phys.N_avogadro/phys.H2O.MolecularWeight  # here: also m_H2O

    # Combine coefficients. Apply rescaling here.
    kappa0 = kappa_self0 * (p_H2O/p_tot) + kappa_foreign0 * ((p_tot-p_H2O)/p_tot)  # self + foreign broadening

    # Now rescale for p,T
    #    (test passed: without this scaling, can't match MTCKD's optical depths at T>>296K. with, it can.)
    kappa0 = kappa0 * (p_tot/1.013e5)*(296./T)

    # now interpolate to whatever wavenumbers I want:
    kappa = numpy.interp( wave, wave0,kappa0, left=numpy.nan,right=numpy.nan )

    return kappa


### -----------------------------------
### Given wavenumber, (T,p,q), return total optical depth for
### continuum. Include both H2O self-continuum plus N2-H2O continuum.
###
### INPUT: wave -> VECTOR;      p,T,q -> SCALARS!
###
### OUTPUT: tau.
###      Assume input is in MKS!

def get_tau_H2OContinuum(wave,T,p_tot,p_h2o,q,params):

    R_mean = (1.-q)*params.R + q*params.Rv
    rho = p_tot / (T * R_mean)
    pathlength = p_tot / (rho * params.g)   # assume hydrostatic ...
    vmr_h2o = p_h2o/p_tot

    ### QUESTION: is path the total path, or scaled for particular constituent only?
    wave0,tau0 = get_tau(p_tot/1e2,T,pathlength*1e2,vmr_h2o)  # note: need [cm],[mbar]

    # now interpolate to whatever wavenumbers I want:
    tau = numpy.interp( wave, wave0,tau0, left=numpy.nan,right=numpy.nan )

    return tau





### -----------------------------------
### Given temperature & pressure, return absorption coefficients
###
### INPUT: p,T
### OUTPUT: kappa_self, kappa_foreign
###
def get_coefficients(p,T,pathlength,vmr_h2o,exe_file=None,clean_files=True):

    if exe_file==None:
        exe_file = mtckd_exe_default

    #
    stamp = time.time()
    tmp_dir = "./tmp_mtckd_%.5f" % stamp   # create a unique tmp dir
    os.mkdir( tmp_dir )
    os.chdir( tmp_dir )
    #print "current working directory: ",os.getcwd()

    run_mtckd(p,T,pathlength,vmr_h2o,exe_file)
    wave,kappa_self,kappa_foreign = read_mtckd_output("WATER.COEF")

    os.chdir('..')  # back to initial dir
    #print "current working directory: ",os.getcwd()

    #
    if clean_files:
        os.remove( tmp_dir + '/INPUT')       # clean up tmp files
        os.remove( tmp_dir + '/WATER.COEF')
        os.remove( tmp_dir + '/CNTNM.OPTDPT')
        os.rmdir( tmp_dir )

    return wave,kappa_self,kappa_foreign


### -----------------------------------
### Given temperature & pressure, return optical depth
###
### INPUT: p,T
### OUTPUT: tau
###
def get_tau(p,T,pathlength,vmr_h2o,exe_file=None,clean_files=True):

    if exe_file==None:
        exe_file = mtckd_exe_default

    #
    stamp = time.time()
    tmp_dir = "./tmp_mtckd_%.5f" % stamp   # create a unique tmp dir
    os.mkdir( tmp_dir )
    os.chdir( tmp_dir )

    run_mtckd(p,T,pathlength,vmr_h2o,exe_file)
    wave,tau = read_mtckd_output("CNTNM.OPTDPT")

    os.chdir('..')  # back to initial dir
    #

    if clean_files:
        os.remove( tmp_dir + '/INPUT')       # clean up tmp files
        os.remove( tmp_dir + '/WATER.COEF')
        os.remove( tmp_dir + '/CNTNM.OPTDPT')
        os.rmdir( tmp_dir )

    return wave,tau





### -----------------------------------
### Run MTCKD for given input parameters
###
### Note: mt_ckd needs [p] = mb, [path] = cm.
###

def run_mtckd(p,T,pathlength,vmr_h2o,exe_file,suppress_stdout=True):

    input_file = "./INPUT"
    f = open(input_file,'w')
    f.write( "%.5f\n" % p )
    f.write( "%.5f\n" % T )
    f.write( "%.5f\n" % pathlength )
    f.write( "%.5f\n" % vmr_h2o )
    f.close()

    ##
    if suppress_stdout:
        #os.system( exe_file + " < " + input_file + " > /dev/null")
        os.system( '"'+exe_file+'"' + " < " + input_file + " > /dev/null")  # add double quotes in case space is in path
    else:
        print( "================================================")
        print( "===   RUN MT_CKD CONTINUUM                   ===")
        print( "===")
        os.system( exe_file + " < " + input_file )
        print( "===")
        print( "===")
        print( "================================================")


### -----------------------------------
### Given output txt file, cut off header and extract data
###
### INPUT:
###  need to specify which output file; 'WATER.COEF' vs 'CNTNM.OPTDPT'
###
### OUTPUT:
###  if 'WATER.COEF', return wavenumber, self and foreign coefficients [units = cm2/molec]
###                   so, the coefficient convolved with radiation field
###  if 'CNTNM.OPTDPT', return wavenumber and optical depth

def read_mtckd_output(file):

    if file=="WATER.COEF":
        f = open(file,'r')
        txt = f.read()
        body = txt.split("MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER \n")[-1]

        # note: genfromtxt needs to be fed a list of strings..
        data = numpy.genfromtxt(body.split('\n'), skip_header=20)  # careful about cutoffs!

        wave = data[:,0]
        kappa_self = data[:,3]
        kappa_foreign = data[:,4]
        return wave,kappa_self,kappa_foreign

    elif file=="CNTNM.OPTDPT":
        f = open(file,'r')
        txt = f.read()
        body = txt.split("MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER \n")[-1]

        data = numpy.genfromtxt(body.split('\n'), skip_header=6)  # careful about cutoffs!

        wave = data[:,0]
        tau = data[:,1]
        return wave,tau
