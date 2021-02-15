# PyRADS
PyRADS is the Python line-by-line RADiation model for planetary atmosphereS. PyRADS is a radiation code that can provide line-by-line spectral resolution, yet is written in Python and so is flexible enough to be useful in teaching.

For Earth-like atmospheres, PyRADS currently uses HITRAN 2016 line lists (http://hitran.org/) and the MTCKD continuum model (http://rtweb.aer.com/continuum_frame.html).

Looking for a version of PyRADS that can deal with shortwave radiation (scattering)?
https://github.com/ddbkoll/PyRADS-shortwave


References:

(1) Koll & Cronin, 2018, https://doi.org/10.1073/pnas.1809868115.

# Installation
1) Download to your own computer.

2) [optional] Install the required libraries using conda:
- cd $PyRADS
- conda env create -f environment.yml
- conda activate pyrads

  NOTE: don't forget to type ``conda activate pyrads'' each time before using PyRADS.

3) Manually compile the MTCKD model:
- cd $PyRADS/DATA/MT_CKD_continuum/cntnm.H2O_N2/build
- (on a Mac if you are using gfortran installed with conda) make -f make_cntnm osxGNUCONDAdbl
- (on a Mac otherwise) make -f make_cntnm osxGNUdbl

4) Run test scripts

To compute outgoing longwave radiation (OLR) in W/m2 for a given surface temperature:
- cd $PyRADS/Test01.olr
- python compute_olr_h2o.py

To compute OLRs for a set of surface temperatures and save the resulting output to txt:
- cd $PyRADS/Test02.runaway
- python compute_olr_h2o.01.100RH.py


NOTE: resolution in test scripts was chosen for relative speed, not accuracy. For research-grade output and model intercomparisons, vertical and spectral resolution need to be increased. For some reference values, see Methods in Koll & Cronin (2018).


# Requirements
Python 2 or 3 with numpy and scipy.

For the MTCKD continuum model: gmake and gfortran.

# Acknowledgements
PyRADS makes use of HITRAN 2016 line lists (http://hitran.org/), AER's MTCKD continuum model (http://rtweb.aer.com/continuum_frame.html), and the PyTran script published by Ray Pierrehumbert as part of the courseware for "Principles of Planetary Climates" (https://geosci.uchicago.edu/~rtp1/PrinciplesPlanetaryClimate/). Brian Rose (http://www.atmos.albany.edu/facstaff/brose/) has helped improve the code.
