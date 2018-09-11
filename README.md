# PyRADS
PyRADS is a Python line-by-line RADiation model for planetary atmosphereS. PyRADS is a radiation code that can provide line-by-line spectral resolution, yet is written in Python and so is flexible enough to be useful in teaching.

For Earth-like atmospheres, PyRADS currently uses HITRAN 2016 line lists (http://hitran.org/) and the MTCKD continuum model (http://rtweb.aer.com/continuum_frame.html). 

Currently, PyRADS is only valid for longwave calculations (no scattering).

References:

(1) Koll & Cronin, 2018. In press.

# Installation
1) Download to your own computer.
2) Manually compile the MTCKD model:
- cd $PyRADS/DATA/MT_CKD_continuum/cntnm.H2O_N2/build
- (on a Mac) gmake -f make_cntnm osxGNUdbl
3) Run test
- cd $PyRADS/Test01.Compute_OLR_H2O
- python compute_olr_h2o.01.100RH.py

# Requirements
Python 2.x with numpy and scipy.

For the MTCKD continuum model: gmake and gfortran.
