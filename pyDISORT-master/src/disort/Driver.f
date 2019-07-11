c     -------------------------------------------------------------------
c     Python wrapper to the DISORT radiative transfer solver
c
c     Author: Sebastian Gimeno Garcia
c
c
c     License:
c
c     Do whatever you want with this piece of code. Enjoy it. If you
c     find it helpful, think about the authors of DISORT and drink to
c     their health, and why not, also to mine.
c
c     If you find any bug, please let me now.
c      
c     Ref:
c     
c     K. Stamnes, SC. Tsay, W. Wiscombe and K. Jayaweera, Numerically
c     stable algorithm for discrete-ordinate-method radiative
c     transfer in multiple scattering and emitting layered media,
c     Appl Opt 27 (1988) (12), pp. 2502–2509.
c     -------------------------------------------------------------------
      
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     Driver based on:
c      
c     RCS version control information:
c     $Header: DISOTEST.f,v 2.1 2000/04/03 21:21:55 laszlo Exp $
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      SUBROUTINE  RUN(
     I                 MAXCLY, DTAUC, SSALB, MAXMOM, TEMPER,       
     I                 IPHAS, GG,
     I                 WVNMLO, WVNMHI, USRTAU, MAXULV, UTAU, NSTR,
     I                 USRANG, MAXUMU, UMU, MAXPHI, PHI, IBCND, FBEAM,
     I                 UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     I                 TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
     O                 RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     O                 ALBMED, TRNMED
     &               )
      
     
c    Runs test problems for DISORT and checks answers. These
c    problems test almost all logical branches in DISORT.

c    It is HIGHLY recommended that you use the code below as a template
c    for creating your own CALLs to DISORT, rather than starting from
c    scratch.  This will prevent a lot of mistakes and ensure that every
c    input argument gets a value.  Note in particular how GETMOM is
c    sometimes called to fill an array section of PMOM (for one layer);
c    several people have done this incorrectly in attempting to write it
c    ab initio (passing array sections for arrays that do not start at
c    element 1 is tricky).

c    Note that the ratio to the 'correct answer' may occasionally be
c    significantly different from unity -- even so different that
c    the ratio just prints as ****** rather than a number.  However,
c    this mostly occurs for values of flux or intensity that are very
c    small compared to the forcing functions (that is, small compared
c    to internal thermal emission and/or radiation incident at the
c    boundaries).  The printed number 'SERIOUSLY NON-UNIT RATIOS'
c    attempts to count just the cases where there is a real disagreement
c    and not those where quantitites are down at their noise level
c    (defined as 10^(-6) times their maximum value).

c    Further documentation can be found in the file DISOTEST.doc.


c  Routines called :

c    DISORT:   The discrete ordinates radiative transfer program

c    BDREF:    Sets bidirectional reflectance of lower boundary

c    GETMOM:   Sets phase function Legendre coefficients

c    PRTFIN:   Prints fluxes and intensities and their ratios to
c              the correct values

c    CHEKDO:   Data block containing correct fluxes and intensities

c    RATIO :   Ratio of calculated to correct value with underflow
c              and overflow protection (kept in file DISORT.f)

c       INPUT: IPHAS   Phase function options
c                      1 : Isotropic
c                      2 : Rayleigh
c                      3 : Henyey-Greenstein with asymmetry factor GG
c                      4 : Haze L as specified by Garcia/Siewert
c                      5 : Cloud C.1 as specified by Garcia/Siewert

c              GG      Asymmetry factor for Henyey-Greenstein case

c              NMOM    Index of highest Legendre coefficient needed
c                        ( = number of streams 'NSTR'  chosen
c                         for the discrete ordinate method)
c+---------------------------------------------------------------------+
c
c                 ** DISORT I/O specifications **

CF2PY INTENT(HIDE) :: MAXCLY, MAXUMU, MAXPHI, MAXULV
CF2PY INTENT(IN)   :: DTAUC, SSALB, TEMPER, MAXMOM
CF2PY INTENT(IN)   :: IPHAS, GG,
CF2PY INTENT(IN)   :: WVNMLO, WVNMHI, USRTAU, UTAU, NSTR
CF2PY INTENT(IN)   :: USRANG, UMU, PHI, IBCND, FBEAM
CF2PY INTENT(IN)   :: UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP
CF2PY INTENT(IN)   :: TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT
CF2PY INTENT(OUT)  :: RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU
CF2PY INTENT(OUT)  :: ALBMED, TRNMED
      
      INTEGER  MAXCLY, MAXULV, MAXUMU, MAXPHI
      INTEGER  MAXMOM
c$$       PARAMETER ( MAXMOM = 299)
      CHARACTER  HEADER*127
      LOGICAL  LAMBER, PLANK, ONLYFL, PRNT(5), USRANG, USRTAU
      INTEGER  IBCND, NMOM, NLYR, NUMU, NSTR, NPHI, NTAU
      REAL     ACCUR, ALBEDO, BTEMP, DTAUC( MAXCLY ), FBEAM, FISOT,
     &         PHI( MAXPHI ), PMOM( 0:MAXMOM, MAXCLY ),
     &         PHI0, SSALB( MAXCLY ), TEMPER( 0:MAXCLY ), TEMIS, TTEMP,
     &         WVNMLO, WVNMHI, UMU( MAXUMU ), UMU0, UTAU( MAXULV )

      REAL     RFLDIR( MAXULV ), RFLDN( MAXULV ), FLUP( MAXULV ),
     &         DFDT( MAXULV ), UAVG( MAXULV ),
     &         UU( MAXUMU, MAXULV, MAXPHI ), ALBMED( MAXUMU ),
     &         TRNMED( MAXUMU )

      INTEGER  IPHAS( MAXCLY )
      REAL     GG( MAXCLY )

c+---------------------------------------------------------------------+

      INTEGER  MXTAU, MXMU, MXPHI
      PARAMETER     ( MXTAU = 5, MXMU = 32, MXPHI = 3 )
      INTEGER       LC
      REAL          PI

c+---------------------------------------------------------------------+
c     .. External Subroutines ..

      EXTERNAL  DISORT, ERRMSG, GETMOM, PRTFIN
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ASIN, FLOAT, INDEX
c     ..


      PI = 2.* ASIN( 1.0 )

      NLYR = MAXCLY
      NMOM = NSTR
      DO LC = 1, NLYR
         CALL  GETMOM( IPHAS( LC ), GG( LC ), NMOM, PMOM(0,LC) )
      END DO 
      NTAU      = MAXULV
      NPHI      = MAXPHI
      NUMU      = MAXUMU
      HEADER = 'Python wrapper to the DISORT radiative transfer solver'

      CALL  DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
     &                 WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
     &                 USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
     &                 UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
     &                 TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
     &                 HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
     &                 MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
     &                 ALBMED, TRNMED )

      
      END   
