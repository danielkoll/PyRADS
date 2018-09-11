C     path:      $Source$
C     author:    $Author: malvarad $
C     revision:  $Revision: 16459 $
C     created:   $Date: 2012-10-22 13:21:52 -0400 (Mon, 22 Oct 2012) $
C 
C
C  --------------------------------------------------------------------------
C |  Copyright ©, Atmospheric and Environmental Research, Inc., 2011         |
C |                                                                          |
C |  All rights reserved. This source code is part of the MT_CKD continuum   |
C |  software and is designed for scientific and research purposes.          |
C |  Atmospheric and Environmental Research, Inc. (AER) grants USER          |
C |  the right to download, install, use and copy this software              |
C |  for scientific and research purposes only. This software may be         |
C |  redistributed as long as this copyright notice is reproduced on any     |
C |  copy made and appropriate acknowledgment is given to AER. This          |
C |  software or any modified version of this software may not be            |
C |  incorporated into proprietary software or commercial software           |
C |  offered for sale.                                                       |
C |                                                                          |
C |  This software is provided as is without any express or implied          |
C |  warranties.                                                             |
C |                       (http://www.rtweb.aer.com/)                        |
C  --------------------------------------------------------------------------
C
                  PROGRAM DRCNTNM
C
C
C
c     The mt_ckd water vapor continuum is a completely new continuum  
c     formulation based on a collision induced  component and a sub-Lorentzian 
c     line wing component.  Both the water vapor continuum coefficients and 
c     those for other molecules are constrained to agree with accurate
c     measurements of continuum absorption in the spectral regions where such
c     measurements exist.
c
c     This is an updated version of the continuum program:
c     this version provides optical depths on file CNTNM.OPTDPT as before:
c     it also provides the continuum coefficients on file  WATER.COEF
c
c     the length of the header records may vary by version:
c         in this version the WATER.COEF header information is 47 records 
c         in this version the CNTNM.OPTDT header information is 34 records 
c
c     presumably the user will want to create an input file to address
c     individual requirements
C
C
      IMPLICIT REAL*8           (V)                                     ! F00030
c
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)              ? 500060
C
      COMMON /CVRCNT/ HNAMCNT,HVRCNT
c
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       F00130
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   F00140
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    F00150
c
      Common /share/ HOLN2
c
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2,GRAV,CPDAIR,AIRMWT,SECDY  
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           A03020
      COMMON /XCONT/  V1C,V2C,DVC,NPTC,C(6000) 
c
c********************************************
      COMMON /cnth2o/ V1h,V2h,DVh,NPTh,Ch(5050),csh2o(5050),cfh2o(5050)
c********************************************
c
      COMMON /IFIL/ IRD,IPRcnt,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,
     *              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,        F00180
     *              NLTEFL,LNFIL4,LNGTH4                                  F00190

      common /cntscl/ XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL

c------------------------------------
c for analytic derivative calculation
c note: ipts  = same dimension as ABSRB
c       ipts2 = same dimension as C
      parameter (ipts=5050,ipts2=6000)
      common /CDERIV/ icflg,iuf,v1absc,v2absc,dvabsc,nptabsc,
     &    dqh2oC(ipts),dTh2oC(ipts),dUh2o,
     &    dqco2C(ipts),dTco2C(ipts),
     &    dqo3C(ipts),dTo3C(ipts),
     &    dqo2C(ipts),dTo2C(ipts),
     &    dqn2C(ipts),dTn2C(ipts)

      real ctmp(ipts2),cself(ipts),cforeign(ipts),ch2o(ipts2)
      real ctmp2(ipts2),ctmp3(ipts2),ctmp4(ipts),ctmp5(ipts)

c------------------------------------
c
      dimension xcnt(7)
c
      equivalence (xself,xcnt(1))
c
      CHARACTER*18 HNAMCNT,HVRCNT
C                                                                         F00100
      CHARACTER*8      XID,       HMOLID,      YID 
      REAL*8               SECANT,       XALTZ
c
      character*8 holn2
C                                                                         F00120
      DATA XLOSMT/2.68675E+19/                                            50024
C
      RADCN1=2.*PLANCK*CLIGHT*CLIGHT*1.E-07                                 3070
      RADCN2=PLANCK*CLIGHT/BOLTZ                                            3080
c
      icflg = -999
C
      do 1, i=1,7
         xcnt(i)=1.
 1    continue

      do 2, i=1,5050
         absrb(i)=0.
 2    continue

      do 3, i=1,60
         wk(i)=0.
 3    continue
c
      ird = 55
      ipr = 66
      ipu = 7
c
      OPEN (ipr,FILE='CNTNM.OPTDPT')
      OPEN (ipu,FILE='WATER.COEF')
c
      print *
      print *, '  This version is limited to 5000 values  '
C
C   THIS PROGRAM CALCULATES THE CONTINUUM OPTICAL DEPTH
C         FOR AN HOMOGENEOUS LAYER
C
C   THE FOLLOWING QUANTITIES MUST BE SPECIFIED: 
C
C          PRESSURE                   PAVE (MB)
C
C          TEMPERATURE                TAVE ( K)
C
C          COLUMN AMOUNT
C            NITROGEN                 WN2    (MOLEC/CM**2)
C            OXYGEN                   WK(7)  (MOLEC/CM**2)
C            CARBON DIOXIDE           WK(2)  (MOLEC/CM**2)
C            WATER VAPOR              WK(1)  (MOLEC/CM**2)
C
C          NUMBER OF MOLECULES        NMOL
C
C          BEGINNING WAVENUMBER       V1ABS (CM-1)
C
C          ENDING WAVENUMBER          V2ABS (CM-1)
C
C          SAMPLING INTERVAL          DVABS (CM-1)
C
C          NUMBER OF VALUES           NPTABS
C
C
C   THE RESULTS ARE IN ARRAY ABSORB
C
C   NOTE THAT FOR AN ATMOSPHERIC LAYER: 
C
C            WTOT   = XLOSMT * (PAVE/1013) * (273/TAVE) * (PATH LENGTH)
C
C            WBROAD = the column amount for all species not explicitly provided
C
C            WK(M)  = (VOLUME MIXING RATIO) * (COLUMN OF DRY AIR)
C
C
      iprcnt = ipr
                   CALL PRCNTM 
      iprcnt = ipu
                   CALL PRCNTM
C
C   THE FOLLOWING IS AN EXAMPLE FOR A ONE CM PATH (SEE CNTNM.OPTDPT FOR RESULTS)
C

      PAVE = 1013.
      TAVE =  296.
c
      VMRH2O = 0.01
C
      xlength = 1.
c
      print *
      print *,' *** For this program, vmr_h2o is taken ',
     * 'with respect to the total column ***'

      print *
      print *,' read: pressure (mb)  if negative use default values'
      read *, press_rd
      
      if (press_rd .gt. 0.) then
         pave = press_rd
         print *,' read:   temperature (K)'
         read *, tave
         print *,' read:   path length (cm)'
         read *, xlength
         print *,' read:   vmr h2o '
         read *, vmrh2o
      endif
      print *,
     * 'Pressure (mb), Temperature (K), Path Length (cm),    VMR H2O'

      print 911, pave,tave,xlength,vmrh2o
 911  format(1x,f13.6,f17.4,f18.4,f12.8)
C
c     It may be preferable to specifiy the water column directly!
c
      WTOT = XLOSMT*(PAVE/1013.)*(273./TAVE)* xlength
C
      W_dry = WTOT * (1.-VMRH2O)
c
c     ww is column of dry air;  vol. mix. ratios are based on dry air
c
c argon:
      WA     = 0.009     * W_dry
c nitrogen:
      WN2    = 0.78      * W_dry
c oxygen:
      WK(7)  = 0.21      * W_dry

c carbon dioxide:
      WK(2)  = 345.E-06  * W_dry

c      WK(2) = 0.

c water vapor:
      if (abs(vmrh2o-1.) .lt. 1.e-05) then
         wk(1) = wtot
      else
         WK(1) = VMRH2O * W_dry
      endif
C
      WBROAD=WN2+WA
c
      NMOL = 7
c
      V1ABS =    0.
      V2ABS = 10000. 
 
      DVABS =    2.
c ..........................................................
c      write (*,*) '  v1abs,  v2abs,  dvabs  '
c      read  (*,*)    v1abs,  v2abs,  dvabs
c ..........................................................

      NPTABS =  1. + (v2abs-v1abs)/dvabs

      do 85 i=1,nptabs
         absrb(i) =0.
 85   continue

cc
      WRITE (IPR,970) PAVE,TAVE
      WRITE (IPR,975) (HMOLID(I),I=1,7),HOLN2                             A23490
      WRITE (IPR,980) (WK(M),M=1,7),WBROAD
c
      WRITE (IPu,970) PAVE,TAVE
      WRITE (IPu,975) (HMOLID(I),I=1,7),HOLN2                             A23490
      WRITE (IPu,980) (WK(M),M=1,7),WBROAD
C
  970 FORMAT (/,29x, 'P(MB)',7X,'T(K)', //,23x,0P,F12.3,F9.2)
  975 FORMAT (/, 9X,'MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER ',//,
     *         8(1X,A6,3X))
  980 FORMAT (/,1P,8E10.3,//)
C
      jrad=1
c
      v1 = v1abs
      v2 = v2abs
c
      CALL CONTNM(JRAD)
C
      DO 100 I=1,NPTABS
      VI=V1ABS+FLOAT(I-1)*DVABS
100   WRITE (ipr, 910) VI, ABSRB(I) 
910   FORMAT(F10.3,1P,E12.3)
C
      WRITE (7,920) tave
  920 FORMAT(//,' self and foreign water vapor continuum coefficients ',/,
     x       'for  ',f8.2,'K - ',   //,
     x       ' the self-continuum scales as ( Rself/Ro ) ',/,
     x       ' the foreign continuum scales as ( (Rtot-Rself)/Ro ) ',/,
     x       ' where R is the density rho [ R = (P/Po)*(To/T) ]. ',//,  
     x   10x,'     without radiation field:  ',
     x   10x,'       with radiation field:   ',      /,
     x   10x,'      self         foreign     ',
     x   10x,'      self         foreign     ',      /,
     x       '    cm-1  ',                    
     x       '       1/(cm-1 molec/cm**2)    ',
     x   10x,'       1/(molec/cm**2) '        ,//)
c
      xkt=tave/radcn2

      do 200 i=1,npth
      vi=v1h+float(i-1)*dvh
      if (vi.ge.v1abs .and. vi.le.v2abs) then
         radfld=radfn(vi,xkt)
         csh2or=csh2o(i) * radfld
         cfh2or=cfh2o(i) * radfld
         write (ipu,930) vi, csh2o(i), cfh2o(i), csh2or, cfh2or
  930    format(f10.2, 1p, 2e15.4,10x, 1p, 2e15.4)
      endif
  200 continue
c
      END 
c**********************************************************************
      Block Data phys_consts
c
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2,GRAV,CPDAIR,AIRMWT,SECDY 
c
      DATA PI /3.1415926535898 /
c
c    Constants from NIST 01/11/2002
c
      DATA PLANCK / 6.62606876E-27 /, BOLTZ  / 1.3806503E-16 /,
     *     CLIGHT / 2.99792458E+10 /, 
     *     AVOGAD / 6.02214199E+23 /, ALOSMT / 2.6867775E+19 /,
     *     GASCON / 8.314472  E+07 /
     *     RADCN1 / 1.191042722E-12 /, RADCN2 / 1.4387752    /
c
c     Pi was obtained from   PI = 2.*ASIN(1.)                             A03980
c
c     units are generally cgs
c
c     The first and second radiation constants are taken from NIST.
c     They were previously obtained from the relations:
c                            RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07      A03990
c                            RADCN2 = PLANCK*CLIGHT/BOLTZ                 A04000
      end
c
      BLOCK DATA                                                          A07600
C
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           A03020
C
      DATA ONEPL/1.001/, ONEMI/0.999/, ARGMIN/34./
C                                                                         A07710

      END                                                                 A07720
      BLOCK DATA cntnm
c
      IMPLICIT REAL*8           (V)                                     ! F00030
c
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),       F00130
     *                WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND,   F00140
     *                EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    F00150
c
      Common /share/ HOLN2
C                                                                         F00100
      CHARACTER*8      XID,       HMOLID,      YID 
      REAL*8               SECANT,       XALTZ
c
      character*8 holn2
c
      DATA HMOlid/ '  H2O   ' , '  CO2   ' , '   O3   ' , '  N2O   ' ,    A12260
     *             '   CO   ' , '  CH4   ' , '   O2   ' , 53*'        '/
c
      DATA HOLN2 / ' OTHER'/                                              A19810
c
      end
      SUBROUTINE XINT (V1A,V2A,DVA,A,AFACT,VFT,DVR3,R3,N1R3,N2R3)         B17520
C                                                                         B17530
C                                                                         B17870
      IMPLICIT REAL*8           (V)                                     ! F00030
C                                                                         B17550
C     THIS SUBROUTINE INTERPOLATES THE A ARRAY STORED                     B17560
C     FROM V1A TO V2A IN INCREMENTS OF DVA USING A MULTIPLICATIVE         B17570
C     FACTOR AFACT, INTO THE R3 ARRAY FROM LOCATION N1R3 TO N2R3 IN       B17580
C     INCREMENTS OF DVR3                                                  B17590
C                                                                         B17600
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           B17610
      DIMENSION A(*),R3(*)                                                B17620
C                                                                         B17630
      RECDVA = 1./DVA                                                     B17640
      ILO = (V1A+DVA-VFT)/DVR3+1.+ONEMI                                   B17650
      ILO = MAX(ILO,N1R3)                                                 B17660
      IHI = (V2A-DVA-VFT)/DVR3+ONEMI                                      B17670
      IHI = MIN(IHI,N2R3)                                                 B17680
C                                                                         B17690
      DO 10 I = ILO, IHI                                                  B17700
         VI = VFT+DVR3*FLOAT(I-1)                                         B17710
         J = (VI-V1A)*RECDVA+ONEPL                                        B17720
         VJ = V1A+DVA*FLOAT(J-1)                                          B17730
         P = RECDVA*(VI-VJ)                                               B17740
         C = (3.-2.*P)*P*P                                                B17750
         B = 0.5*P*(1.-P)                                                 B17760
         B1 = B*(1.-P)                                                    B17770
         B2 = B*P                                                         B17780
         CONTI = -A(J-1)*B1+A(J)*(1.-C+B2)+A(J+1)*(C+B1)-A(J+2)*B2        B17790
         R3(I) = R3(I)+CONTI*AFACT                                        B17800
   10 CONTINUE                                                            B17810
C                                                                         B17820
      RETURN                                                              B17830
C                                                                         B17840
      END                                                                 B17850
      FUNCTION RADFN (VI,XKT)                                             B17860
C                                                                         B17870
      IMPLICIT REAL*8           (V)                                     ! F00030
C                                                                         B17890
C     FUNCTION RADFN CALCULATES THE RADIATION TERM FOR THE LINE SHAPE     B17900
C                                                                         B17910
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   B17920
C                                                                         B17930
C               LAST MODIFICATION:    12 AUGUST 1991                      B17940
C                                                                         B17950
C                  IMPLEMENTATION:    R.D. WORSHAM                        B17960
C                                                                         B17970
C             ALGORITHM REVISIONS:    S.A. CLOUGH                         B17980
C                                     R.D. WORSHAM                        B17990
C                                     J.L. MONCET                         B18000
C                                                                         B18010
C                                                                         B18020
C                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         B18030
C                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          B18040
C                                                                         B18050
C----------------------------------------------------------------------   B18060
C                                                                         B18070
C               WORK SUPPORTED BY:    THE ARM PROGRAM                     B18080
C                                     OFFICE OF ENERGY RESEARCH           B18090
C                                     DEPARTMENT OF ENERGY                B18100
C                                                                         B18110
C                                                                         B18120
C      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             B18130
C                                                                         B18140
C                                             FASCOD3                     B18150
C                                                                         B18160
C                                                                         B18170
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   B18180
C                                                                         B18190
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           B18200
C                                                                         B18210
C      IN THE SMALL XVIOKT REGION 0.5 IS REQUIRED                         B18220
C                                                                         B18230
      XVI = VI                                                            B18240
C                                                                         B18250
      IF (XKT.GT.0.0) THEN                                                B18260
C                                                                         B18270
         XVIOKT = XVI/XKT                                                 B18280
C                                                                         B18290
         IF (XVIOKT.LE.0.01) THEN                                         B18300
            RADFN = 0.5*XVIOKT*XVI                                        B18310
C                                                                         B18320
         ELSEIF (XVIOKT.LE.10.0) THEN                                     B18330
            EXPVKT = EXP(-XVIOKT)                                         B18340
            RADFN = XVI*(1.-EXPVKT)/(1.+EXPVKT)                           B18350
C                                                                         B18360
         ELSE                                                             B18370
            RADFN = XVI                                                   B18380
         ENDIF                                                            B18390
C                                                                         B18400
      ELSE                                                                B18410
         RADFN = XVI                                                      B18420
      ENDIF                                                               B18430
C                                                                         B18440
      RETURN                                                              B18450
C                                                                         B18460
      END                                                                 B18470
c*******
c*******
c*******
      Include 'contnm.f90'
