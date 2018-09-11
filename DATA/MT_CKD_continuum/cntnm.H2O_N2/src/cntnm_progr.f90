!     path:      $Source$
!     author:    $Author: jdelamer $
!     revision:  $Revision: 11007 $
!     created:   $Date: 2011-03-29 16:16:57 -0400 (Tue, 29 Mar 2011) $
! 
!
!  --------------------------------------------------------------------------
! |  Copyright ©, Atmospheric and Environmental Research, Inc., 2011         |
! |                                                                          |
! |  All rights reserved. This source code is part of the MT_CKD continuum   |
! |  software and is designed for scientific and research purposes.          |
! |  Atmospheric and Environmental Research, Inc. (AER) grants USER          |
! |  the right to download, install, use and copy this software              |
! |  for scientific and research purposes only. This software may be         |
! |  redistributed as long as this copyright notice is reproduced on any     |
! |  copy made and appropriate acknowledgment is given to AER. This          |
! |  software or any modified version of this software may not be            |
! |  incorporated into proprietary software or commercial software           |
! |  offered for sale.                                                       |
! |                                                                          |
! |  This software is provided as is without any express or implied          |
! |  warranties.                                                             |
! |                       (http://www.rtweb.aer.com/)                        |
!  --------------------------------------------------------------------------
!
                  PROGRAM DRCNTNM
!
!
!
!     The mt_ckd water vapor continuum is a completely new continuum  
!     formulation based on a collision induced  component and a sub-Lorentzian 
!     line wing component.  Both the water vapor continuum coefficients and 
!     those for other molecules are constrained to agree with accurate
!     measurements of continuum absorption in the spectral regions where such
!     measurements exist.
!
!     This is an updated version of the continuum program:
!     this version provides optical depths on file CNTNM.OPTDPT as before:
!     it also provides the continuum coefficients on file  WATER.COEF
!
!     the length of the header records may vary by version:
!         in this version the WATER.COEF header information is 47 records 
!         in this version the CNTNM.OPTDT header information is 34 records 
!
!     presumably the user will want to create an input file to address
!     individual requirements
!
!
      Use lblparams
      USE phys_consts
      USE planet_consts
      IMPLICIT REAL*8           (V)                                     ! F00030
!
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)              ! 500060
!
      COMMON /CVRCNT/ HNAMCNT,HVRCNT
!
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4),  &    ! F00130
                      WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &  ! F00140
                      EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF   ! F00150
!
      Common /share/ HOLN2
!
      COMMON /CONSTS/ CPDAIR,SECDY  
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                          ! A03020
      COMMON /XCONT/  V1C,V2C,DVC,NPTC,C(6000) 
!
!********************************************
      COMMON /cnth2o/ V1h,V2h,DVh,NPTh,Ch(5050),csh2o(5050),cfh2o(5050)
!********************************************
!
      COMMON /IFIL/ IRD,IPRcnt,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL, &
                   NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,   &    ! F00180
                   NLTEFL,LNFIL4,LNGTH4                                 ! F00190

      common /cntscl/ XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL

!------------------------------------
! for analytic derivative calculation
! note: ipts  = same dimension as ABSRB
!       ipts2 = same dimension as C
!      parameter (ipts=5050,ipts2=6000)
      common /CDERIV/ icflg,iuf,v1absc,v2absc,dvabsc,nptabsc, &
          dqh2oC(ipts),dTh2oC(ipts),dUh2o, &
          dqco2C(ipts),dTco2C(ipts), &
          dqo3C(ipts),dTo3C(ipts), &
          dqo2C(ipts),dTo2C(ipts), &
          dqn2C(ipts),dTn2C(ipts)

      real ctmp(ipts2),cself(ipts),cforeign(ipts),ch2o(ipts2)
      real ctmp2(ipts2),ctmp3(ipts2),ctmp4(ipts),ctmp5(ipts)

!------------------------------------
!
      dimension xcnt(7)
!
      equivalence (xself,xcnt(1))
!
      CHARACTER*18 HNAMCNT,HVRCNT
!                                                                         F00100
      CHARACTER*8      XID,       HMOLID,      YID 
      REAL*8               SECANT,       XALTZ
!
      character*8 holn2
!                                                                         F00120
!      DATA XLOSMT/2.68675E+19/                                            !50024
!
!      RADCN1=2.*PLANCK*CLIGHT*CLIGHT*1.E-07                                ! 3070
!      RADCN2=PLANCK*CLIGHT/BOLTZ                                           ! 3080
!
      
      icflg = -999
!
! DKOLL: do I want to use all continua?
!        'xcnt' are flags for different calculations, see "ICNTNM Value" in contnm.f90
!       do 1, i=1,7
!          xcnt(i)=1.
!  1    continue

! DKOLL:
      xcnt(1)=1.   ! h2o self broadening?
      xcnt(2)=1.   ! h2o foreign broadening?
      xcnt(3)=0.   ! co2
      xcnt(4)=0.   ! o3
      xcnt(5)=0.   ! o2
      xcnt(6)=0.   ! n2
      xcnt(7)=0.   ! rayleigh scattering

      do 2, i=1,5050
         absrb(i)=0.
 2    continue

      do 3, i=1,60
         wk(i)=0.
 3    continue
!
      ird = 55
      ipr = 66
      ipu = 7
!
      OPEN (ipr,FILE='CNTNM.OPTDPT')
      OPEN (ipu,FILE='WATER.COEF')
!
      print *
      print *, '  This version is limited to 5000 values  '
!
!   THIS PROGRAM CALCULATES THE CONTINUUM OPTICAL DEPTH
!         FOR AN HOMOGENEOUS LAYER
!
!   THE FOLLOWING QUANTITIES MUST BE SPECIFIED: 
!
!          PRESSURE                   PAVE (MB)
!
!          TEMPERATURE                TAVE ( K)
!
!          COLUMN AMOUNT
!            NITROGEN                 WN2    (MOLEC/CM**2)
!            OXYGEN                   WK(7)  (MOLEC/CM**2)
!            CARBON DIOXIDE           WK(2)  (MOLEC/CM**2)
!            WATER VAPOR              WK(1)  (MOLEC/CM**2)
!
!          NUMBER OF MOLECULES        NMOL
!
!          BEGINNING WAVENUMBER       V1ABS (CM-1)
!
!          ENDING WAVENUMBER          V2ABS (CM-1)
!
!          SAMPLING INTERVAL          DVABS (CM-1)
!
!          NUMBER OF VALUES           NPTABS
!
!
!   THE RESULTS ARE IN ARRAY ABSORB
!
!   NOTE THAT FOR AN ATMOSPHERI! LAYER: 
!
!            WTOT   = ALOSMT * (PAVE/1013) * (273/TAVE) * (PATH LENGTH)
!
!            WBROAD = the column amount for all species not explicitly provided
!
!            WK(M)  = (VOLUME MIXING RATIO) * (COLUMN OF DRY AIR)
!
!
      iprcnt = ipr
                   CALL PRCNTM 
      iprcnt = ipu
                   CALL PRCNTM
!
!   THE FOLLOWING IS AN EXAMPLE FOR A ONE CM PATH (SEE CNTNM.OPTDPT FOR RESULTS)
!

      PAVE = 1013.
      TAVE =  296.
!
      VMRH2O = 0.01
!
      xlength = 1.
!
      print *
      print *,' *** F      USE phys_constsor this program, vmr_h2o is taken ', &
       'with respect to the total column ***'

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
      print *, &
       'Pressure (mb), Temperature (K), Path Length (cm),    VMR H2O'

      print 911, pave,tave,xlength,vmrh2o
 911  format(1x,f13.6,f17.4,f18.4,f12.8)
!
!     It may be preferable to specifiy the water column directly!
!
      WTOT = ALOSMT*(PAVE/1013.)*(273./TAVE)* xlength
!
      W_dry = WTOT * (1.-VMRH2O)
!
!     ww is column of dry air;  vol. mix. ratios are based on dry air
!
! argon:
      WA     = 0.0    ! DKOLL
! nitrogen:
      WN2    = W_dry  ! DKOLL
! oxygen:
      WK(7)  = 0.0    ! DKOLL

! carbon dioxide:
      WK(2)  = 0.0    ! DKOLL

! ozone:
      WK(3) = 0.

! water vapor:
!       if (abs(vmrh2o-1.) .lt. 1.e-05) then
!          wk(1) = wtot     ! DKOLL: in case atm is (almost) pure H2O... ?
!       else
!          WK(1) = VMRH2O * W_dry
!       endif

      WK(1) = WTOT * VMRH2O  ! DKOLL: ??

!
      WBROAD=WN2+WA
!
      NMOL = 7
!
      V1ABS =     0.
      V2ABS = 10000. 
 
      DVABS =    2.
! ..........................................................
!      write (*,*) '  v1abs,  v2abs,  dvabs  '
!      read  (*,*)    v1abs,  v2abs,  dvabs
! ..........................................................

      NPTABS =  1. + (v2abs-v1abs)/dvabs

      do 85 i=1,nptabs
         absrb(i) =0.
 85   continue

!!
      WRITE (IPR,970) PAVE,TAVE
      WRITE (IPR,975) (HMOLID(I),I=1,7),HOLN2                             !A23490
      WRITE (IPR,980) (WK(M),M=1,7),WBROAD
!
      WRITE (IPu,970) PAVE,TAVE
      WRITE (IPu,975) (HMOLID(I),I=1,7),HOLN2                             !A23490
      WRITE (IPu,980) (WK(M),M=1,7),WBROAD
!
  970 FORMAT (/,29x, 'P(MB)',7X,'T(K)', //,23x,0P,F12.3,F9.2)
  975 FORMAT (/, 9X,'MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER ',//, &
              8(1X,A6,3X))
  980 FORMAT (/,1P,8E10.3,//)
!
      jrad=1
!
      v1 = v1abs
      v2 = v2abs
!
      CALL CONTNM(JRAD)
!
      DO 100 I=1,NPTABS
      VI=V1ABS+FLOAT(I-1)*DVABS
100   WRITE (ipr, 910) VI, ABSRB(I) 
910   FORMAT(F10.3,1P,E12.3)
!
      WRITE (7,920) tave

  920 FORMAT(//,' self and foreign water vapor continuum coefficients ',/, &
             'for  ',f8.2,'K - ',   //, &
             ' the self-continuum scales as ( Rself/Ro ) ',/, &
             ' the foreign continuum scales as ( (Rtot-Rself)/Ro ) ',/, &
             ' where R is the density rho [ R = (P/Po)*(To/T) ]. ',//,  &
         10x,'     without radiation field:  ',  &
         10x,'       with radiation field:   ',      /, &
         10x,'      self         foreign     ', &
         10x,'      self         foreign     ',      /, &
             '    cm-1  ',                     &
             '       1/(cm-1 molec/cm**2)    ', &
         10x,'       1/(molec/cm**2) '        ,//) 
!
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
!
      END 

!
      BLOCK DATA                                                          !A07600
!
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           !A03020
!
      DATA ONEPL/1.001/, ONEMI/0.999/, ARGMIN/34./
!                                                                         A07710

      END                                                                 !A07720
      BLOCK DATA cntnm
!
      IMPLICIT REAL*8           (V)                                     ! F00030
!
      COMMON /FILHDR/ XID(10),SECANT,PAVE,TAVE,HMOLID(60),XALTZ(4), &      !F00130
                      WK(60),PZL,PZU,TZL,TZU,WBROAD,DV ,V1 ,V2 ,TBOUND, &   !F00140
                      EMISIV,FSCDID(17),NMOL,LAYER ,YI1,YID(10),LSTWDF    !F00150
!
      Common /share/ HOLN2
!                                                                         F00100
      CHARACTER*8      XID,       HMOLID,      YID 
      REAL*8               SECANT,       XALTZ
!
      character*8 holn2
!
      DATA HMOlid/ '  H2O   ' , '  CO2   ' , '   O3   ' , '  N2O   ' ,  &  !A12260
                   '   CO   ' , '  CH4   ' , '   O2   ' , 53*'        '/
!
      DATA HOLN2 / ' OTHER'/                                              !A19810
!
      end
      SUBROUTINE XINT (V1A,V2A,DVA,A,AFACT,VFT,DVR3,R3,N1R3,N2R3)         !B17520
!                                                                         B17530
!                                                                         B17870
      IMPLICIT REAL*8           (V)                                     ! F00030
!                                                                         B17550
!     THIS SUBROUTINE INTERPOLATES THE A ARRAY STORED                     B17560
!     FROM V1A TO V2A IN INCREMENTS OF DVA USING A MULTIPLICATIVE         B17570
!     FACTOR AFACT, INTO THE R3 ARRAY FROM LOCATION N1R3 TO N2R3 IN       B17580
!     INCREMENTS OF DVR3                                                  B17590
!                                                                         B17600
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           !B17610
      DIMENSION A(*),R3(*)                                                !B17620
!                                                                         B17630
      RECDVA = 1./DVA                                                     !B17640
      ILO = (V1A+DVA-VFT)/DVR3+1.+ONEMI                                   !B17650
      ILO = MAX(ILO,N1R3)                                                 !B17660
      IHI = (V2A-DVA-VFT)/DVR3+ONEMI                                      !B17670
      IHI = MIN(IHI,N2R3)                                                 !B17680
!                                                                         B17690
      DO 10 I = ILO, IHI                                                  !B17700
         VI = VFT+DVR3*FLOAT(I-1)                                         !B17710
         J = (VI-V1A)*RECDVA+ONEPL                                        !B17720
         VJ = V1A+DVA*FLOAT(J-1)                                          !B17730
         P = RECDVA*(VI-VJ)                                               !B17740
         C = (3.-2.*P)*P*P                                                !B17750
         B = 0.5*P*(1.-P)                                                 !B17760
         B1 = B*(1.-P)                                                    !B17770
         B2 = B*P                                                         !B17780
         CONTI = -A(J-1)*B1+A(J)*(1.-C+B2)+A(J+1)*(C+B1)-A(J+2)*B2        !B17790
         R3(I) = R3(I)+CONTI*AFACT                                        !B17800
   10 CONTINUE                                                            !B17810
!                                                                         B17820
      RETURN                                                              !B17830
!                                                                         B17840
      END                                                                 !B17850
      FUNCTION RADFN (VI,XKT)                                             !B17860
!                                                                         B17870
      IMPLICIT REAL*8           (V)                                       !F00030
!                                                                         B17890
!     FUNCTION RADFN CALCULATES THE RADIATION TERM FOR THE LINE SHAPE     B17900
!                                                                         B17910
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    B17920
!                                                                         B17930
!               LAST MODIFICATION:    12 AUGUST 1991                      B17940
!                                                                         B17950
!                  IMPLEMENTATION:    R.D. WORSHAM                        B17960
!                                                                         B17970
!             ALGORITHM REVISIONS:    S.A. CLOUGH                         B17980
!                                     R.D. WORSHAM                        B17990
!                                     J.L. MONCET                         B18000
!                                                                         B18010
!                                                                         B18020
!                     ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         B18030
!                     840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          B18040
!                                                                         B18050
!----------------------------------------------------------------------   B18060
!                                                                         B18070
!               WORK SUPPORTED BY:    THE ARM PROGRAM                     B18080
!                                     OFFICE OF ENERGY RESEARCH           B18090
!                                     DEPARTMENT OF ENERGY                B18100
!                                                                         B18110
!                                                                         B18120
!      SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             B18130
!                                                                         B18140
!                                             FASCOD3                     B18150
!                                                                         B18160
!                                                                         B18170
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   B18180
!                                                                         !B18190
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           !B18200
!                                                                         !B18210
!      IN THE SMALL XVIOKT REGION 0.5 IS REQUIRED                         !B18220
!                                                                         !B18230
      XVI = VI                                                            !B18240
!                                                                         !B18250
      IF (XKT.GT.0.0) THEN                                                !B18260
!                                                                         !B18270
         XVIOKT = XVI/XKT                                                 !B18280
!                                                                         !B18290
         IF (XVIOKT.LE.0.01) THEN                                         !B18300
            RADFN = 0.5*XVIOKT*XVI                                        !B18310
!                                                                         !B18320
         ELSEIF (XVIOKT.LE.10.0) THEN                                     !B18330
            EXPVKT = EXP(-XVIOKT)                                         !B18340
            RADFN = XVI*(1.-EXPVKT)/(1.+EXPVKT)                           !B18350
!                                                                         !B18360
         ELSE                                                             !B18370
            RADFN = XVI                                                   !B18380
         ENDIF                                                            !B18390
!                                                                         !B18400
      ELSE                                                                !B18410
         RADFN = XVI                                                      !B18420
      ENDIF                                                               !B18430
!                                                                         !B18440
      RETURN                                                              !B18450
!                                                                         !B18460
      END                                                                 !B18470
!*******
!*******
!*******
!      Include 'contnm.f90'
