*  
* $Id: my_tmshel.F,v 1.2 2008-06-26 09:35:56 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/07/11 15:15:48  aplin
* Initial version
*
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
c      SUBROUTINE TMSHEL(NAME,NPT,LIST,IOPT,VV0,EE0,CHI2)
      SUBROUTINE MY_TMSHEL(NPT,LIST,IOPT,VV0,EE0,CHI2)
C****************************************************************************
c* Version number: 0.0      Date of last modification  26.08.99             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
C* Main modifiers from ALEPH to BRAHMS code:                                *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
C-----**************************************************---------C
C! FAST FITTING ROUTINE FOR SCATTERED HELICES IN THE TPC         C
C!    AUTHORS:  N. CHERNOV, G. OSOSKOV & M. POPPE                C
C!    Modified: M.Cattaneo 97/12/05
C!              Protect against divide by zero when DENOM=0.0
C!    Modified: M.Cattaneo 98/09/28
C!              Protect against divide by zero when DENOM=0.0, again...
C!    INPUT: NAME  = NAME OF THE BANK CONTAINING COORDINATES     C
C!            NPT  = NUMBER OF POINTS TO BE FITTED               C
C!            LIST = LIST OF POINTS TO BE FITTED                 C
C!            IOPT = 0 -> OBSOLETE                               C
C!                   1 -> WEIGHTED WITH 1/SIMA(R*PHI)**2         C
C!                   2 -> ERROR MATRIX CALCULATED                C
C!                   3 -> 3-DIMENSIONAL ITERATION                C
C!  OUTPUT:   VV0 = 1/R*CHARGE   [1/CM]  POS. IF CLOCKWISE       C
C!                  TAN(LAMBDA)  {=DZ/DS}TAN(ANGLE TO X,Y PLANE) C
C!                  PHI0         {0,2PI} ANGLE TO X-AXIS         C
C!                  D0*SIGN      [CM]    MINIMAL DIST. TO Z-AXIS,C
C!                                       POS. IF AXIS ENCIRCLED  C
C!                  Z0           [CM]    Z POS AT R=D0           C
C!            EE0 = ERROR MATRIX IN TRIANG. FORM                 C
C!            CHI2= CHI SQUARED = SUM (DEVIATIONS/ERRORS)**2     C
C!  NOTE: DEGREES OF FREEDOM = 2*NPT-5                           C
C!**WARNING***WARNING***WARNING***WARNING***WARNING***WARNING****C
C!THIS SUBROUTINE IS NOT PROTECTED AGAINNST MISSING INPUT BANKS  C
C!**WARNING***WARNING***WARNING***WARNING***WARNING***WARNING****C
C                                                                C
C----------------------------------------------------------------C
C     BASED ON  SUBROUTINE CIRCLE                                C
C     REFERENCE:  COMPUTER PHYSICS COMMUNICATIONS VOL 33,P329    C
C?
C!======================================================================
C----------------------------------------------------------------------
C my stuff:
#include "include/padrow.inc"
#include "include/fndparms.inc"
#include "include/tpco.inc"
c#include "gcflag.inc"
C GAB:
C(KH) +SEQ,VXPGEO.
C(KH) +SEQ,VXDCHAIN.
C----------------------------------------------------------------------
C GAB:
C      PARAMETER (EPS = 1.0E-16, ITMAX =15, MPT=40)
      PARAMETER (EPS = 1.0E-16, ITMAX =15, MPT=LTPDRO)
      REAL   PF(MPT),RF(MPT),PC(MPT),SP2(MPT),VV0(*),EE0(*),VV1(5),
     1       DEL(MPT),SXY(MPT),ZF(MPT),WZF(MPT),SS0(MPT),EEE(MPT),
     2       DELZ(MPT),GRAD(5),COV(15),DV(5),GGRAD(5),DRPHI(*),
     3       DZ(*),DCHI2(*),DELN(MPT),DELZN(MPT)
      DOUBLE PRECISION XF(MPT),YF(MPT),WF(MPT),XMID,YMID
      DOUBLE PRECISION ALF,ALM,A0,A1,A2,A22,BEM,BET,CUR,
     1   DD,DEN,DET,DY,D2,F,FACT,FG,F1,G,GAM,GAM0,GMM,G1,
     2   H,H2,P2,Q2,RM,RN,ROOT,
     3   XA,XB,XD,XI,XM,XX,XY,X1,X2,DEN2,
     4   YA,YB,YD,YI,YM,YY,   Y1,Y2,WN,SA2B2,DD0,CU2,PHIC
      INTEGER LIST(*)
c gab:
c      CHARACTER*4 NAME
      LOGICAL FIRST

      SAVE

      DATA FIRST/.TRUE./
      DATA SCACO /0.0160/
C     MEAN ANGLE OF DEFLECTION FOR ARGON, B=1.5 TESLA, AND
C     ALL RADII IN CM. ( <- MULTIPLE SCATTERING )
C
C     OFFSET FOR COORDINATE "IJ" IN 'TPCO'
C GAB:
C      KKTPCO(IJ) = KTPCO+2+(IJ-1)*(IW(KTPCO+1))
C     OFFSET FOR COORDINATE "IJ" IN 'TPCO'


      IF(FIRST) THEN
c gab: dummy code to avoild compile-time warnings errors
       GMM=GMM
       ALM=ALM
       XMID=XMID
       YMID=YMID
       BEM=BEM
       CU2=CU2
       GGRAD(1)=GGRAD(1)
       PC(1)=PC(1)

          PI    = 2.0*ASIN(1.0)
          PIO2  = 0.5*PI
          PIT2  = 2.0*PI
C GAB:
C          NTPCO = NAMIND('TPCO')
          FIRST =.FALSE.
      END IF
      DO 1   I=1,15
  1     EE0(I)=0.0
      DO 2   I=1,5
        GRAD(I)=0.0
        VV0(I) =0.0
  2   CONTINUE
      CHI2=0.0
      N   = NPT
      IF(N.GT.MPT) RETURN
      IF(N.LT.3)   RETURN
C
C-----> INPUT DATA
C
C GAB
C          KTPCO=IW(NTPCO)
      DO 10 I=1,N
c GAB
c        KSTRT = KKTPCO(LIST(I))
C        XF(I)  = RW(KSTRT+2)*COS(RW(KSTRT+3))
C        YF(I)  = RW(KSTRT+2)*SIN(RW(KSTRT+3))
C        RF(I)  = RW(KSTRT+2)
C        PF(I)  = RW(KSTRT+3)
C        WF(I)  = (RW(KSTRT+5)+0.000000001)**(-1)
C        SP2(I) = WF(I)*(RF(I)*RF(I))
C        ZF(I)  = RW(KSTRT+4)
C        WZF(I) = 1.0/(RW(KSTRT+6)+0.000001)
        IHIT=LIST(I)
        XF(I)  = DBLE(TPCO_XV(IHIT))
        YF(I)  = DBLE(TPCO_YV(IHIT))
        RF(I)  = TPCO_RV(IHIT)
        PF(I)  = TPCO_PH(IHIT)
        WF(I)  = (DBLE(TPCO_SR(IHIT))+0.000000001D0)**(-1)
        SP2(I) = WF(I)*(RF(I)*RF(I))
        ZF(I)  = TPCO_ZV(IHIT)
        WZF(I) = 1.0/(TPCO_SZ(IHIT)+0.000001)
   10 CONTINUE

C(KH)--------removed
CC GAB: Allow for VXD hits to be included in the fit:
C      DO I=1,NHVX
C       N=N+1
C       XF(N)  = DBLE(XH(I))
C       YF(N)  = DBLE(YH(I))
C       RF(N)  = SQRT(XH(I)**2+YH(I)**2)
C       PF(N)  = PHI(I)
C       WF(N)  = (DBLE(W1(I)**2)+0.000000001D0)**(-1)
C       SP2(N) = WF(I)*(RF(I)*RF(I))
C       ZF(N)  = ZH(I)
C       WZF(N) = 1.0/(W2(I)**2+0.000001)
C      ENDDO
C(KH)--------end

C***************************************************************C
C                                                               C
C           CIRCLE FIT IN THE [X,Y] PLANE                       C
C           =============================                       C
C                                                               C
C***************************************************************C
      WSUM= 0.0
      RSS = 0.0
      PRO = 0.0
      XM = 0.
      YM = 0.
      WN=0.0
      DO 100 I= 1, N
        XM = XM + XF(I)*WF(I)
        YM = YM + YF(I)*WF(I)
        WN = WN + WF(I)
  100 CONTINUE
      RN = 1.D0/WN
C **
      XM = XM * RN
      YM = YM * RN
      X2 = 0.
      Y2 = 0.
      XY = 0.
      XD = 0.
      YD = 0.
      D2 = 0.
      DO 102 I= 1, N
         XI = XF(I) - XM
         YI = YF(I) - YM
         XX = XI**2
         YY = YI**2
         X2 = X2 + XX*WF(I)
         Y2 = Y2 + YY*WF(I)
         XY = XY + XI*YI*WF(I)
         DD = XX + YY
         XD = XD + XI*DD*WF(I)
         YD = YD + YI*DD*WF(I)
         D2 = D2 + DD**2*WF(I)
  102 CONTINUE
C **
      X2 = X2*RN
      Y2 = Y2*RN
      XY = XY*RN
      D2 = D2*RN
      XD = XD*RN
      YD = YD*RN
      F = 3.D0*X2 + Y2
      G = 3.D0*Y2 + X2
      FG = F*G
      H = XY + XY
      H2 = H**2
      P2 = XD**2
      Q2 = YD**2
      GAM0 = X2 + Y2
         FACT = GAM0**2
      A2 = (FG-H2-D2)/FACT
         FACT = FACT*GAM0
      A1 = (D2*(F+G) - 2.D0*(P2+Q2))/FACT
         FACT = FACT*GAM0
      A0 = (D2*(H2-FG) + 2.D0*(P2*G + Q2*F) - 4.D0*XD*YD*H)/FACT
      A22 = A2 + A2
      YB = 1.0E30
      ITER = 0
      XA = 1.D0
C **                MAIN ITERATION
  103 YA = A0 + XA*(A1 + XA*(A2 + XA*(XA-4.D0)))
      IF (ITER .GE. ITMAX)                      GO TO 105
      DY = A1 + XA*(A22 + XA*(4.D0*XA - 12.D0))
      XB = XA - YA/DY
      IF (ABS(YA).GT.ABS(YB)) XB=0.5D0*(XB+XA)
      IF (ABS(XA-XB) .LT. EPS)                  GO TO 105
      XA = XB
      YB = YA
      ITER = ITER + 1
      GO TO 103
C **
  105 CONTINUE
      ROOT = XB
      GAM = GAM0*XB
      F1 = F - GAM
      G1 = G - GAM
      X1 = XD*G1 - YD*H
      Y1 = YD*F1 - XD*H
      DET = F1*G1 - H2
      DEN2= 1.D0/(X1**2 + Y1**2 + GAM*DET**2)
      IF(DEN2.LE.0.D0)                GO TO 999
      DEN = DSQRT(DEN2)
      CUR = DET*DEN                  + 0.0000000001D0
      ALF = -(XM*DET + X1)*DEN
      BET = -(YM*DET + Y1)*DEN
      RM = XM**2 + YM**2
      GAM = ((RM-GAM)*DET + 2.D0*(XM*X1 + YM*Y1))*DEN*0.5D0
C
C--------> CALCULATION OF STANDARD CIRCLE PARAMETERS
C          NB: CUR IS ALWAYS POSITIVE
      RR0=CUR
      ASYM = BET*XM-ALF*YM
      SST = 1.0
      IF(ASYM.LT.0.0) SST=-1.0
      RR0 = SST*CUR
      IF((ALF*ALF+BET*BET).LE.0.D0)              GO TO 999
      SA2B2 = 1.D0/DSQRT(ALF*ALF+BET*BET)
      DD0 = (1.D0-1.D0/SA2B2)/CUR
      VV0(1) = RR0
      VV0(4) = DD0
      CHECK=SST*RR0*DD0
      IF(CHECK.EQ.1.) THEN
        DD0=DD0-.007
        VV0(4)=DD0
      ENDIF
C***************************************************************C
C                                                               C
C           STRAIGHT LINE FIT IN THE [S,Z] PLANE                C
C           ====================================                C
C                                                               C
C***************************************************************C
      DO 120 I=1,N
        EEE(I) = 0.5*VV0(1)
     1          *SQRT(ABS( (RF(I)*RF(I)-VV0(4)*VV0(4))
     2                /(1.0-SST*VV0(1)*VV0(4))     ) )
        SXY(I)=2.0*EEE(I)*(1.0+0.166667*EEE(I)*EEE(I))/RR0
CD      IF(EEE(I).GT.0.99990)  EEE(I)= 0.99990
CD      IF(EEE(I).LT.-0.99990) EEE(I)= -0.99990
CD      SXY(I)=2.0*ASIN(EEE(I))/RR0
  120 CONTINUE
      SUMS  = 0.0
      SUMSS = 0.0
      SUMZ  = 0.0
      SUMSZ = 0.0
      SUMW  = 0.0
      DO 130 I=1,N
        SUMW  = SUMW  +                 WZF(I)
        SUMS  = SUMS  + SXY(I)        * WZF(I)
        SUMSS = SUMSS + SXY(I)*SXY(I) * WZF(I)
        SUMZ  = SUMZ  + ZF(I)         * WZF(I)
        SUMSZ = SUMSZ + ZF(I)*SXY(I)  * WZF(I)
  130 CONTINUE
      DENOM = SUMW*SUMSS - SUMS*SUMS
      IF( DENOM .EQ. 0.0 ) GOTO 999
      DZDS  = (SUMW*SUMSZ-SUMS*SUMZ) /DENOM
      VV0(2)= DZDS
C***************************************************************C
C                                                               C
C              MULTIPLE SCATTERING ERRORS                       C
C              ==========================                       C
C                                                               C
C***************************************************************C
       PSI0  = SCACO * COS( ATAN(DZDS) ) * CUR
       XL = 0.0
       DO 140 I=2,N
         DDR = RF(I) - RF(I-1)
         IF(DDR.EQ.0.) GOTO 140
         DDX = XF(I) - XF(I-1)
         DDY = YF(I) - YF(I-1)
         DDZ = ZF(I) - ZF(I-1)
         XL  = XL + SQRT( DDX*DDX + DDY*DDY + DDZ*DDZ )
         SNA = DDR/SQRT( DDX*DDX + DDY*DDY )
         SNB = DDR/SQRT( DDR*DDR + DDZ*DDZ )
         WF(I)  = WF(I) / (1.D0 + WF(I)*XL*(XL*PSI0/SNA)**2 )
         SP2(I) = WF(I)*RF(I)*RF(I)
         WZF(I) = WZF(I) / (1.D0 + WZF(I)*XL*(XL*PSI0/SNB)**2 )
  140 CONTINUE
C***************************************************************C
C                                                               C
C          SECOND CIRCLE FIT IN THE [X,Y] PLANE                 C
C          ====================================                 C
C                                                               C
C***************************************************************C
      WSUM= 0.0
      RSS = 0.0
      PRO = 0.0
      XM = 0.
      YM = 0.
      WN=0.0
      DO 201 I= 1, N
        XM = XM + XF(I)*WF(I)
        YM = YM + YF(I)*WF(I)
        WN = WN + WF(I)
  201 CONTINUE
      RN = 1.D0/WN
C **
      XM = XM * RN
      YM = YM * RN
      X2 = 0.
      Y2 = 0.
      XY = 0.
      XD = 0.
      YD = 0.
      D2 = 0.
      DO 202 I= 1, N
         XI = XF(I) - XM
         YI = YF(I) - YM
         XX = XI**2
         YY = YI**2
         X2 = X2 + XX*WF(I)
         Y2 = Y2 + YY*WF(I)
         XY = XY + XI*YI*WF(I)
         DD = XX + YY
         XD = XD + XI*DD*WF(I)
         YD = YD + YI*DD*WF(I)
         D2 = D2 + DD**2*WF(I)
  202 CONTINUE
C **
      X2 = X2*RN
      Y2 = Y2*RN
      XY = XY*RN
      D2 = D2*RN
      XD = XD*RN
      YD = YD*RN
      F = 3.D0*X2 + Y2
      G = 3.D0*Y2 + X2
      FG = F*G
      H = XY + XY
      H2 = H**2
      P2 = XD**2
      Q2 = YD**2
      GAM0 = X2 + Y2
         FACT = GAM0**2
      A2 = (FG-H2-D2)/FACT
         FACT = FACT*GAM0
      A1 = (D2*(F+G) - 2.D0*(P2+Q2))/FACT
         FACT = FACT*GAM0
      A0 = (D2*(H2-FG) + 2.D0*(P2*G + Q2*F) - 4.D0*XD*YD*H)/FACT
      A22 = A2 + A2
      YB = 1.0E30
      ITER = 0
      XA = 1.D0
C **                MAIN ITERATION
  203 YA = A0 + XA*(A1 + XA*(A2 + XA*(XA-4.D0)))
      IF (ITER .GE. ITMAX)                      GO TO 205
      DY = A1 + XA*(A22 + XA*(4.D0*XA - 12.D0))
      XB = XA - YA/DY
      IF (ABS(YA).GT.ABS(YB)) XB=0.5D0*(XB+XA)
      IF (ABS(XA-XB) .LT. EPS)                  GO TO 205
      XA = XB
      YB = YA
      ITER = ITER + 1
      GO TO 203
C **
  205 CONTINUE
      ROOT = XB
      GAM = GAM0*XB
      F1 = F - GAM
      G1 = G - GAM
      X1 = XD*G1 - YD*H
      Y1 = YD*F1 - XD*H
      DET = F1*G1 - H2
      DEN2= 1.D0/(X1**2 + Y1**2 + GAM*DET**2)
      IF(DEN2.LE.0.D0)                GO TO 999
      DEN = DSQRT(DEN2)
      CUR = DET*DEN                  + 0.0000000001D0
      ALF = -(XM*DET + X1)*DEN
      BET = -(YM*DET + Y1)*DEN
      RM = XM**2 + YM**2
      GAM = ((RM-GAM)*DET + 2.D0*(XM*X1 + YM*Y1))*DEN*0.5D0
C
C--------> CALCULATION OF STANDARD CIRCLE PARAMETERS
C          NB: CUR IS ALWAYS POSITIVE
C
      RR0=CUR
      ASYM = BET*XM-ALF*YM
      SST = 1.0
      IF(ASYM.LT.0.0) SST=-1.0
      RR0 = SST*CUR
      IF((ALF*ALF+BET*BET).LE.0.D0)              GO TO 999

      SA2B2 = 1.D0/DSQRT(ALF*ALF+BET*BET)
      DD0 = (1.D0-1.D0/SA2B2)/CUR
      PHIC = DASIN(ALF*SA2B2)+PIO2
      IF(BET.GT.0)    PHIC=PIT2-PHIC
      PH0 = PHIC+PIO2
      IF(RR0.LE.0)    PH0=PH0-PI
      IF(PH0.GT.PIT2) PH0=PH0-PIT2
      IF(PH0.LT.0.0)  PH0=PH0+PIT2
      VV0(1) = RR0
      VV0(3) = PH0
      VV0(4) = DD0
      CHECK=SST*RR0*DD0
      IF(CHECK.EQ.1.) THEN
        DD0=DD0-.007
        VV0(4)=DD0
      ENDIF
C
C-----> CALCULATE PHI DISTANCES TO MEASURED POINTS
C
      AA0 =SST
      OME =RR0
      GG0 = OME*DD0-AA0
      HH0 = 0.0000000001
      IF(ABS(GG0).LT.(1.0/HH0)) HH0=1.0/GG0
      HH0=1.0/GG0
      DO 210 I=1,N
        ASYM   = BET*XF(I)-ALF*YF(I)
        SS0(I) =1.0
        IF(ASYM.LT.0.0) SS0(I)=-1.0
        FF0   = OME*(RF(I)*RF(I)-DD0*DD0)/(2.0*RF(I)*GG0) + DD0/RF(I)
C
        IF(FF0.LT.-1.0) FF0 = -1.0
        IF(FF0.GT.1.0)  FF0 = 1.0
C
        DEL(I)= PH0 + (SS0(I)-AA0)*PIO2 + SS0(I)*ASIN(FF0) - PF(I)
        IF(DEL(I).GT.PI) DEL(I)=DEL(I)-PIT2
        IF(DEL(I).LT.-PI)DEL(I)=DEL(I)+PIT2
  210 CONTINUE
C***************************************************************C
C                                                               C
C          SECOND SRAIGHT LINE FIT IN THE [S,Z] PLANE           C
C          ==========================================           C
C                                                               C
C***************************************************************C
      DO 220 I=1,N
        EEE(I) = 0.5*VV0(1)
     1          *SQRT(ABS( (RF(I)*RF(I)-VV0(4)*VV0(4))
     2                /(1.0-AA0*VV0(1)*VV0(4))     ) )
C
        IF(EEE(I).GT.0.99990)  EEE(I)= 0.99990
        IF(EEE(I).LT.-0.99990) EEE(I)= -0.99990
C
        SXY(I)=2.0*ASIN(EEE(I))/OME
  220 CONTINUE
      SUMS  = 0.0
      SUMSS = 0.0
      SUMZ  = 0.0
      SUMZZ = 0.0
      SUMSZ = 0.0
      SUMW  = 0.0
      DO 230 I=1,N
        SUMW  = SUMW  +                 WZF(I)
        SUMS  = SUMS  + SXY(I)        * WZF(I)
        SUMSS = SUMSS + SXY(I)*SXY(I) * WZF(I)
        SUMZ  = SUMZ  + ZF(I)         * WZF(I)
        SUMZZ = SUMZZ + ZF(I)*ZF(I)   * WZF(I)
        SUMSZ = SUMSZ + ZF(I)*SXY(I)  * WZF(I)
  230 CONTINUE
      DENOM = SUMW*SUMSS - SUMS*SUMS
      IF( DENOM .EQ. 0.0 ) GOTO 999
      DZDS  = (SUMW*SUMSZ-SUMS*SUMZ) /DENOM
      ZZ0   = (SUMSS*SUMZ-SUMS*SUMSZ)/DENOM
      VV0(2)= DZDS
      VV0(5)= ZZ0
C
C-----> CALCULATION CHI**2
C
      DO 370 I=1,N
        DELZ(I)= ZZ0+DZDS*SXY(I)-ZF(I)
        CHI2   = CHI2  + SP2(I)*DEL(I)*DEL(I)
     1                 + WZF(I)*DELZ(I)*DELZ(I)
  370 CONTINUE
      IF(IOPT.LT.2)                          RETURN
CCC   IF(CHI2.GT.5000.0)                     RETURN
C***************************************************************C
C                                                               C
C                     ERROR MATRIX                              C
C                     ============                              C
C                                                               C
C***************************************************************C
      DO 380 I=1,N
        FF0 = OME*(RF(I)*RF(I)-DD0*DD0)/(2.0*RF(I)*GG0) + DD0/RF(I)
        IF(FF0.GT. 0.99990)  FF0=  0.99990
        IF(FF0.LT.-0.99990)  FF0= -0.99990
        ETA = SS0(I)/SQRT(ABS((1.0+FF0)*(1.0-FF0)))
        DFD = (1.0+HH0*HH0*(1.0-OME*OME*RF(I)*RF(I)))/(2.0*RF(I))
        DFO = -AA0*(RF(I)*RF(I)-DD0*DD0)*HH0*HH0/(2.0*RF(I))
        DPD = ETA*DFD
        DPO = ETA*DFO
C-----> DERIVATIVES OF Z COMPONENT
        GGG = EEE(I)/SQRT(ABS( (1.0+EEE(I))*(1.0-EEE(I))))
        DZA = SXY(I)
        CHECK=RF(I)*RF(I)-VV0(4)*VV0(4)
        IF(CHECK.EQ.0.) CHECK=2.*.007
        DZD = 2.0*( VV0(2)/VV0(1) ) * ABS( GGG )
     1       *(  0.5*AA0*VV0(1)/( 1.0-AA0*VV0(4)*VV0(1) )
     2           -VV0(4)/CHECK    )
        DZO = -VV0(2)*SXY(I)/VV0(1)
     1        +VV0(2)* GGG/( VV0(1)*VV0(1))
     2         *( 2.0+ AA0*VV0(1)*VV0(4)/(1.0-AA0*VV0(1)*VV0(4)) )
C-----> ERROR MARTIX
        EE0(1) = EE0(1) + SP2(I)*  DPO*DPO  + WZF(I) * DZO*DZO
        EE0(2) = EE0(2)                     + WZF(I) * DZA*DZO
        EE0(3) = EE0(3)                     + WZF(I) * DZA*DZA
        EE0(4) = EE0(4) + SP2(I)*  DPO
        EE0(5) = 0.0
        EE0(6) = EE0(6) + SP2(I)
        EE0(7) = EE0(7) + SP2(I)*  DPO*DPD  + WZF(I) * DZO*DZD
        EE0(8) = EE0(8)                     + WZF(I) * DZA*DZD
        EE0(9) = EE0(9) + SP2(I)*      DPD
        EE0(10)= EE0(10)+ SP2(I)*  DPD*DPD  + WZF(I) * DZD*DZD
        EE0(11)= EE0(11)                    + WZF(I) * DZO
        EE0(12)= EE0(12)                    + WZF(I) * DZA
        EE0(13)= 0.0
        EE0(14)= EE0(14)                    + WZF(I) * DZD
        EE0(15)= EE0(15)                    + WZF(I)
C-----> GRADIENT VECTOR
        GRAD(1)=GRAD(1) - DEL(I) *SP2(I)*DPO - DELZ(I)*WZF(I)*DZO
        GRAD(2)=GRAD(2) -                      DELZ(I)*WZF(I)*DZA
        GRAD(3)=GRAD(3) - DEL(I) *SP2(I)
        GRAD(4)=GRAD(4) - DEL(I) *SP2(I)*DPD - DELZ(I)*WZF(I)*DZD
        GRAD(5)=GRAD(5) -                      DELZ(I)*WZF(I)
  380 CONTINUE
C
      IF(IOPT.LT.3)                          RETURN
C***************************************************************C
C                                                               C
C         NEWTONIAN ITERATION IN 5 PARAMETERS                   C
C         ===================================                   C
C                                                               C
C***************************************************************C
      DO 401 I=1,15
  401   COV(I)=EE0(I)
      CALL SMINV(COV,VV1,5,0,NRANK)
      CALL SMAV(DV,COV,GRAD,5)
      DO 402 I=1,5
  402   VV1(I)=VV0(I)+DV(I)
C------>  NEW DIFFERENCES IN PHI AND Z
      GG0 = VV1(1)*VV1(4)-AA0
      DO 410 I=1,N
        FF0    = VV1(1)*(RF(I)*RF(I)-VV1(4)**2)
     1          /(2.0*RF(I)*GG0) + VV1(4)/RF(I)
C
        IF(FF0.GT.1.0)  FF0 = 1.0
        IF(FF0.LT.-1.0) FF0 = -1.0
C
        DELN(I) = VV1(3) + (SS0(I)-AA0)*PIO2+SS0(I)*ASIN(FF0)-PF(I)
        IF(DELN(I).GT.PI) DELN(I)=DELN(I)-PIT2
        IF(DELN(I).LT.-PI)DELN(I)=DELN(I)+PIT2
        EEE(I) = 0.5*VV1(1)
     1          *SQRT(ABS( (RF(I)*RF(I)-VV1(4)*VV1(4))
     2                /(1.0-AA0*VV1(1)*VV1(4))     ))
        IF(EEE(I).GT.0.99990)  EEE(I)= 0.99990
        IF(EEE(I).LT.-0.99990) EEE(I)= -0.99990
        SXY(I) = 2.0*ASIN(EEE(I))/VV1(1)
        DELZN(I)= VV1(5)+VV1(2)*SXY(I)-ZF(I)
  410 CONTINUE
C-----> CALCULATION CHI**2
      CHI1 = 0.0
      DO 420 I=1,N
        CHI1   = CHI1  + SP2(I)*DELN(I)*DELN(I)
     1                 + WZF(I)*DELZN(I)*DELZN(I)
  420 CONTINUE
      IF(CHI1.LT.CHI2) THEN
         DO 430 I=1,5
  430    VV0(I) = VV1(I)
         CHI2 = CHI1
      END IF
      RETURN
  999 CONTINUE
C------> JUMP HERE IF SOMETHING GOES CRAZY
      CHI2 = 1.0E30
      RETURN
      ENTRY TMSHED(DRPHI,DZ)
C-----**********************-------------------------------------C
C      OUTPUTS THE DIFFERENCES IN R*PHI AND Z FOR THE            C
C      COORDINATES FITTED BY THE LAST CALL OF TFTHEL  [CM]       C
C----------------------------------------------------------------C
      DO 500 I=1,N
        DRPHI(I) = RF(I)*DEL(I)
  500   DZ(I)    = DELZ(I)
      RETURN
      ENTRY TMSHEC(DCHI2)
C-----**********************-------------------------------------C
C      OUTPUTS THE CHI**2 CONTRIBUTIONS OF ALL POINTS            C
C----------------------------------------------------------------C
      DO 600 I=1,N
        DCHI2(I)  = SP2(I)*DEL(I)*DEL(I)
     1             +WZF(I)*DELZ(I)*DELZ(I)
  600 CONTINUE
      RETURN
      ENTRY TMSCHI(DRPHI,DZ)
C-----**********************-------------------------------------C
C      OUTPUTS THE DIFFERENCES IN R*PHI AND Z NORMALISED         C
C      TO THE ERRORS ON COORDINATES + MULTIPLE SCATTERING  OF    C
C      ALL COORDINATES FITTED BY THE LAST CALL OF TFTHEL  [CM]   C
C----------------------------------------------------------------C
      DO 700 I=1,N
        DRPHI(I) = DEL(I)  * SQRT(SP2(I))
  700   DZ(I)    = DELZ(I) * SQRT(WZF(I))
      RETURN
      END
c---------------------------------------------------------------------
