*  
* $Id: trslnk.F,v 1.6 2006-02-24 16:36:55 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.5  2005/08/08 07:09:13  aplin
* Made f77 tracking code use GEAR to define the geomtery of the TPC. LTPDRO now defines the maximum number of rows is used to define the size of arrays, this is limited to 224 due the use of 7 '32 bit' bit registers in trkfnd.F increased, though at present it is not likely that anybody would want more. The number of TPC padrows is defined at run time by NRTPC which should of course not exceed LTPDRO, although this is checked and the programe exits with a verbose error message. A wrapper function gettpcgeom is used to pass the GEAR TPC parameters from C++ to f77. MarlinUtil/include/marlin_tpcgeom.h have MarlinUtil/src/marlin_tpcgeom.cc consequently been removed as they are no longer needed.
*
* Revision 1.4  2005/08/04 11:18:31  aplin
* corrected comment
*
* Revision 1.3  2005/08/04 11:09:41  aplin
* added comment
*
* Revision 1.2  2005/08/04 11:08:42  aplin
* drop event for to many link, a more ellegant solution is needed
*
* Revision 1.1.1.1  2005/07/11 15:15:48  aplin
* Initial version
*
* Revision 1.2  2005/06/20 18:49:46  aplin
* change tpc geometery so that it is read from  same header file as cpp i.e. marlin_tpcgeom.h
*
* Revision 1.1.1.1  2005/06/18 10:05:35  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      SUBROUTINE TRSLNK(IR3,IR2,IR1,IPTR3,NLK,IER)
      IMPLICIT NONE
C****************************************************************************
c* Version number: 0.0      Date of last modification  26.08.99             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
C* Main modifiers from ALEPH to BRAHMS code:                                *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
C
C ---------------------------------------------------------------------
C! Search for link across 3  padrows
C!
C!  Author:   M. Mermikides   6-10-86
C!  Modified: R. Johnson     19-01-88
C!
C!  Input:
C!      - IR3       /I      TPC row number
C!      - IPTR3     /I      Hit in this row at which to begin.  This
C!                          points into the workbank INDBIN
C!      - IR2       /I      Middle row number
C!      - IR1       /I      Inner row number
C!  Output:
C!      - NLK       /I      No of valid links found
C!      - IER       /I      Error code =0 for good return
C!                          Set to 1 if BOS runs out of space
C!                          Set to -1 if 
C!  Called by TRKFND
C!
C!  Description:
C!  -----------
C!  Condition for accepting a link are:
C!        1. Phi difference between hits on successive padrows < DPNEXT
C!        2. z difference between hits on successive padrows   < DZNEXT
C!        3. Radius of curvature of link  > RHOMIN
C!        4. Consistency in z (|pred.z -msd z|) < TDZTOL
C!        5. Closest approach to origin in r-phi, D0 < TD0TOL
C!        6. Z at D0  < TZ0TOL
C!  Rows, hit ref and link parameters are entered in workbanks INDLNK
C!---------------------------------------------------------------------
C
#include "include/fndparms.inc"
#include "include/alcons.inc"
#include "include/tpgeom.inc"
c#include "marlin_tpcgeom.h"
#include "include/trfdct.inc"
#include "include/tcbins.inc"
#include "include/tpco.inc"
#include "include/tpcdims.inc"
#include "include/my_jndb.inc"
#include "include/my_link.inc"
C
C
C----------------------------------------------------------------------
c GAB:
      INTEGER IR1,IR2,IR3,IPTR1,IPTR2,IPTR3,NLK,IER,JJ
      INTEGER NXLNK,IC1,IC2,IC3,IBST1,IBST2,II,IBIN1,IBIN2
      REAL PADFAC,PHI0,PHI1,PHI2,PHI3,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3
      REAL XC,YC,A1,B1,C1,A2,B2,C2,DET,RHOSQ,RHO,CH10,CH12,CH13
      REAL DFI21,DFI32,SGZ2,Z2PR,ZDIFF,CHI2Z,ARG,RC,PHIC,CP,SP
      REAL Y1R,Y2R,Y3R,QQ,BETA,COSB,SINB,DW2,S12,S13,TANL
      REAL X01,X02,Y01,Y02,D1SQ,D2SQ,D0,X0,Y0,Z0,PSI0


      REAL TRPHDI
      EXTERNAL TRPHDI

      LOGICAL FIRST
      DATA FIRST/.TRUE./
      DATA  NXLNK/100/

c#include "bkgmode.inc"
c#include "gcbank.inc"

#include "include/tpcfunc2.inc"

      SAVE
#include "include/tpcfunc.inc"

C----------------------------------------------------------------------
      IF (FIRST) THEN
        PADFAC=TPCPADR**2/12.
        FIRST=.FALSE.
c        NTPCO=NAMIND('TPCO')
c        NTCRL=NAMIND('TCRL')
      ENDIF
C
      IER=0
c      KTPCO = IW(NTPCO)
c      KTCRL = IW(NTCRL)
      NLK = 0


C GAB:
C      IC3= ITABL(INDBIN,IPTR3,JNDBCO)
C     I3= IC3 - ITABL(KTCRL,IR3,JTCROC)
      IC3= MY_JNDBCO(IPTR3)
C      PHI3  = RTABL(KTPCO,IC3,JTPCPH)
C      PHI3=TPCPHI(IC3)
      PHI3=TPCO_PH(IC3)

C GAB:
C      X3 = RTABL(INDCXY,IC3,1)
C      Y3 = RTABL(INDCXY,IC3,2)
C      Z3 = RTABL(KTPCO,IC3,JTPCZV)
C      X3 = TPCHITS(1,IC3)
C      Y3 = TPCHITS(2,IC3)
C      Z3 = TPCHITS(3,IC3)
      X3 = TPCO_XV(IC3)
      Y3 = TPCO_YV(IC3)
      Z3 = TPCO_ZV(IC3)
C
C  Hits on middle padrow.  First find the correct phi bins.
C
      IBST2= IXWRAP(INT(TCOFFS+PHI3/TC2DEL)) + 1

      DO 20 II=IBST2,IBST2+4*NMXWRP,2
         IBIN2= IXWRAP(II)
C
C  Loop over unused hits in this bin, using the
C  linked list in bank INDBIN.
C
         IPTR2=IBNPTR(IBIN2,IR2)
c(kh)  212    IF (IPTR2.EQ.0) GO TO 213
  212    IF (IPTR2.LE.0) GO TO 213
c GAB:
c            IC2= ITABL(INDBIN,IPTR2,JNDBCO)
C            PHI2 = RTABL(KTPCO,IC2,JTPCPH)
            IC2= MY_JNDBCO(IPTR2)
C            PHI2=TPCPHI(IC2)
            PHI2=TPCO_PH(IC2)
            DFI32 = TRPHDI(PHI3,PHI2)
C
C Test on max phi and z diff.
C
            IF (ABS(DFI32).GT.DPNEXT) THEN
              IPTR2=MY_JNDBFW(IPTR2)
              GO TO 212
            ENDIF
C GAB:
C            Z2 = RTABL(KTPCO,IC2,JTPCZV)
C            Z2 = TPCHITS(3,IC2)
            Z2 = TPCO_ZV(IC2)
            IF (ABS(Z2-Z3).GT.DZNEXT) THEN
              IPTR2=MY_JNDBFW(IPTR2)
              GO TO 212
            ENDIF

C GAB:
C            SGZ2 = RTABL(KTPCO,IC2,JTPCSZ)
C            X2= RTABL(INDCXY,IC2,1)
C            Y2= RTABL(INDCXY,IC2,2)
C GAB: I2 Doesn't seem to be used again, so ignore it ...
C            I2= IC2 - ITABL(KTCRL,IR2,JTCROC)
C            SGZ2 = (TPCHITS(7,IC2))**2
C            X2= TPCHITS(1,IC2)
C            Y2= TPCHITS(2,IC2)
            SGZ2 = TPCO_SZ(IC2)
            X2= TPCO_XV(IC2)
            Y2= TPCO_YV(IC2)
C
C Hits on inner layer.  First find the correct bin.
C

            IBST1= IXWRAP(INT(TCOFFS+PHI2/TC2DEL)) + 1
            DO 10 JJ=IBST1,IBST1+4*NMXWRP,2
               IBIN1= IXWRAP(JJ)
C
C  Loop over unused hits in this bin, using the
C  linked list in bank INDBIN.
C
               IPTR1=IBNPTR(IBIN1,IR1)
c(kh) 112          IF (IPTR1.EQ.0) GO TO 113
  112          IF (IPTR1.LE.0) GO TO 113

C GAB:
C                  IC1= ITABL(INDBIN,IPTR1,JNDBCO)
                  IC1= MY_JNDBCO(IPTR1)

C                  PHI1 = RTABL(KTPCO,IC1,JTPCPH)
C                  PHI1=TPCPHI(IC1)
                  PHI1=TPCO_PH(IC1)

                  DFI21 = TRPHDI(PHI2,PHI1)
                  IF (ABS(DFI21).GT.DPNEXT) THEN
C GAB:
c                    IPTR1=ITABL(INDBIN,IPTR1,JNDBFW)
                    IPTR1=MY_JNDBFW(IPTR1)
                    GO TO 112
                  ENDIF

C GAB:
C                  Z1 = RTABL(KTPCO,IC1,JTPCZV)
C                  Z1 = TPCHITS(3,IC1)
                  Z1 = TPCO_ZV(IC1)
                  IF (ABS(Z2-Z1).GT.DZNEXT) THEN
C GAB:
C                    IPTR1=ITABL(INDBIN,IPTR1,JNDBFW)
                    IPTR1=MY_JNDBFW(IPTR1)
                    GO TO 112
                  ENDIF
C
C GAB:
C                  X1 = RTABL(INDCXY,IC1,1)
C                  Y1 = RTABL(INDCXY,IC1,2)
C                  X1 = TPCHITS(1,IC1)
C                  Y1 = TPCHITS(2,IC1)
                  X1 = TPCO_XV(IC1)
                  Y1 = TPCO_YV(IC1)
                  A1 = X2-X1
                  B1 = Y2-Y1
                  C1 = 0.5*(X1*X1-X2*X2+Y1*Y1-Y2*Y2)
                  A2 = X3-X2
                  B2 = Y3-Y2
                  C2 = 0.5*(X2*X2-X3*X3+Y2*Y2-Y3*Y3)
                  DET = A1*B2-A2*B1
C
C Centre of track orbit; protect against 3 perfectly aligned points
C
                  IF (DET.EQ.0.) THEN
C GAB:
C                    IPTR1=ITABL(INDBIN,IPTR1,JNDBFW)
                    IPTR1=MY_JNDBFW(IPTR1)
                    GO TO 112
                  ENDIF
                  XC = (B1*C2-B2*C1)/DET
                  YC =-(A1*C2-A2*C1)/DET
                  RHOSQ = (X1-XC)**2 + (Y1-YC)**2

C
C  Check on radius of curvature
                  IF (RHOSQ.LT.RHOMIN**2) THEN
C GAB:
C                    IPTR1=ITABL(INDBIN,IPTR1,JNDBFW)
                    IPTR1=MY_JNDBFW(IPTR1)
                    GO TO 112
                  ENDIF
                  RHO = SQRT(RHOSQ)
                  CH12 = SQRT(A1*A1+B1*B1)
                  ARG = CH12/(2.*RHO)
                  IF (ARG.GE.1.0) THEN

C GAB:
C                    IPTR1=ITABL(INDBIN,IPTR1,JNDBFW)
                    IPTR1=MY_JNDBFW(IPTR1)
                    GO TO 112
                  ENDIF

C
C  If chord length is less than 15 cm and radius > 20 cm approximate
C  path length. The error is less than 1 mm under these conditions.
C
                  IF (CH12.LT.15. .AND. RHO. GT.20.) THEN
                     S12 = CH12*(1. + ARG*ARG/8.)
                  ELSE
                     S12 = 2.*RHO*ASIN(ARG)
                  ENDIF

C
                  CH13 = SQRT((X1-X3)**2+(Y1-Y3)**2)
                  ARG = CH13/(2.*RHO)
                  IF (ARG.GE.1.0) THEN
C GAB:
C                    IPTR1=ITABL(INDBIN,IPTR1,JNDBFW)
                    IPTR1=MY_JNDBFW(IPTR1)
                    GO TO 112
                  ENDIF
                  IF (CH13.LT.15. .AND. RHO. GT.20.) THEN
                     S13 = CH13*(1. + ARG*ARG/8.)
                  ELSE
                     S13 = 2.*RHO*ASIN(ARG)
                  ENDIF
C
C Get tangent of dip angle and predict z of middle point
C
                  TANL = (Z3-Z1)/S13
                  Z2PR = Z1 + S12*TANL

C
C Apply test on difference of Z2 from predicted Z
C
                  ZDIFF = Z2PR - Z2
                  CHI2Z = ZDIFF*ZDIFF/SGZ2

                  IF (CHI2Z.GT.TDZTOL) THEN
C GAB:
C                    IPTR1=ITABL(INDBIN,IPTR1,JNDBFW)
                    IPTR1=MY_JNDBFW(IPTR1)
                    GO TO 112
                  ENDIF
C
                  RC = SQRT(XC*XC + YC*YC)
                  X01 = XC*(1. + RHO/RC)
                  X02 = XC*(1. - RHO/RC)

C
C Protect against circle centre on Y-axis
C
                  IF (ABS(XC).LT.0.00001) THEN
                     Y01 = YC*(1. + RHO/RC)
                     Y02 = YC*(1. - RHO/RC)
                  ELSE
                     Y01 = YC*X01/XC
                     Y02 = YC*X02/XC
                  ENDIF

C
C   Find closest point to origin (D0)
C   D0 is positive if track encloses origin, otherwise negative
C
                  D1SQ = X01*X01 + Y01*Y01
                  D2SQ = X02*X02 + Y02*Y02
                  IF (D1SQ.LT.D2SQ) THEN
                     D0=SQRT(D1SQ)
                     X0 = X01
                     Y0 = Y01
                  ELSE
                     D0=SQRT(D2SQ)
                     X0 = X02
                     Y0 = Y02
                  ENDIF
C
C  Cut on D0 and Z0
C
                  IF (D0.GT.TD0TOL) THEN
C GAB:
C                    IPTR1=ITABL(INDBIN,IPTR1,JNDBFW)
                    IPTR1=MY_JNDBFW(IPTR1)
                    GO TO 112
                  ENDIF
                  IF (RC.GT.RHO)   D0 = -D0
                  CH10 = SQRT((X1-X0)**2 + (Y1-Y0)**2)
                  PSI0 = 2.*ASIN(CH10/(2.*RHO))
                  Z0 = Z1 - RHO*PSI0*TANL

                  IF (ABS(Z0).GT.TZ0TOL) THEN
C GAB:
C                    IPTR1=ITABL(INDBIN,IPTR1,JNDBFW)
                    IPTR1=MY_JNDBFW(IPTR1)
                    GO TO 112
                  ENDIF

C  Find charge (assuming track direction is away from origin)
C  Get rotation angle to bring centre of track circle on x-axis
C
                  PHIC = ATAN2(YC,XC)
                  CP = COS(PHIC)
                  SP = SIN(PHIC)
C
C  All 3 points must be on the same semicircle
C
                  Y1R= Y1*CP - X1*SP
                  Y2R= Y2*CP - X2*SP
                  Y3R= Y3*CP - X3*SP
                  IF (Y2R.GE.0.) THEN
                     IF (Y1R.LT.0. .OR. Y3R.LT.0.) THEN
C GAB:
C                        IPTR1=ITABL(INDBIN,IPTR1,JNDBFW)
                        IPTR1=MY_JNDBFW(IPTR1)
                        GO TO 112
                     ENDIF
                     QQ = 1.
                  ELSE
                     IF (Y1R.GT.0. .OR. Y3R.GT.0.) THEN
C GAB:
C                        IPTR1=ITABL(INDBIN,IPTR1,JNDBFW)
                        IPTR1=MY_JNDBFW(IPTR1)
                        GO TO 112
                     ENDIF
                     QQ = -1.
                  ENDIF
                  PHI0 = PHIC + QQ*PIBY2
                  IF (PHI0.LT.0.) PHI0 = PHI0 + TWOPI

C
C  Now we have a link which passes the tests.
C  Store the new link.  Note that we no longer try to resolve any
C  ambiguities here resulting from more than one link sharing a
C  single hit.  These will be resolved only after trying to extrapolate
C  the links further.
C

                  NLK = NLK + 1
c                  NDONE = IW(INDLNK+1)*IW(INDLNK+2) + 2
c                  IF (NDONE + IW(INDLNK+1).GT.IW(INDLNK)) THEN
c                     NDATA = IW(INDLNK) +  NXLNK
c                     CALL WBANK(IW,INDLNK,NDATA,*800)
c                  ENDIF

c GAB:
c                  KLNK = INDLNK + NDONE
c                  IW(KLNK + 1) = IR3
c                  IW(KLNK + 2) = IR2
c                  IW(KLNK + 3) = IR1
c                  IW(KLNK + 4) = IPTR3
c                  IW(KLNK + 5) = IPTR2
c                  IW(KLNK + 6) = IPTR1
c                  RW(KLNK + 7) = Q*RHO
c                  RW(KLNK + 8) = TANL
c                  RW(KLNK + 9) = PHI0
c                  RW(KLNK +10) = D0
c                  RW(KLNK +11) = Z0
c                  RW(KLNK +12) = S13
c                  RW(KLNK +13) = CHI2Z


C GAB:
                  NLNK = NLNK + 1
C(KH)-------begin
c FIXME: SJA: drop event for to many link, 
c a more ellegant solution is needed   
                  
                  IF (NLNK.GT.NMAXLNK) THEN
                     WRITE(6,*) 'TRSLNK: NMAXLNK TOO SMALL, NLNK =',NLNK 
                     WRITE(6,*) 'THIS IS FATAL UNTIL A BETTER' 
                     WRITE(6,*) 'METHOD IS FOUND FOR CURLERS'
                     WRITE(6,*) 'EVENT DROPPED'
                     CALL TKVETO(1) 
                     IER = 1
                     GO TO 999
                     NLNK=NLNK-1
                  END IF
C(KH)-------end

                  RLNK(NLNK , 1) = REAL(IR3)
                  RLNK(NLNK , 2) = REAL(IR2)
                  RLNK(NLNK , 3) = REAL(IR1)
                  RLNK(NLNK , 4) = REAL(IPTR3)
                  RLNK(NLNK , 5) = REAL(IPTR2)
                  RLNK(NLNK , 6) = REAL(IPTR1)
                  RLNK(NLNK , 7) = QQ*RHO
                  RLNK(NLNK , 8) = TANL
                  RLNK(NLNK , 9) = PHI0
                  RLNK(NLNK ,10) = D0
                  RLNK(NLNK ,11) = Z0
                  RLNK(NLNK ,12) = S13
                  RLNK(NLNK ,13) = CHI2Z
C
C   Compute sagitta errors for the three points. We take account
C   of the uncertainty in the radial position which is significant
C   for tracks at oblique incidence (large beta). The factor 0.1875
C   corresponds to (half pad length)**2/12.
C
                  BETA = ATAN2(-(X3-XC), Y3-YC )
                  IF (BETA.LT.0.) BETA=BETA+TWOPI
                  BETA = TRPHDI(PHI3,BETA)
                  COSB = COS(BETA)
                  SINB = SIN(BETA)
C GAB:
C                  DW2 = RTABL(KTPCO,IC3,JTPCSR)
C                  RW(KLNK +14) = DW2*COSB*COSB + 0.1875*SINB*SINB
C                  DW2 = TPCHITS(6,IC3)**2
                  DW2 = TPCO_SR(IC3)
                  RLNK(NLNK ,14) = DW2*COSB*COSB + PADFAC*SINB*SINB

                  BETA = ATAN2(-(X2-XC), Y2-YC )
                  IF (BETA.LT.0.) BETA=BETA+TWOPI
                  BETA = TRPHDI(PHI2,BETA)
                  COSB = COS(BETA)
                  SINB = SIN(BETA)
C GAB:
C                  DW2 = RTABL(KTPCO,IC2,JTPCSR)
C                  RW(KLNK +15) = DW2*COSB*COSB + 0.1875*SINB*SINB
C                  DW2 = TPCHITS(6,IC2)**2
                  DW2 = TPCO_SR(IC2)
                  RLNK(NLNK ,15) = DW2*COSB*COSB + PADFAC*SINB*SINB
                  BETA = ATAN2(-(X1-XC), Y1-YC )
                  IF (BETA.LT.0.) BETA=BETA+TWOPI
                  BETA = TRPHDI(PHI1,BETA)
                  COSB = COS(BETA)
                  SINB = SIN(BETA)
C GAB:
C                  DW2 = RTABL(KTPCO,IC1,JTPCSR)
C                  RW(KLNK +16) = DW2*COSB*COSB + 0.1875*SINB*SINB
C                  IW(INDLNK + 2) = IW(INDLNK + 2) + 1
C                  DW2 = TPCHITS(6,IC1)**2
                  DW2 = TPCO_SR(IC1)
                  RLNK(NLNK ,16) = DW2*COSB*COSB + PADFAC*SINB*SINB
C
C  Get pointer to the next coordinate in the list for this bin
C
C GAB:
C                  IPTR1=ITABL(INDBIN,IPTR1,JNDBFW)
                  IPTR1=MY_JNDBFW(IPTR1)
                  GO TO 112
  113          CONTINUE
   10       CONTINUE

C
C  Get pointer to the next coordinate in the list for this bin
C
C GAB:
C            IPTR2=ITABL(INDBIN,IPTR2,JNDBFW)
            IPTR2=MY_JNDBFW(IPTR2)
            GO TO 212
  213    CONTINUE
   20 CONTINUE

      GO TO 999


C
c  800 CONTINUE
cc      CALL RERROR('TRSLNK',1,'Error extending link work bank')
c      IER=1

  999 CONTINUE

      RETURN
      END

c------------------------------------------------------------------------------
