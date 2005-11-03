*  
* $Id: trseek.F,v 1.2 2005-08-08 07:09:13 aplin Exp $
*  
* $Log: not supported by cvs2svn $
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
      SUBROUTINE TRSEEK
      IMPLICIT NONE
C****************************************************************************
c* Version number: 0.1      Date of last modification   6.11.99             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
c* Variable inner row for use in unmerging routine                          *
C* Main modifiers from ALEPH to BRAHMS code:                                *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
C
C--------------------------------------------------------------------*
C! Propagate current 3-point link to form chain
C!
C!   Author:     M. Mermikides    06-10-86
C!   Modified:   R. Johnson       19-01-88
C!   Modified:   M.Mermikides     17-08-89
C!
C!   Called by TRKFND
C!
C!    Description
C!    ===========
C! Extract hits from adjacent padrows extrapolating from current
C! starting link.  Given a link joining hits on padrows j, j-1,
C! j-2, we form new  3-point links using the previously assigned
C! hits on rows j-1, j-2 and each hit on row j-3.  A link is
C! acceptable if it passes the following tests:
C!
C!       1.Difference in phi between hits in successive padrows
C!                                   < DPNEXT
C!       2.Difference in z between hits in successive padrows
C!                                   < DZNEXT
C!       3.Radius of curvature of link > RHOMIN
C!       4.Difference between predicted z value of middle point
C!         of new link and old z value < TDZTOL (take account
C!                                             of dip angle)
C!       5.D0 of new link < TD0TOL
C!       6.Z0 of new link < TZ0TOL
C!       7.Consistency in curvature with previous link < CURTOL
C!
C!    The process is repeated until one of the following occurs :
C!             1. Innermost padrow reached
C!             2. More than one link passes above tests
C!             3. No link is found over MAXGAP successive padrows
C!
C -----------------------------------------------------------------*
C
#include "include/fndparms.inc"
#include "include/alcons.inc"
#include "include/tpgeom.inc"
c#include "marlin_tpcgeom.h"
#include "include/tchain.inc"
#include "include/trfdct.inc"
#include "include/tcbins.inc"
#include "include/tpco.inc"
#include "include/tpcdims.inc"
#include "include/my_jndb.inc"

      INTEGER IROW,IR2,IR3,IPTR,IPTR2,IPTR3,IC,IC2,IC3
      INTEGER IBIN,NGAP,NFOUN
      REAL PHI1,PHI2,PHI3,Z1,Z2,Z3,X1,X2,X3,DZ2,Y1,Y2,Y3,XC,YC,RC,DET
      REAL A1,B1,C1,A2,B2,C2,CURVL,COUNT,SMCUR,DPHI,RHOSQ,RHO
      REAL C12SQ,S13,TANL,Z2PR,ZDIFF,X01,X02,Y01,Y02,D1SQ,D2SQ
      REAL CH10,CH12,CH13,ARG,S12,CH2Z,D0,X0,Y0,Z0,PSI0,CP,SP,PHIC
      REAL Y1R,Y2R,Y3R,QQ,PHI0,ALPHA,BETA,SINB,COSB,DSAG2,DCUR2,CDIFF,CH
      REAL QORHO

      REAL PADFAC
      REAL TRPHDI
      EXTERNAL TRPHDI
c GAB:
#include "include/unmrg.inc"

C----------------------------------------------------------------------
cC
      LOGICAL FIRST
      SAVE PADFAC
      DATA FIRST/.TRUE./

C
c#include "bkgmode.inc"
c#include "gcbank.inc"
#include "include/tpcfunc2.inc"


#include "include/tpcfunc.inc"
C----------------------------------------------------------------------

C
c GAB:
      IF (FIRST) THEN
        FIRST=.FALSE.
        PADFAC=TPCPADR**2/12.
c        NTPCO=NAMIND('TPCO')
c        NTCRL=NAMIND('TCRL')
      ENDIF
C
c      KTPCO = IW(NTPCO)
c      KTCRL = IW(NTCRL)
C

      IF(LUNMRG)THEN
       MINROW=MINMRG
      ELSE
       MINROW=1
      ENDIF

C Starting row (outer) is the middle point of old link
C
      IR3 = IRCHAI(2)
      IR2 = IRCHAI(3)
      IPTR3 = IHCHAI(2)
      IPTR2 = IHCHAI(3)
C
      IROW = IR2 - 1

c GAB:
c      IC3 = ITABL(INDBIN,IPTR3,JNDBCO)
      IC3 = MY_JNDBCO(IPTR3)
C GAB:
C      J3 =  KROW(KTPCO,IC3)
C      IC2 = ITABL(INDBIN,IPTR2,JNDBCO)
C      J2 =  KROW(KTPCO,IC2)
      IC2 = MY_JNDBCO(IPTR2)

C GAB:
C      PHI3 = RW(J3 + JTPCPH)
C      PHI2 = RW(J2 + JTPCPH)
C      Z3   = RW(J3 + JTPCZV)
C      Z2   = RW(J2 + JTPCZV)
C      DZ2 = RW(J2 + JTPCSZ)
C      JX3 = KROW(INDCXY,IC3)
C      JX2 = KROW(INDCXY,IC2)
C      X3   = RW(JX3 + 1)
C      Y3   = RW(JX3 + 2)
C      X2   = RW(JX2 + 1)
C      Y2   = RW(JX2 + 2)
C      A2   = X3-X2
C      B2   = Y3-Y2
C      C2   = 0.5*(X2*X2-X3*X3+Y2*Y2-Y3*Y3)

C      PHI3 = TPCPHI(IC3)
C      PHI2 = TPCPHI(IC2)
C      Z3   = TPCHITS(3,IC3)
C      Z2   = TPCHITS(3,IC2)
C      DZ2 =  TPCHITS(7,IC2)**2
C      X3   = TPCHITS(1,IC3)
C      Y3   = TPCHITS(2,IC3)
C      X2   = TPCHITS(1,IC2)
C      Y2   = TPCHITS(2,IC2)

      PHI3 = TPCO_PH(IC3)
      PHI2 = TPCO_PH(IC2)
      Z3   = TPCO_ZV(IC3)
      Z2   = TPCO_ZV(IC2)
      DZ2 =  TPCO_SZ(IC2)
      X3   = TPCO_XV(IC3)
      Y3   = TPCO_YV(IC3)
      X2   = TPCO_XV(IC2)
      Y2   = TPCO_YV(IC2)
      A2   = X3-X2
      B2   = Y3-Y2
      C2   = 0.5*(X2*X2-X3*X3+Y2*Y2-Y3*Y3)

C
C  Get initial link parameters
C
      IF (ABS(PCHAIN(1)) .LE. 0.0001) GOTO 999
      CURVL = 1./PCHAIN(1)
      COUNT = 1.
      SMCUR = CURVL
      NGAP = 0


C
 100  NFOUN = 0
C
C  Find the appropriate phi bin and loop over all hits within it
C
      IBIN= IXWRAP(INT(TCOFFS+PHI2/TC2DEL))
C(KH) IPTR= IBNPTR(IBIN,IROW)
C(KH)-----begin
      IF (IROW.GT.0) THEN
        IPTR= IBNPTR(IBIN,IROW)
      ELSE
        IPTR=0
      END IF
C(KH)-----end

  440 IF (IPTR.EQ.0) GO TO 441

C GAB:
C         IC= ITABL(INDBIN,IPTR,JNDBCO)
C         J1=  KROW(KTPCO,IC)
C         JX1= KROW(INDCXY,IC)
C         PHI1 = RW(J1 + JTPCPH)
C         DPHI = TRPHDI(PHI2,PHI1)
         IC= MY_JNDBCO(IPTR)
C         PHI1 =TPCPHI(IC)
         PHI1 =TPCO_PH(IC)
         DPHI = TRPHDI(PHI2,PHI1)

C
C  Check on difference in phi and z
C
         IF (ABS(DPHI).GT.DPNEXT) THEN
C GAB:
C            IPTR= ITABL(INDBIN,IPTR,JNDBFW)
            IPTR= MY_JNDBFW(IPTR)
            GO TO 440
         ENDIF

C GAB:
C         Z1 = RW(J1 + JTPCZV)
         Z1 = TPCO_ZV(IC)

         IF (ABS(Z2-Z1).GT.DZNEXT) THEN
C GAB:
C            IPTR= ITABL(INDBIN,IPTR,JNDBFW)
            IPTR= MY_JNDBFW(IPTR)
            GO TO 440
         ENDIF

C
C GAB:
C         X1=RW(JX1 + 1)
C         Y1=RW(JX1 + 2)
         X1=TPCO_XV(IC)
         Y1=TPCO_YV(IC)
         A1=X2-X1
         B1=Y2-Y1
         C1=0.5*(X1*X1-X2*X2+Y1*Y1-Y2*Y2)
         DET= A1*B2-A2*B1
         IF (ABS(DET) .LE. .00001) THEN
C GAB:
C            IPTR= ITABL(INDBIN,IPTR,JNDBFW)
            IPTR= MY_JNDBFW(IPTR)
            GO TO 440
         ENDIF

C
C Centre and radius of curvature (rho)
C
         XC= (B1*C2-B2*C1)/DET
         YC=-(A1*C2-A2*C1)/DET
         RHOSQ=(X1-XC)**2+(Y1-YC)**2
         RHO = SQRT(RHOSQ)
C
C  Check on radius of curvature
C
         IF (RHO.LT.RHOMIN) THEN
C GAB:
C            IPTR= ITABL(INDBIN,IPTR,JNDBFW)
            IPTR= MY_JNDBFW(IPTR)
            GO TO 440
         ENDIF
C
C  Get turning angle to middle and inner padrows
C
         C12SQ = A1*A1+B1*B1
         CH12 = SQRT(C12SQ)
         ARG = CH12/(2.*RHO)
         IF (ARG.GE.1.0) THEN
C GAB:
C            IPTR= ITABL(INDBIN,IPTR,JNDBFW)
            IPTR= MY_JNDBFW(IPTR)
            GO TO 440
         ENDIF
C
C  Approximate arc length in the case rho > 20 cm and chord < 15 cm
C
         IF (CH12.LT.15. .AND. RHO. GT.20.) THEN
            S12 = CH12*(1. + ARG*ARG/8.)
         ELSE
            S12 = 2.*RHO*ASIN(ARG)
         ENDIF
C
         CH13=SQRT((X1-X3)**2+(Y1-Y3)**2)
         ARG = CH13/(2.*RHO)

         IF (ARG.GE.1.0) THEN
C GAB:
C            IPTR= ITABL(INDBIN,IPTR,JNDBFW)
            IPTR= MY_JNDBFW(IPTR)
            GO TO 440
         ENDIF

         IF (CH13.LT.15. .AND. RHO. GT.20.) THEN
            S13 = CH13*(1. + ARG*ARG/8.)
         ELSE
            S13 = 2.*RHO*ASIN(ARG)
         ENDIF
C
C  Predict z. Test on distance of Z2 from track orbit
C
         TANL=(Z3-Z1)/S13
         Z2PR=Z1 + S12*TANL
         ZDIFF = Z2PR-Z2
         CH2Z = ZDIFF*ZDIFF/DZ2

         IF (CH2Z.GT.TDZTOL) THEN
C GAB:
C            IPTR= ITABL(INDBIN,IPTR,JNDBFW)
            IPTR= MY_JNDBFW(IPTR)
            GO TO 440
         ENDIF

C
         RC=SQRT(XC*XC+YC*YC)
         X01=XC*(1.+RHO/RC)
         X02=XC*(1.-RHO/RC)
C
C Protect against circle centre on Y-axis
C
         IF (ABS(XC).LT.0.00001) THEN
            Y01 = YC*(1.+ RHO/RC)
            Y02 = YC*(1.- RHO/RC)
         ELSE
            Y01=YC*X01/XC
            Y02=YC*X02/XC
         ENDIF
C
C Find closest point to origin (D0)
C D0 is positive if track encloses origin, otherwise negative
C
         D1SQ=X01*X01+Y01*Y01
         D2SQ=X02*X02+Y02*Y02
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
C  D0 test
C
         IF (D0.GT.TD0TOL) THEN
C GAB:
C            IPTR= ITABL(INDBIN,IPTR,JNDBFW)
            IPTR= MY_JNDBFW(IPTR)
            GO TO 440
         ENDIF

         IF (RC.GT.RHO)   D0=-D0
         CH10=SQRT((X1-X0)**2+(Y1-Y0)**2)
         PSI0=2.*ASIN(CH10/(2.*RHO))
C
C  Z0 test
C
         Z0=Z1-RHO*PSI0*TANL

         IF (ABS(Z0).GT.TZ0TOL) THEN
C GAB:
C            IPTR= ITABL(INDBIN,IPTR,JNDBFW)
            IPTR= MY_JNDBFW(IPTR)
            GO TO 440
         ENDIF
C
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
C               IPTR= ITABL(INDBIN,IPTR,JNDBFW)
               IPTR= MY_JNDBFW(IPTR)
               GO TO 440
            ENDIF
            QQ = 1.
         ELSE
            IF (Y1R.GT.0. .OR. Y3R.GT.0.) THEN
C GAB:
C               IPTR= ITABL(INDBIN,IPTR,JNDBFW)
               IPTR= MY_JNDBFW(IPTR)
               GO TO 440
            ENDIF
            QQ = -1.
         ENDIF

         PHI0 = PHIC + QQ*PIBY2
         IF (PHI0.LT.0.) PHI0 = PHI0 + TWOPI
C
C  Check on curvature (must be consistent with last link)
C  d(curv) = 8.*d(sagitta)/chord**2.  Use r-phi error of last
C  point for d(sagitta)
C
C  We include a radial component of sagitta error to account for the
C  uncertainty in the radial hit coordinate for tracks
C (such as curlers and off-centre cosmics) which cross the pads
C  transversly.  The factor 0.1875 corresponds to the variance
C  for a uniform distribution over the interval of half a pad
C  length ( 1.5**2/12)
C GAB: Replace 0.1875 factor by PADFAC = rpad**2/12
C
         ALPHA = ATAN2((Y3-Y1),(X3-X1))
         IF (ALPHA.LT.0.0) ALPHA = ALPHA + TWOPI
         BETA = TRPHDI(PHI2,ALPHA)
         COSB = COS(BETA)
         SINB = SIN(BETA)
c GAB:
c         DSAG2 = RW(J1+JTPCSR)*COSB*COSB + 0.1875*SINB*SINB
C         DSAG2 = (TPCHITS(6,IC)**2)*COSB*COSB + PADFAC*SINB*SINB
         DSAG2 = TPCO_SR(IC)*COSB*COSB + PADFAC*SINB*SINB
         DCUR2 = 64.*DSAG2/(C12SQ*C12SQ)
         CDIFF = QQ/RHO - CURVL
         CH = CDIFF*CDIFF/DCUR2

         IF (CH.GT.CURTOL) THEN
C GAB:
C            IPTR= ITABL(INDBIN,IPTR,JNDBFW)
            IPTR= MY_JNDBFW(IPTR)
            GO TO 440
         ENDIF
C
C Save the curvature in order to make a running average if this hit
C is indeed added to the chain
C
         QORHO=QQ/RHO
C
C No of consistent hits found
C
         NFOUN =  NFOUN + 1

         IF (NFOUN.EQ.1) THEN
            NHCHAI = NHCHAI + 1
            IRCHAI(NHCHAI) = IROW
            IHCHAI(NHCHAI) = IPTR
            D2CHAI(NHCHAI) = DSAG2
C  If more than one hit is found, delete last one in chain and exit
         ELSE
            NHCHAI = NHCHAI - 1
            GO TO 500
         ENDIF

C
C  Get pointer to the next coordinate in this bin
C
C GAB:
C         IPTR= ITABL(INDBIN,IPTR,JNDBFW)
         IPTR= MY_JNDBFW(IPTR)
         GO TO 440
  441 CONTINUE

C New root.
C Exit when chain is broken by more than MAXGAP gaps
C
   50 IF (NFOUN.EQ.0) THEN
         NGAP = NGAP + 1

         IF (NGAP.GT.MAXGAP)then
           GO TO 500
         endif
         IROW = IROW - 1
C GAB:
C         IF (IROW.EQ.0) GO TO 500
         IF (IROW.EQ.MINROW) GO TO 500
C
C No hit in this row.  Try the next row inward.
C
         GO TO 100
      ELSE
C
C Get running mean of curvature
C
         SMCUR = SMCUR + QORHO
         COUNT = COUNT + 1.
         CURVL = SMCUR/COUNT
C
         X3 = X2
         Y3 = Y2
         Z3 = Z2
         IR2  = IROW
         IPTR2  = IHCHAI(NHCHAI)

C GAB:
C         IC2 = ITABL(INDBIN,IPTR2,JNDBCO)
         IC2 = MY_JNDBCO(IPTR2)
C GAB:
C         J2  = KROW(KTPCO,IC2)
C         JX2 = KROW(INDCXY,IC2)
C         PHI2 = RW(J2 + JTPCPH)
C         Z2   = RW(J2 + JTPCZV)
C         X2   = RW(JX2 + 1)
C         Y2   = RW(JX2 + 2)
C         PHI2 = TPCPHI(IC2)
C         Z2   = TPCHITS(3,IC2)
C         X2   = TPCHITS(1,IC2)
C         Y2   = TPCHITS(2,IC2)
         PHI2 = TPCO_PH(IC2)
         Z2   = TPCO_ZV(IC2)
         X2   = TPCO_XV(IC2)
         Y2   = TPCO_YV(IC2)

         A2   = X3-X2
         B2   = Y3-Y2
         C2   = 0.5*(X2*X2-X3*X3+Y2*Y2-Y3*Y3)
         IROW = IROW - 1
         IF (IROW.EQ.0)           GO TO 500

         NGAP = 0
C
C Try to extend chain again inward by one row
C
         GO TO 100
      ENDIF
C
C  We have now followed the chain as far as we can.
C
 500  CONTINUE

 999  RETURN
      END

c-----------------------------------------------------------------------------
