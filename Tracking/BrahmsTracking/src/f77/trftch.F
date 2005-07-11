*  
* $Id: trftch.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
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
      SUBROUTINE TRFTCH
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
C----------------------------------------------------------------
C! Perform "helix" fit to coordinates of current chain
C!
C!    Author:    M. Mermikides
C!    Modified:  M. Mermikides 22-04-87
C!    Modified:  R. Johnson    19-01-88
C!    Modified:  F. Weber       8-6-89
C!    Modified:  M.Mermikides  17-08-89
C!
C!   Called by TRKFND
C!
C!    Description
C!    ===========
C!    Use internal errors on r-phi and z coordinates
C!    Add m/s error and external r-phi error for oblique tracks
C----------------------------------------------------------------

C
#include "include/fndparms.inc"
#include "include/alcons.inc"
#include "include/tpco.inc"
c#include "include/tpgeom.inc"
#include "marlin_tpcgeom.h"
#include "include/tchajj.inc"
#include "include/rflags.inc"
#include "include/tchain.inc"
#include "include/tpccon.inc"
#include "include/tpcirc.inc"
#include "include/tplinc.inc"
#include "include/trfdct.inc"
#include "include/rconds.inc"
#include "include/my_jndb.inc"
#include "include/tpcdbg.inc"

      LOGICAL FIRST
C-----------------------------------------------------------------------
C   These are needed for the Fast Helix Fitting Routine:
C
      INTEGER MPT
      PARAMETER (MPT = LTPDRO+2)
      REAL RCIRCL,PHCIRC,ERRMAT
      DIMENSION RCIRCL(MPT),PHCIRC(MPT),ERRMAT(15)
      DOUBLE PRECISION XFIT(MPT),YFIT(MPT),WFIT(MPT)
      REAL HELPAR
      DIMENSION HELPAR(5)
      REAL SCATC,P2CUT,SSIGX,SSIGZ,SIGRP,SIGZ,Z1,Z2,CHIRP,CHISZ,CHITOT
      REAL ROFP,SGN,RHO,TANL,PHI0,D0,Z0,X0,Y0,PHIC,RR,DPSI,SECL,PP,P2
      REAL DX,DY,SS,SIGMS,SCALX,SCALZ,PULMX,CURV,COSP0,SINP0,TERM,DXY
      REAL PULL,ZZ,DZ

      INTEGER I,J,NCIRCL,NTRY,IR1,IR2,IOPT,KCIR,KS,JMAX
C----------------------------------------------------------------------
C MY STUFF:
      INTEGER IHIT
      SAVE
C----------------------------------------------------------------------
C
C  SCATC  = m/s constant = .0141**2 / (3. * X0), X0 = 10983 cm for Ar/CH
C  P2CUT  = Momentum squared above which we ignore m/s
C
      DATA SCATC/6.034E-9/
      DATA P2CUT/9.0/
      DATA FIRST/.TRUE./

c#include "bkgmode.inc"
c#include "gcbank.inc"
#include "include/tpcfunc2.inc"


#include "include/tpcfunc.inc"
C      ICPTR(M1)= ITABL(INDBIN,M1,JNDBCO)
C
C----------------------------------------------------------------------
C
c      IF (FIRST) THEN
c        FIRST=.FALSE.
c        NTPCO=NAMIND('TPCO')
c        NTCRL=NAMIND('TCRL')
c      ENDIF
C
C GAB:
C      KTPCO = IW(NTPCO)
C      IF (KTPCO.EQ.0) GO TO 999
C      KTCRL = IW(NTCRL)
C      IF (KTCRL.EQ.0) GO TO 999

C
C  Clear dropped hit list
C
      CALL VZERO(MCHAIN,LTPDRO)
C
      NTRY = 0
C
  500 NCIRCL = 0
      NTRY = NTRY + 1
      IF (NTRY.GT.MAXTRY)THEN
       IF (TPCDBG.GT.0) THEN
         WRITE(*,*)'TRFTCH: NTRY exceeds MAXTRY=',MAXTRY
       END IF
       GO TO 999
      ENDIF
      SSIGX = 0.
      SSIGZ = 0.
C
C  Transfer hits to fitter array in increasing R order
C
      DO 10 J=NHCHAI,1,-1
C
         IF (MCHAIN(J).EQ.0) THEN
            NCIRCL = NCIRCL+1
C GAB:
C            XCIRCL(NCIRCL) = RTABL(INDCXY,ICPTR(IHCHAI(J)),1)
C            YCIRCL(NCIRCL) = RTABL(INDCXY,ICPTR(IHCHAI(J)),2)
C            ZCIRCL(NCIRCL) = RTABL(KTPCO,ICPTR(IHCHAI(J)),JTPCZV)
C            RCIRCL(NCIRCL) = RTABL(KTPCO,ICPTR(IHCHAI(J)),JTPCRV)
C            PHCIRC(NCIRCL) = RTABL(KTPCO,ICPTR(IHCHAI(J)),JTPCPH)

            IHIT=MY_JNDBCO(IHCHAI(J))

C            XCIRCL(NCIRCL) = TPCHITS(1,IHIT)
C            YCIRCL(NCIRCL) = TPCHITS(2,IHIT)
C            ZCIRCL(NCIRCL) = TPCHITS(3,IHIT)
C            RCIRCL(NCIRCL) = TPCR(IHIT)
C            PHCIRC(NCIRCL) = TPCPHI(IHIT)
            XCIRCL(NCIRCL) = TPCO_XV(IHIT)
            YCIRCL(NCIRCL) = TPCO_YV(IHIT)
            ZCIRCL(NCIRCL) = TPCO_ZV(IHIT)
            RCIRCL(NCIRCL) = TPCO_RV(IHIT)
            PHCIRC(NCIRCL) = TPCO_PH(IHIT)
            
C     SJAWRT
c            write(90,*) "XCIRCL(NCIRCL) = ", XCIRCL(NCIRCL)  
c            write(90,*) "YCIRCL(NCIRCL) = ", YCIRCL(NCIRCL)            
c            write(90,*) "ZCIRCL(NCIRCL) = ", ZCIRCL(NCIRCL) 
c            write(90,*) "RCIRCL(NCIRCL) = ", RCIRCL(NCIRCL) 
c            write(90,*) "PHCIRC(NCIRCL) = ", PHCIRC(NCIRCL) 

            IF (NCIRCL.EQ.1) THEN
               Z1 = ZCIRCL(NCIRCL)
               IR1 = IRCHAI(J)
            ENDIF
C
            SIGRP = D2CHAI(J)
C GAB:
C            SIGZ =  RTABL(KTPCO,ICPTR(IHCHAI(J)),JTPCSZ)
C            SIGZ =  TPCHITS(7,IHIT)**2
            SIGZ =  TPCO_SZ(IHIT)
            IR2 = IRCHAI(J)
            Z2 = ZCIRCL(NCIRCL)
            WCIRCL(NCIRCL) = 1./SIGRP
            WLINEA(NCIRCL) = 1./SIGZ
            SSIGX = SSIGX + SIGRP
            SSIGZ = SSIGZ + SIGZ
         ENDIF
   10 CONTINUE

C
C  Fit if more than 3 hits
C
      IF (NCIRCL.LE.3)  GO TO 999
C
C  Mean sigma**2
C
      SSIGX = SSIGX/FLOAT(NCIRCL)
      SSIGZ = SSIGZ/FLOAT(NCIRCL)


C
C-----------------------------------------------------------------------
C   Call the modified Fast Helix Fitting Routine                       !
C                                                                      !
C  least sqrs fit w/error mat calc.
   25 IOPT = 2
      DO 26 I = 1,NCIRCL
          XFIT(I) = XCIRCL(I)
          YFIT(I) = YCIRCL(I)
          WFIT(I) = WCIRCL(I)
   26 CONTINUE

      CALL TFITHL(NCIRCL,XFIT,YFIT,RCIRCL,PHCIRC,WFIT,
     1            ZCIRCL,WLINEA,IOPT,
     2            HELPAR,ERRMAT,CHIRP,CHISZ)

C
C  Convert the CHI2 values to the ones this routine expects, calculate C
C  as a check vs. TFITHL failure
C
      CHITOT = CHIRP + CHISZ
      CHIRP=CHIRP/FLOAT(NCIRCL-3)
      CHISZ=CHISZ/FLOAT(NCIRCL-2)
C
      IF (CHITOT.GT.5000.0) THEN
C
C   TFITHL failed.  Borrowing code from L. Garrido's UFITMS
C   routine, I replace any zero diagonal elements of ERRMAT with the
C   average errors on a 3 GeV particle.  This assumes that the other
C   elements of the inverted ERRMAT are o.k.:
C
C     conversion radius of track <=> momentum
C     radius in meter , B in Tesla, p in GeV/c  q in units of e
C
C      p = 0.29979 * q * B * r
C
C     R[cm] = ROFP * P[Gev/c]:
C
         ROFP=1./(0.29979*FIELRC/10.)*100.
         IF (ERRMAT(1).EQ.0.0) ERRMAT(1)=1.0/(.005/ROFP)**2
         IF (ERRMAT(3).EQ.0.0) ERRMAT(3)=1.0/.002**2
         IF (ERRMAT(6).EQ.0.0) ERRMAT(6)=1.0/.003**2
         IF (ERRMAT(10).EQ.0.0) ERRMAT(10)=1.0/.007**2
         IF (ERRMAT(15).EQ.0.0) ERRMAT(15)=1.0/.007**2

      END IF
C
C   HELPAR contains the ALEPH-standard helix parameters; must be loaded
C   into internal variables
C
      SGN = SIGN(1.,HELPAR(1))
      RHO = ABS(1./HELPAR(1))
      TANL = HELPAR(2)
      PHI0 = HELPAR(3)
      D0 = HELPAR(4)
      Z0 = HELPAR(5)
      X0 = ABS(D0)*COS(PHI0)
      Y0 = ABS(D0)*SIN(PHI0)
C
C   Now compute SCIRCL values
C
      PHIC = PHI0 - PIBY2
      IF(SGN*D0.GT.0.0)   PHIC = PHI0 + PIBY2
      X0 = ABS(D0)*COS(PHIC)
      Y0 = ABS(D0)*SIN(PHIC)

C
C Find S-values measured from the point (X0,Y0); intercept
C in s-z plane then gives Z0. (Assume hits are ordered in S)
C
      SCIRCL(1)=0.
      KCIR=NCIRCL
      RR = RHO
      DPSI=SQRT((XCIRCL(1)-X0)**2
     1            +(YCIRCL(1)-Y0)**2)/(2.*RR)
      SCIRCL(1)= 2.*RR*ASIN(DPSI)

      DO 137 KS=2,KCIR
         DPSI=SQRT((XCIRCL(KS)-XCIRCL(KS-1))**2
     1            +(YCIRCL(KS)-YCIRCL(KS-1))**2)/(2.*RR)
C
C Check if we turn through pi
C

         IF (DPSI.GT.1.0) THEN
            NCIRCL=KS
            GO TO 138
         ENDIF
         SCIRCL(KS)=SCIRCL(KS-1) + 2.*RR*ASIN(DPSI)
  137 CONTINUE
  138 CONTINUE
C                                                                      !
C-----------------------------------------------------------------------
C
C  Take m/s into account in chisq. For now we just rescale
C  by average error.
C
      SECL = SQRT(1. + TANL*TANL)
      PP = RHO*BCFGEV*SECL
      P2 = PP*PP
C
C Skip track with p > P2CUT
C
      IF (P2.GT.0.0 .AND.P2.LT.P2CUT) THEN
         DX = XCIRCL(1) - XCIRCL(NCIRCL)
         DY = YCIRCL(2) - YCIRCL(NCIRCL)
C
C Approx. arc by chord
C
         SS = SQRT(DX*DX + DY*DY)*SECL/4.
         SIGMS = SCATC*SS*SS*SS/P2
         SCALX = SSIGX/(SIGMS +SSIGX)
         SCALZ = SSIGZ/(SIGMS +SSIGZ)
      ELSE
         SCALX = 1.
         SCALZ = 1.
      ENDIF

C
C  In case of failed circle fit, reject hit with worst pull
C
      IF (CHIRP*SCALX.GT.CHLIM1) THEN
         PULMX = 0.
         CURV = - SGN/(RHO-D0)
         COSP0 = COS(PHI0)
         SINP0 = SIN(PHI0)
C
C  Compute residuals
C
         DO 20 J = 1,NCIRCL
            TERM = XCIRCL(J)*XCIRCL(J) + YCIRCL(J)*YCIRCL(J) - D0*D0
            DXY = 0.5*CURV*TERM + XCIRCL(J)*SINP0
     &           -YCIRCL(J)*COSP0 + SGN*D0
            PULL = ABS(DXY)*SQRT(WCIRCL(J))
            IF (PULL.GT.PULMX) THEN
               PULMX = PULL
               JMAX = J
            ENDIF
C            PRINT 556, XCIRCL(J),YCIRCL(J),ZCIRCL(J),DXY,PULL
C 556         FORMAT ( 2(5X,3F10.3,2X,G10.3,2F7.2) )
   20    CONTINUE
         MCHAIN(NHCHAI + 1 - JMAX) = 1
         GO TO 500
C
      ENDIF

C
C  Do the same in case of bad s-z fit
C
      IF (CHISZ*SCALZ.GT.CHLIM2) THEN
         PULMX = 0.
         DO 22 J = 1,NCIRCL
            ZZ = Z0 + SCIRCL(J)*TANL
            DZ = ZCIRCL(J) - ZZ
            PULL = ABS(DZ)*SQRT(WLINEA(J))
            IF (PULL.GT.PULMX) THEN
               PULMX = PULL
               JMAX = J
            ENDIF
   22    CONTINUE
         MCHAIN(NHCHAI + 1 - JMAX) = 1
         GO TO 500
      ENDIF

C
C Update chain parameters
C
      PCHAIN(1) = RHO*SGN
      PCHAIN(JTCHTL) = TANL
      PCHAIN(JTCHP0) = PHI0
      PCHAIN(JTCHD0) = D0
      PCHAIN(JTCHZ0) = Z0
      PCHAIN(JTCHLC) = SCIRCL(NCIRCL) - SCIRCL(1)
      PCHAIN(JTCHC1) = CHIRP
      PCHAIN(JTCHC2) = CHISZ
      ISTCHA = 10000*IR1 + 100*IR2 + NCIRCL
      CALL UCOPY (ERRMAT(1),PCHAIN(JTCHER),15)

C     SJAWRT
c      write(*,*) "RHO*SGN = ",PCHAIN(1)
c      write(*,*) "TANL   =  ",PCHAIN(JTCHTL)
c      write(*,*) "PHI0   =  ",PCHAIN(JTCHP0)
c      write(*,*) "D0     =  ",PCHAIN(JTCHD0)
c      write(*,*) "Z0     =  ",PCHAIN(JTCHZ0)
c      write(*,*) "SCIRCL =  ",PCHAIN(JTCHLC)
c      write(*,*) "CHIRP  =  ",PCHAIN(JTCHC1)
c      write(*,*) "CHISZ  =  ",PCHAIN(JTCHC2)
c
c      DO I = 1, 15
c         write(*,*) "ERRMAT(", i ,") = ",ERRMAT(i)
c      enddo
c      

C
  999 RETURN
      END
C------------------------------------------------------------------
