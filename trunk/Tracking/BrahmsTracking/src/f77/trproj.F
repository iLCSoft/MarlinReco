*  
* $Id: trproj.F,v 1.2 2005-08-08 07:09:13 aplin Exp $
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
      SUBROUTINE TRPROJ
      IMPLICIT NONE
C****************************************************************************
c* Version number: 0.1      Date of last modification   6.11.99             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
c* Variable inner row, for unmerging routine                                *
C* Main modifiers from ALEPH to BRAHMS code:                                *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
C
C---------------------------------------------------------------*
C! Projects chain to find more hits                             *
C!                                                              *
C!   Author:     M. Mermikides    06-10-86                      *
C!   Modified:   R. Johnson       19-01-88                      *
C!   Modified:   M. Mermikides    19-08-89
C!                                                              *
C!   Called by TRKFND                                           *
C!                                                              *
C!  Road width is Min(TDZSTR,TDZMAX)*Min(TDWSTR,TDWMAX)         *
C!  The chain is extended iteratively until there are no more   *
C!  hits that can be added.                                     *
C---------------------------------------------------------------*

C
#include "include/fndparms.inc"
#include "include/alcons.inc"
#include "include/tchajj.inc"
#include "include/tpgeom.inc"
c#include "marlin_tpcgeom.h"
#include "include/tchain.inc"
#include "include/tpgeop.inc"
#include "include/trfdct.inc"
#include "include/tcbins.inc"
#include "include/tpco.inc"
#include "include/tpcdims.inc"
#include "include/my_tcrl.inc"
#include "include/my_jndb.inc"

C
C  JROAD holds hits for each semicircle
C
      INTEGER JROAD,IRSAV,IHSAV,MHSAV
      REAL D2SAV,PSAVE,SROAD,V
      DIMENSION JROAD(LTPDRO+2),IRSAV(LTPDRO+2),IHSAV(LTPDRO+2)
      DIMENSION PSAVE(LTPDRO+4)
      DIMENSION D2SAV(LTPDRO+2),MHSAV(LTPDRO+2),SROAD(LTPDRO+2)
      DIMENSION V(4,2)
C----------------------------------------------------------------------
      INTEGER JSPAN,IR1,IR2,NSAVE,ISTSV,M,K,NXTR1,NCOR,IERR,NFND1
      INTEGER IBIN,IPTR,IC,ISAV1,NEW
      REAL RHO,QQ,TANL,PHI0,D0,DMAX,DMIN,DELZ,PHIC,RC,XC,YC,RAD
      REAL ZP1,PHI1,BETA,SINB,COSB,DELW,PHI,Z,DPHI,DW,DZ,SAGER2

C(KH) ---begin
      INTEGER NLOOP,NTRY
C(KH) ---end

      REAL TRPHDI
      EXTERNAL TRPHDI

      REAL PADFAC
      SAVE PADFAC
C----------------------------------------------------------------------
      LOGICAL TURN,FIRST
      DATA FIRST/.TRUE./
C----------------------------------------------------------------------
c GAB:
      REAL CHI2,CHI2MIN
      INTEGER NSKPMX,NSKIP,MMAX,MSTART,CHNSIZ,ISIZ
#include "include/unmrg.inc"
C----------------------------------------------------------------------
c#include "bkgmode.inc"
c#include "gcbank.inc"
#include "include/tpcfunc2.inc"


      SAVE
#include "include/tpcfunc.inc"

C
C GAB:
      IF (FIRST) THEN
        PADFAC=TPCPADR**2/12.
        FIRST=.FALSE.
C GAB: experimental code:
        NSKPMX=4

c        NTPCO=NAMIND('TPCO')
c        NTCRL=NAMIND('TCRL')
      ENDIF
C
      IF (ISTCHA.EQ.0) GO TO 999
C(KH) --- begin
      NTRY=1
C(KH) --- end
C
C   Skip if we already have full track
C
C GAB:
C  200 IF(NHCHAI.EQ.NRTPC) GO TO 999
      IF(LUNMRG)THEN
       MINROW=MINMRG
      ELSE
       MINROW=0
      ENDIF
200   IF(NHCHAI.GE.NRTPC-MINROW) GO TO 999

C
c GAB:
c      KTPCO = IW(NTPCO)
c      IF (KTPCO.EQ.0) GO TO 999
c      KTCRL = IW(NTCRL)
c      IF (KTCRL.EQ.0) GO TO 999
C

      RHO = ABS(PCHAIN(1))
      QQ = SIGN(1.,PCHAIN(1))
      TANL = PCHAIN(JTCHTL)
      PHI0 = PCHAIN(JTCHP0)
      D0 = PCHAIN(JTCHD0)
      DMAX = 2.*RHO - D0
      TURN = .FALSE.
      IF(DMAX.LE.TPCACRO) TURN=.TRUE.
C
C  Get first and last padrow to be looked at.  For straight-through track
C  we consider range covered by extending present chain length in both
C  directions
C
      IF(.NOT.TURN) THEN
         JSPAN = IRCHAI(1)-IRCHAI(NHCHAI)+1
         IF (JSPAN.GE.8) THEN
C GAB:
           IF(LUNMRG)THEN
            IR1 = MINMRG
           ELSE
            IR1=1
           ENDIF
           IR2 = NRTPC
         ELSE
c gab:
           IF(LUNMRG)THEN
            IR1 = MAX(IRCHAI(NHCHAI)-JSPAN, MINMRG)
           ELSE
            IR1 = MAX(IRCHAI(NHCHAI)-JSPAN, 1)
           ENDIF
           IR2 = MIN(IRCHAI(1) + JSPAN, NRTPC)
         ENDIF

C
C  Case of spiral. Go up to furthest point of approach.
C
      ELSE
         DMIN = ABS(D0)
         IR1 = INT((DMIN-TPPROW(1))/TPDRST(1)) + 1
c gab:
         IF(LUNMRG)THEN
          IR1 = MAX(IR1,MINMRG)
         ELSE
          IR1 = MAX(IR1,1)
         ENDIF

         IR2 = INT((DMAX-TPPROW(1))/TPDRST(1)) + 1
         IR2 = MIN(IR2,NRTPC)
C
C  Return if we have all hits
C
         IF (NHCHAI.GT.IR2-IR1+1) GO TO 999
      ENDIF
C
C  Save initial chain in case fit fails
C
      NSAVE = NHCHAI
      CALL UCOPY(IRCHAI(1),IRSAV,NHCHAI)
      CALL UCOPY(IHCHAI(1),IHSAV,NHCHAI)
      CALL UCOPY(MCHAIN(1),MHSAV,NHCHAI)
      CALL UCOPY(D2CHAI(1),D2SAV,NHCHAI)
      CALL UCOPY(PCHAIN(1),PSAVE,NPARCH)
      ISTSV = ISTCHA
C
C  Get z road width allowing for track inclination
C

      DELZ = TDZSTR*SQRT(1. + TANL*TANL)
      IF(DELZ.GT.TDZMAX)    DELZ=TDZMAX
      PHIC = PHI0 - QQ*PIBY2
      RC = RHO - D0
      XC = RC*COS(PHIC)
      YC = RC*SIN(PHIC)
C
C  Get hits in initial chain
C
      CALL VFILL(JROAD,LTPDRO+2,-1)
      DO 40 M=1,NHCHAI
         JROAD(IRCHAI(M))=IHCHAI(M)
         SROAD(IRCHAI(M))=D2CHAI(M)
 40   CONTINUE
      NXTR1 = 0
C
C   Look at unused padrows for more hits
C

      DO 20 K = IR1,IR2
C
c GAB:
c         NCOR = ITABL(KTCRL,K,JTCRNC)
         NCOR=MY_JTCRNC(K)         

         IF (NCOR.LE.0) GO TO 20
C
C  Don't look at padrow where we already have hit
C
         IF (JROAD(K).NE.-1)  GO TO 20
C
C   Find intersection of orbit with padrow K
C
         RAD = TPPROW(K)
         CALL TRHLX2(RAD,PCHAIN(1),V,IERR)
         IF (IERR.NE.0) GO TO 20
         ZP1 = V(3,1)
C
C   Skip padrow if intersection is unphysical or outside in Z
C
         IF (ABS(ZP1).GT.ZDRIFT+DELZ) GO TO 20
C
C   Get angle-dependent phi tolerance
C

         PHI1 = V(4,1)
         BETA = ATAN2(-(V(1,1)-XC), V(2,1)-YC )
         IF (BETA.LT.0.) BETA=BETA+TWOPI
         BETA = TRPHDI(PHI1,BETA)
         COSB = COS(BETA)
         SINB = SIN(BETA)
         DELW = TDWSTR/ABS(COSB)
         IF (DELW.GT.TDWMAX)  DELW=TDWMAX
c
         NFND1 = 0

C
C   Find the appropriate phi bin, and loop over all hits in it
C
         IBIN= IXWRAP(INT(TCOFFS+PHI1/TC2DEL))
         IPTR= IBNPTR(IBIN,K)
C
C   Search for hits over 1st semicircle
C
C(KH) --- begin
C(KH)  814    IF (IPTR.EQ.0) GO TO 815
         NLOOP=0
c(kh)  814    IF (IPTR.EQ.0) GO TO 815
  814    IF (IPTR.LE.0) GO TO 815
         NLOOP=NLOOP+1
         IF (NLOOP.GT.10000) then
           write(6,*) 'TRPROJ: NLOOP limit reached. abort.'
           GOTO 815
         end if
C(KH) --- end
C GAB:
C            IC= ITABL(INDBIN,IPTR,JNDBCO)
C            JJ = KROW(KTPCO,IC)
C            PHI = RW(JJ + JTPCPH)
C            Z   = RW(JJ + JTPCZV)

            IC= MY_JNDBCO(IPTR)
C            PHI = TPCPHI(IC)
C            Z   = TPCHITS(3,IC)
            PHI = TPCO_PH(IC)
            Z   = TPCO_ZV(IC)

            DPHI = TRPHDI(PHI1,PHI)
            DW = ABS(RAD*DPHI)
            IF (DW.GT.DELW) THEN
C GAB:
C              IPTR=ITABL(INDBIN,IPTR,JNDBFW)
              IPTR=MY_JNDBFW(IPTR)
              GO TO 814
            ENDIF

            DZ = ABS(ZP1-Z)

            IF (DZ.GT.DELZ) THEN
C GAB:
C              IPTR=ITABL(INDBIN,IPTR,JNDBFW)
              IPTR=MY_JNDBFW(IPTR)
              GO TO 814
            ENDIF

            NFND1 = NFND1 + 1
C-----------------------------------------------------------

C GAB: Experimental code
C GAB: Take the best hit, rather than the first one
          CHI2=DW**2/TPCO_SR(IC) +DZ**2/TPCO_SZ(IC)

         IF(NFND1.EQ.1)THEN
          CHI2MIN=CHI2
          ISAV1 = IPTR
C GAB:
C            SAGER2=RW(JJ+JTPCSR)*COSB*COSB + 0.1875*SINB*SINB
            SAGER2=TPCO_SR(IC)*COSB*COSB + PADFAC*SINB*SINB
c          ENDIF


C GAB: Experimental code
         ELSE IF((NFND1.GT.1).AND.(CHI2.LT.CHI2MIN))THEN
          CHI2MIN=CHI2
          ISAV1=IPTR
C GAB:
C            SAGER2=RW(JJ+JTPCSR)*COSB*COSB + 0.1875*SINB*SINB
            SAGER2=TPCO_SR(IC)*COSB*COSB + PADFAC*SINB*SINB
         ENDIF
C------------------------------------

C GAB:
c            IF (NFND1.GT.1) GO TO 815

C
C GAB:
C            IPTR=ITABL(INDBIN,IPTR,JNDBFW)
            IPTR=MY_JNDBFW(IPTR)
            GO TO 814
  815    CONTINUE
C
cC  Accept only 1 candidate in each semicircle
cC
C GAB: Accept the best candidate rather than the first:
C         IF (NFND1.EQ.1) THEN
         IF (NFND1.GT.0) THEN
            NXTR1 = NXTR1 + 1
            JROAD(K) = ISAV1
            SROAD(K) = SAGER2
         ENDIF
C
 20   CONTINUE
C-------------------------------------------------------------------------
c GAB experimental code:
      IF(LUNMRG)THEN
c Only allow a maximum of NSKPMX sequential holes in any chain
       IF (NXTR1.GT.0) THEN
        NSKIP=0
        MMAX=0
        MSTART=0
        CHNSIZ=0
        ISIZ=0
C First get the starting row of the longest ok chain:
        DO M=1,NRTPC
         IF(NSKIP.GT.NSKPMX)THEN
          IF(ISIZ.GT.CHNSIZ)THEN
           CHNSIZ=ISIZ
           MMAX=MSTART
          ENDIF
          MSTART=M
          ISIZ=0
          NSKIP=0
         ENDIF
         IF (JROAD(M).LT.0)THEN
          NSKIP=NSKIP+1
         ELSE
          NSKIP=0
          ISIZ=ISIZ+1
         ENDIF
        ENDDO
        IF(ISIZ.GT.CHNSIZ)THEN
         CHNSIZ=ISIZ
         MMAX=MSTART
        ENDIF

C Now delete any hits which are too disjoint from the biggest chain part:
        DO M=1,MMAX-1
         JROAD(M)=-1
        ENDDO
        DO M=MMAX+CHNSIZ+1,NRTPC
         JROAD(M)=-1
        ENDDO
       ENDIF
      ENDIF
C---------------------------------------------------------------
      IF (NXTR1.GT.0) THEN
         NEW = 0
         DO 50 M = NRTPC,1,-1
            IF (JROAD(M).LT.0) GO TO 50
            NEW = NEW + 1
            IHCHAI(NEW) = JROAD(M)
            IRCHAI(NEW) = M
            D2CHAI(NEW) = SROAD(M)
   50    CONTINUE
         NHCHAI = NEW
C
C  Refit
C
         ISTCHA = 0

         CALL TRFTCH
C
C  If fit is unsuccessful, restore original chain
C
         IF (ISTCHA.EQ.0) THEN
            CALL UCOPY(IRSAV,IRCHAI(1),NSAVE)
            CALL UCOPY(IHSAV,IHCHAI(1),NSAVE)
            CALL UCOPY(MHSAV,MCHAIN(1),NSAVE)
            CALL UCOPY(D2SAV,D2CHAI(1),NSAVE)
            CALL UCOPY(PSAVE,PCHAIN(1),NPARCH)
            ISTCHA = ISTSV
            NHCHAI=NSAVE
            NXTR1 = 0
         ENDIF
      ENDIF
C
C  If extra hits have been found, project again until full length is
C  covered.  Skip curlers.
C
      IF (TURN) GO TO 999
      IF (ISTSV.EQ.ISTCHA) GO TO 999
C
C(KH) --- begin
C(KH)      IF(JSPAN.LT.8)  GO TO 200
      NTRY=NTRY+1
      IF ((jspan.lt.8).and.(ntry.lt.100)) goto 200
      if (ntry.ge.100) write(6,*) 'TRPROJ: NTRY limit. abort.'
C(KH) --- end
C
  999 RETURN
      END

C--------------------------------------------------------------------------

