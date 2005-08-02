

*     
* $Id: tpcrun.F,v 1.5 2005-07-29 15:06:17 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.4  2005/07/26 17:56:21  aplin
* Removed unecessary initialisation of hit arrays which is now done in
* tktrev. Corrected (finaly) the calculation of d0 in LEPTrackingProcessor.
*
* Revision 1.3  2005/07/21 16:37:00  aplin
* The itedat bank has been added with read and write interface.
*
* TKADDH has now been added so that the silicon hits are included as
* TE's.
*
* Added FSFSTR from delsearch to steer the creation of
* TS's
*
* The following include files need to be added for FSFSTR
* include/fsfdebug.inc
* include/fsfreject.inc
* include/fsfflags.inc
*
* FSFSTR is now called by TKTREV although all action is commented out.
*
* Revision 1.2  2005/07/19 17:02:02  aplin
* Removed using namespace from .h files. Commented out remove hit from tkbank, as well as adding methods to give information about the number and positions of subdetector hit collections with in the bank structure. Common block CNTPC no longer needed so commented out in tpchitbank.h, and else where it was accessed. Added protocalls for cfortran.h for the new tkbank methods. Replaced TPCRUN as the main tracking steering routine by TKTREV which itself calls TPCRUN.
*
* Revision 1.1.1.1  2005/07/11 15:15:48  aplin
* Initial version
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

      SUBROUTINE TPCRUN(IER)

C****************************************************************************
c* Version number: 1.0      Date of last modification  03.10.99             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
c* TPC track finder, based on the ALEPH track finding code.                 *
c* (many authors from ALEPH involved - see the subroutines for citations)   *
c*                                                                          *
c* Includes background hits and two-hit effects in TPC via routine TPCDIGI  *
c*------------------------------------------------------------------------- *
c*  INPUT: TPCHITS via /TPCHITS/ common block                               *
c*  OUTPUT:  DELPHI-type array                                              *
C* Main Authors:                                                            *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
C#include "include/fndparms.inc"
C#include "turnon.inc"
C#include "include/patrec.inc"
C#include "detnums.inc"
C
C#include "include/my_tcha.inc"
C#include "include/my_ttcc.inc"
C#include "include/my_t1ft.inc"

C
C#include "bkgmode.inc"
C#include "gcbank.inc"
C#include "include/tpcfunc2.inc"
C#include "tpcdims.inc"
C
      IMPLICIT NONE

#include "include/fndparms.inc"
#include "include/tpcdims.inc"
#include "include/tpcfunc2.inc"
#include "include/detnums.inc"
#include "include/tpco.inc"
#include "include/fkparm.inc"
#include "include/bkgmode.inc"
#include "include/chgtrak.inc"
#include "include/tkconst.inc" 
#include "include/tkbank.inc"
#include "include/my_tgft.inc"

      INTEGER ITPTRA(NTPMAX)
      INTEGER ITPDGT(NTPCOMX)

c --- hit information array given back by TKREAD
      REAL HIT(20)
      INTEGER IHIT(20)
      EQUIVALENCE (HIT,IHIT)

      INTEGER TKNUMB

      INTEGER IER,I

      EXTERNAL TKNUMB

      EXTERNAL numofsubdethits
      INTEGER  numofsubdethits

c-----------------------------------------------------------


#include "include/tpcfunc.inc"

      p_flag = 0

c     Initialise the number of tracks found to 0
      NTGFT = 0

c --- get the number of hits in TPC      
CCC SJA the number of hits is set from the MARLIN side

      NTPHITS=numofsubdethits("TPC")

CCC      NTPHITS=TKNUMB('HIT',ID_TPC)
      IF(NTPHITS.GT.NTPMAX)THEN
        IER = 1
        WRITE(*,*)'TPCRUN: Too many hits - increase NTPMAX'
        CALL TKVETO(1)
        GOTO 999
      ENDIF

CCC SJA the hits are filled on the MARLIN side
Cc --- read hits from TK banks and store them locally
C      DO I=1,NTPHITS
C        CALL TKREAD('HIT',ID_TPC,I,HIT,IER)
C        IF ( HIT(4) .EQ. 0. ) HIT(4) = 100.
C        Q(LBRTP1 + 7*(I-1)+1+1) = HIT(1)           ! x
C        Q(LBRTP1 + 7*(I-1)+1+2) = HIT(2)           ! y
C        Q(LBRTP1 + 7*(I-1)+1+3) = HIT(3)           ! z
C        Q(LBRTP1 + 7*(I-1)+1+4) = HIT(4)           ! energy
C        Q(LBRTP1 + 7*(I-1)+1+5) = REAL(ID_TPC)     ! detector ID
C        Q(LBRTP1 + 7*(I-1)+1+6) = HIT(10)          ! rphi resolution
C        Q(LBRTP1 + 7*(I-1)+1+7) = HIT(11)          ! z resolution
C        ITPTRA(I)               = IHIT(6)          ! MC track number
C      END DO
C
     

c --- assign hits to pixels and add noise
c      CALL TPCDIGI(ITPTRA,ITPDGT,IER)
c      IF(IER.NE.0)THEN
c        WRITE(*,*)'TPCRUN: TPCDIGI returned IER=',IER
c        GOTO 999
c      ENDIF

c --- order TPC hits according to (a) padrow and (b) phi within padrow
      CALL THITSRT(ITPDGT,IER)
      IF(IER.NE.0)THEN
        WRITE(*,*)'TPCRUN: THITSRT returned IER=',IER
        GOTO 999
      ENDIF

c --- attempt to remove hits which are probably merged:
      CALL UNMERGE(IER)
      IF(IER.NE.0)THEN
        WRITE(*,*)'TPCRUN: UNMERGE returned IER=',IER
        GOTO 999
      ENDIF

c --- now that merges have been flagged, return to normal track finder:
c --- find links and chains:
      CALL TRKFND(IER)
      IF(IER.NE.0)THEN
        WRITE(*,*)'TPCRUN: TRKFND returned IER=',IER
        GOTO 999
      ENDIF

c --- join up chains where possible:
      CALL TFLNKS(IER)
      IF(IER.NE.0)THEN
        WRITE(*,*)'TPCRUN: TFLNKS returned IER=',IER
        GOTO 999
      ENDIF

c --- assemble track candidates:
      CALL TFCAND(IER)
      IF(IER.NE.0)THEN
        WRITE(*,*)'TPCRUN: TFCAND returned IER=',IER
        GOTO 999
      ENDIF

c --- fit the candidates:
      CALL TFITTK(IER)
      IF(IER.NE.0)THEN
        WRITE(*,*)'TPCRUN: TFITTK returned IER=',IER
        GOTO 999
      ENDIF

c --- link up any broken tracks:
      CALL TFLNK2(IER)
      IF(IER.NE.0)THEN
        WRITE(*,*)'TPCRUN: TFLNK2 returned IER=',IER
        GOTO 999
      ENDIF

c --- transfer track candidates to TK banks
      CALL TPCSAV

 

 999  RETURN
      END

C------------------------------------------------------------------

