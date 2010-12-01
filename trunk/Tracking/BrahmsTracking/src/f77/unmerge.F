*  
* $Id: unmerge.F,v 1.5 2006-02-24 16:36:55 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.4  2006/02/07 10:52:07  owendt
* changed NMAXBL from 4000 to 10000 to ensure stable tracking with Mokka 5.04
*
* Revision 1.3  2005/11/03 15:16:14  aplin
* Added the Trackstring creation and the biulding of full Track candiates (TK's) which have passed the Delphi Ambiguity resolver fxambi. The material description of the vtx detector, as for the TPC, is hard coded in setmat. Presently the VTX and SIT resolutions are hard coded in LEPTrackingProcessor. The debug output has been reduced and can be controlled via TKSTDBG etc. in tkinit.F. delsolve contains the delphi ambuguity resolver written in C and is contained in the directory named C. The Tk's are written back into the C++ side in tktrev. The corresponding Tk bank structure analogous to the TE bank structure has been added in tktkbank whilst the access wrapper functions are contained in LEPTracking.
*
* Revision 1.2  2005/08/08 07:09:13  aplin
* Made f77 tracking code use GEAR to define the geomtery of the TPC. LTPDRO now defines the maximum number of rows is used to define the size of arrays, this is limited to 224 due the use of 7 '32 bit' bit registers in trkfnd.F increased, though at present it is not likely that anybody would want more. The number of TPC padrows is defined at run time by NRTPC which should of course not exceed LTPDRO, although this is checked and the programe exits with a verbose error message. A wrapper function gettpcgeom is used to pass the GEAR TPC parameters from C++ to f77. MarlinUtil/include/marlin_tpcgeom.h have MarlinUtil/src/marlin_tpcgeom.cc consequently been removed as they are no longer needed.
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
      SUBROUTINE UNMERGE(IER)
      IMPLICIT NONE
C****************************************************************************
c* Version number: 0.1      Date of last modification  07.03.00             *
c*------------------------------------------------------------------------- *
c* Main track finding routine                                               *
c* List of modifications since last version:                                *
c*  Modified X/Y/ZBLOT arrays and structure to save on memory               *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
C
#include "include/fndparms.inc"
#include "include/alcons.inc"
#include "include/tpco.inc"
#include "include/tpgeom.inc"
c#include "marlin_tpcgeom.h"
#include "include/tpcdims.inc"
#include "include/tpcres.inc"
#include "include/tphtcns.inc"

#include "include/my_tctc.inc"
#include "include/my_tcha.inc"
C#include "include/geantrk.inc"

C
C----------------------------------------------------------------------
C
C GAB:
      INTEGER IER,I,J,IERR,JJ,IROW,NCMIN
      REAL RAD,V1(4,2),V2(4,2),DPHI,RMIN(MAXCND)
      REAL DZ,DRPH2,TWSG_RP,TWSG_Z,PIX_RP2
      LOGICAL FIRST
      REAL PCHAIN1(5),PCHAIN2(5)


      SAVE TWSG_RP,TWSG_Z
      REAL EOMES,EALPS,EPH0S,EDD0S,EZZ0S,CHI2
      REAL TRPHDI
      EXTERNAL TRPHDI
#include "include/unmrg.inc"

C gab new 7.3.00:
      INTEGER NTOTBL,NMAXBL
C     PARAMETER(NMAXBL=4000)
      PARAMETER(NMAXBL=10000)
      REAL XBLOT(NMAXBL),YBLOT(NMAXBL),ZBLOT(NMAXBL)

      INTEGER NBLOTS(LTPDRO),IBLOT
      REAL BLOTTER,BLOTTER2

c(KH) add index array to find blots for specified pad row
      INTEGER IBLOTS(LTPDRO)

      REAL RMAX1,RMIN1,RMAX2,RMIN2
      real chimin

      SAVE

C-----> error weights scaled with mean chi**2/freedom
c Note: these factors are inverted compared to TFLNKS because
c the inverse error matrix was used there.

      EOMES(JJ)= (MY_JTCHER(JJ,1))/(MY_JTCHC2(JJ)+1.5)
C                           ... FOR 1/R
      EALPS(JJ)= (MY_JTCHER(JJ,3))/(MY_JTCHC2(JJ)+1.5)
C                           ... FOR TAN(LAM)
      EPH0S(JJ)= (MY_JTCHER(JJ,6))/(MY_JTCHC2(JJ)+1.5)
C                           ... FOR PHI(0)
      EDD0S(JJ)= (MY_JTCHER(JJ,10))/(MY_JTCHC2(JJ)+1.5)
C                           ... FOR D(0)
      EZZ0S(JJ)= (MY_JTCHER(JJ,15))/(MY_JTCHC2(JJ)+1.5)
C                           ... FOR Z(0)

      DATA FIRST/.TRUE./

C---------------------------------------------------------------------------

      IF(FIRST)THEN
       FIRST=.FALSE.
       TWSG_RP=2.*TPRP
       TWSG_Z=2.*TPZ
       PIX_RP2=PIX_RP**2
       MINMRG=120
       BLOTTER=1.0
       BLOTTER2=BLOTTER**2
       NCMIN=20
C TMP:
c Dummy code to prevent compile time warning:
       MINROW=MINROW
      ENDIF
      IER=0

C gab new 7.3.00:
      NTOTBL=0

C Quick chain find, without full track finding
      LUNMRG=.TRUE.

c      write(*,*) "about to call trkfnd"

      CALL TRKFND(IER)


      LUNMRG=.FALSE.

CSJAWRT
c      write(*,*) "the number of hits in unmerge = " ,  NTPCO     


      DO I=1,NTPCO
       LBLOT(I)=.FALSE.
      ENDDO

      CALL VZERO(RMIN,MY_NCHAIN)

      chimin=9999.

      DO 5 IROW=1,NRTPC
       RAD = TPPROW(IROW)
       NBLOTS(IROW)=0
c(KH)  memorise where to find blots for this pad row later
       IBLOTS(IROW)=NTOTBL
       DO 10 I=1,MY_NCHAIN-1
        IF(MY_JTCHNC(I).GE.NCMIN)THEN
         RMIN1=TPCO_RV(MY_TCTC(MY_JTCHOF(I)+1))
         RMAX1=TPCO_RV(MY_TCTC(MY_JTCHOF(I)+MY_JTCHNC(I)))

         PCHAIN1(1)=1.0/MY_JTCHIR(I)
         PCHAIN1(2)=MY_JTCHTL(I)
         PCHAIN1(3)=MY_JTCHP0(I)
         PCHAIN1(4)=MY_JTCHD0(I)
         PCHAIN1(5)=MY_JTCHZ0(I)
         CALL TRHLX2(RAD,PCHAIN1(1),V1,IERR)
         IF (IERR.NE.0) GO TO 10
         DO 20 J= I+1,MY_NCHAIN
          IF(MY_JTCHNC(J).GE.NCMIN)THEN
           RMIN2=TPCO_RV(MY_TCTC(MY_JTCHOF(J)+1))
           RMAX2=TPCO_RV(MY_TCTC(MY_JTCHOF(J)+MY_JTCHNC(J)))
c Ensure that there is potential overlap of the chains (to avoid
c blotting out chains belonging to the same track)
           IF((RMIN1.LT.RMAX2).AND.(RMIN2.LT.RMAX1))THEN
            PCHAIN2(1)=1.0/MY_JTCHIR(J)
            PCHAIN2(2)=MY_JTCHTL(J)
            PCHAIN2(3)=MY_JTCHP0(J)
            PCHAIN2(4)=MY_JTCHD0(J)
            PCHAIN2(5)=MY_JTCHZ0(J)
            DPHI=MY_JTCHP0(I)-MY_JTCHP0(J)
            IF (ABS(DPHI) .GT. PI) DPHI=TWOPI-ABS(DPHI)

C Check compatibility of track parameters:
            CHI2 = (MY_JTCHIR(I)
     1                          -MY_JTCHIR(J))**2
     2            *EOMES(I)*EOMES(J)/(EOMES(I)+EOMES(J))
     3            + (MY_JTCHTL(I)
     4                          -MY_JTCHTL(J))**2
     5            *EALPS(I)*EALPS(J)/(EALPS(I)+EALPS(J))
     6             + DPHI**2
     7            *EPH0S(I)*EPH0S(J)/(EPH0S(I)+EPH0S(J))
     8            + (MY_JTCHD0(I)
     9                          -MY_JTCHD0(J))**2
     1            *EDD0S(I)*EDD0S(J)/(EDD0S(I)+EDD0S(J))
     2            + (MY_JTCHZ0(I)
     3                          -MY_JTCHZ0(J))**2
     4            *EZZ0S(I)*EZZ0S(J)/(EZZ0S(I)+EZZ0S(J))

            IF(CHI2.GT.15.)THEN
             CALL TRHLX2(RAD,PCHAIN2(1),V2,IERR)
             IF (IERR.NE.0) GO TO 20

             DZ = ABS(V1(3,1)-V2(3,1))
             DRPH2=(V1(1,1)-V2(1,1))**2+(V1(2,1)-V2(2,1))**2

             IF((DRPH2.LT.PIX_RP2*BLOTTER2).AND.
     &          (DZ.LT.PIX_Z*BLOTTER)) THEN
cc Candidates for hit merging:
              NBLOTS(IROW)=NBLOTS(IROW)+1
C gab new 7.3.00:
              NTOTBL=NTOTBL+1
              IF (NTOTBL.GT.NMAXBL)THEN
               WRITE(*,*)'UNMERGE: too many blots => increase NMAXBL'
c               CALL TKVETO(1)
               GO TO 20
              ENDIF
C gab new 7.3.00:
              XBLOT(NTOTBL)=(V1(1,1)+V2(1,1))/2.
              YBLOT(NTOTBL)=(V1(2,1)+V2(2,1))/2.
              ZBLOT(NTOTBL)=(V1(3,1)+V2(3,1))/2.
             ENDIF
            ENDIF
           ENDIF
          ENDIF
20       ENDDO
        ENDIF
10     ENDDO
5     ENDDO


      DO I=1,NTPCO
        IROW=1+ INT((TPCO_RV(I)-TPCACRI)/TPCPADR )
        IF ( IROW .LE. 0 ) IROW = 1
        IF ( IROW .GE. NRTPC ) IROW = NRTPC
c(KH) pointer just before first blot in current pad row
        NTOTBL=IBLOTS(IROW)
        DO IBLOT=1,NBLOTS(IROW)
C gab new 7.3.00:
        NTOTBL=NTOTBL+1
c(kh)
        if (ntotbl.le.nmaxbl) then
          IF(ABS(TPCO_ZV(I)-ZBLOT(NTOTBL)).LT.blotter*PIX_Z/2.)THEN
           IF((TPCO_XV(I)-XBLOT(NTOTBL))**2+
     &        (TPCO_YV(I)-YBLOT(NTOTBL))**2
     &                      .LT.blotter2*PIX_RP2/4.)THEN
            IF(ABS(TPCO_ZV(I)-ZBLOT(NTOTBL)).LT.blotter*PIX_Z/2.)THEN
                  LBLOT(I)=.TRUE.
            ENDIF
           ENDIF
          ENDIF
        endif
       ENDDO
      ENDDO

  999 RETURN
      END

C------------------------------------------------------------------

