*  
* $Id: thitsrt.F,v 1.6 2008-01-24 17:21:50 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.5  2006/02/24 16:36:55  aplin
* increased NMXHTRW again by a factor of 2 in thitsrt.F
* increased the minimum number of hits to form a chain in trkfnd.F from MNLNK=4 to MNLNK=7
* set the minimum radius/pad ring for inward search in first track finding to MINMRG=120 in unmerge.F
* reduced the maximum number of links NMAXLNK to 25000
*
* Revision 1.4  2006/02/07 10:50:57  owendt
* changed NMXHTRW from 600 to 1200 to ensure stable tracking with Mokka 5.04
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
* Revision 1.2  2003/12/15 14:57:08  hvogt
* lcio I/O included
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      SUBROUTINE THITSRT(ITPTRA,IER)
C****************************************************************************
c* Version number: 0.1      Date of last modification   6.11.99             *
c*------------------------------------------------------------------------- *
c* Routine to sort the TPC hits by row and then by phi within the row       *
c* List of modifications since last version:                                *
C* Main Authors:                                                            *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c* K.Moenig      Zeuthen                                                    *
c*                                                                          *
c****************************************************************************

      IMPLICIT NONE

      INTEGER IER,ITPTRA(*)

#include "include/fndparms.inc"
#include "include/alcons.inc"
#include "include/tpcdims.inc"
#include "include/tpco.inc"
#include "include/tphtcns.inc"
c#include "include/geantrk.inc"
#include "include/tpgeom.inc"
c#include "marlin_tpcgeom.h"
#include "include/tpcfunc2.inc"

      LOGICAL LUNSORT(NTPCOMX)

      INTEGER I,J,IHIT,IROW,ISORT,IOLD
      INTEGER IUCOMP
      REAL RHIT(NTPCOMX),PHI(NTPCOMX),MINPHI
      INTEGER NMAXRW,NMXHTRW
      PARAMETER(NMAXRW=LTPDRO)
ccc      PARAMETER(NMXHTRW=300)
* increase by a factor of 2 (Mokka input)
ccc      PARAMETER(NMXHTRW=600)
* increase again by a factor of 2 because of the relase of the energy cut in Mokka
c SJA * increased again by a factor of 2 because of the relase of the energy cut in Mokka
      PARAMETER(NMXHTRW=2400)
      INTEGER NINROW(NMAXRW),IHITRW(NMAXRW,NMXHTRW),NSORTED

      REAL ZRES2_PIX
      REAL RPRES2_PIX
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      SAVE

C------------------------------------------
c#include "bkgmode.inc"
c#include "gcbank.inc"
#include "include/tpcfunc.inc"

      IF(FIRST)THEN
       FIRST=.FALSE.
       ZRES2_PIX=(2.*PIX_Z)**2/12.
       RPRES2_PIX=(2.*PIX_RP)**2/12.
      ENDIF

      IER=0

C quick checks:
      IF(NTPHITS.GT.NTPCOMX)THEN
       WRITE(*,*)'THITSRT: too many TPC hits'
CCC       CALL TKVETO(1)
       IER=1
       RETURN
      ENDIF
      IF (NRTPC.GT.LTPDRO)THEN
       WRITE(*,*)'THITSRT: the number of TPC padrows specified in the' 
       WRITE(*,*)'GEAR file is larger than the max number of rows '
       WRITE(*,*)'that the code can handle i.e. LTPDRO'
       WRITE(*,*)'increase N32BITREG in ./src/f77/include/padrow.inc'
       WRITE(*,*)'gmake clean and recompile'
       IER=99
       RETURN
      ENDIF

*
* count number of GEANT tracks
*
* K Moenig improvement:
CCC      ngntrk = 0
CCC      iold   = itptra(1)
CCC      if (iold.gt.0) then
CCC         ngntrk   = 1
CCC         labtn(1) = itptra(1)
CCC      endif
CCC      do ihit=2,ntphits
CCC         if (itptra(ihit).gt.0) then
CCC            if (itptra(ihit).ne.iold) then
CCC               iold = itptra(ihit)
CCC               if (iucomp(itptra(ihit),labtn,ngntrk).eq.0) then
CCC                  ngntrk        = ngntrk+1
CCC                  labtn(ngntrk) = itptra(ihit)
CCC               endif
CCC            endif
CCC         endif
CCC      enddo

      DO I=1,NTPHITS
       LUNSORT(I)=.TRUE.
      ENDDO

      DO I=1,NRTPC
       NINROW(I)=0
       DO J=1,NMXHTRW
        IHITRW(I,J)=0
       ENDDO
      ENDDO

      NTPCO=0
      DO IHIT=1,NTPHITS
c         IF((ABS(TPCHITS(4,IHIT)).GT.TPTHRS).AND.
c     &        (ABS(TPCHITS(4,IHIT)).LT.TPSAT))THEN
         RHIT(IHIT)=SQRT(TPCHITS(1,IHIT)**2+TPCHITS(2,IHIT)**2)
         PHI(IHIT)=ATAN2(TPCHITS(2,IHIT),TPCHITS(1,IHIT))
         IF (PHI(IHIT).LT.0.)PHI(IHIT)=PHI(IHIT)+TWOPI
         IROW=1+INT((RHIT(IHIT)-TPCACRI)/TPCPADR)
         IF (IROW.LT.1) IROW=1
         IF (IROW.GT.NRTPC) IROW=NRTPC
C     first sort according to rows:
         NINROW(IROW)=NINROW(IROW)+1
         IF(NINROW(IROW).GT.NMXHTRW)THEN
            WRITE(6,*) 'THITSRT error: too many TPC hits in one row'
CCC   CALL TKVETO(1)
            NINROW(IROW)=NMXHTRW
         ENDIF
         IHITRW(IROW,NINROW(IROW))=IHIT
C     SJAWRT
c      else
c            write(*,*) "####### threshold not met IHIT =" , IHIT,
c     +           "threshold = ", TPCHITS(4,IHIT)
c         ENDIF
      ENDDO
      
c Now sort within rows according to phi:
      


      DO IROW=1,NRTPC
         NSORTED=0
         DO WHILE(NSORTED.LT.NINROW(IROW))
            MINPHI=TWOPI
            DO I=1,NINROW(IROW)
               IHIT=IHITRW(IROW,I)
               IF(LUNSORT(IHIT))THEN
                  IF(PHI(IHIT).LT.MINPHI)THEN
                     ISORT=IHIT
                     MINPHI=PHI(IHIT)
                  ENDIF
               ENDIF
            ENDDO

            NTPCO=NTPCO+1
            TPCO_XV(NTPCO)=TPCHITS(1,ISORT)
            TPCO_YV(NTPCO)=TPCHITS(2,ISORT)
            TPCO_ZV(NTPCO)=TPCHITS(3,ISORT)
            TPCO_RV(NTPCO)=RHIT(ISORT)
            TPCO_PH(NTPCO)=PHI(ISORT)
c     Note that merged hits (flagged by negative amplitudes) receive
c     a resolution given by the 'pixel' size.
            IF(TPCHITS(4,ISORT).GT.0.)THEN
               TPCO_SR(NTPCO)=TPCHITS(6,ISORT)**2
               TPCO_SZ(NTPCO)=TPCHITS(7,ISORT)**2
            ELSE IF(TPCHITS(4,ISORT).LT.0.)THEN
               TPCO_SR(NTPCO)=RPRES2_PIX
               TPCO_SZ(NTPCO)=ZRES2_PIX
            ENDIF
C     FIXME: relplaced zebra bank with array HITLNK
            HITLNK(NTPCO) = ISORT
c     IQ(LBRTP2+NTPCO) = ISORT
c     IGNTRK(NTPCO)=ITPTRA(ISORT)
            LUNSORT(ISORT)=.FALSE.
            NSORTED=NSORTED+1
         ENDDO
      ENDDO
      
C     SJAWRT
c      write(*,*) "the number of hits at end of thitsrt = " ,  NTPCO     


      RETURN
      END

C----------------------------------------------------------------------

