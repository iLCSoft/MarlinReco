*  
* $Id: thitsrt.F,v 1.2 2005-08-08 07:09:13 aplin Exp $
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
      PARAMETER(NMXHTRW=600)
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
      IF (NRTPC.GT.NMAXRW)THEN
       WRITE(*,*)'THITSRT: NMAXRW is too small'
       IER=2
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
         IF((ABS(TPCHITS(4,IHIT)).GT.TPTHRS).AND.
     &        (ABS(TPCHITS(4,IHIT)).LT.TPSAT))THEN
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
         else
c            write(*,*) "####### threshold not met IHIT =" , IHIT,
c     +           "threshold = ", TPCHITS(4,IHIT)
         ENDIF
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
      write(*,*) "the number of hits at end of thitsrt = " ,  NTPCO     


      RETURN
      END

C----------------------------------------------------------------------

