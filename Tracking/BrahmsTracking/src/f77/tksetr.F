*  
* $Id: tksetr.F,v 1.3 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.2  2005/08/02 10:32:02  aplin
* removed initialisation of BFIELD from here
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
* Revision 1.1.1.1  2003/05/23 13:17:51  hvogt
* Brahms V308
*
*  
      SUBROUTINE TKSETR
**********************************************************************
*
*     subroutine TKSETR
*
*     track finder initialisation part 2; called by UGINIT after FFGO
*
*     Kristian Harder, September 1999
*
**********************************************************************
      IMPLICIT NONE
#include "include/tkffread.inc"
#include "include/brcntl.inc"
#include "include/fkdebug.inc"
#include "include/turnon.inc"
#include "include/tknoisepar.inc"
#include "include/coildims.inc"

      INTEGER     I,ISTAT
      CHARACTER*20 CHDET(3)
      DATA CHDET  /'Silicon Detectors','TPC','FCH'/



C.....check pattern recognition flags
      DO I=1,3
         
        IF (IPATRC(I).LT.1.OR.IPATRC(I).GT.2) THEN
          WRITE(6,*) 'TKSETR: Wrong mode in PATREC ffread card.'
          stop
        END IF
      END DO

C.....issue pattern recognition status messages
c      WRITE(6,*)
c      WRITE(6,*) '******************************************'
c      WRITE(6,*) 'LOCAL PATTERN RECOGNITION IN SUBDETECTORS:'
c      DO I=1,3
c        IF (IPATRC(I).EQ.1) THEN
c          WRITE(6,*) CHDET(I),': realistic pattern recognition'
c        ELSE
c          WRITE(6,*) CHDET(I),': fake (ideal) pattern recognition'
c        END IF
c      END DO
c      WRITE(6,*) '******************************************'
c      WRITE(6,*)
c
C.....create ntuple and histograms for tracking info if desired
      IF (FSFHIST.GT.0.OR.IHIS.GT.0.OR.TKSTHST.GT.0) THEN
cc        CALL HROPEN(99,'TRKHST','tracking.hst','N',1024,ISTAT)
        IF (ISTAT.NE.0) THEN
          WRITE(6,*) 'error opening histogram file.'
          WRITE(*,*) 'histogram filling for tracking disabled. '
          TKSTHST = 0
          FSFHIST=0
          IHIS=0
        END IF
      END IF

C.....create the relevant subdirectories for tracking histograms
cc      CALL HMDIR('//PAWC/TRACKING',' ')
cc      IF ( TKSTHST .GT. 0 ) CALL HMDIR('//PAWC/TRACKING/TKSTEER',' ')

C.....initialize track fit
      CALL FKINI

C.....initialize track search
      CALL FSFINI

C.....init ambiguity processor
      CALL FXAINI

C.....initialisation of performance analysis
cc      IF (ITKPERF.GT.0) CALL TKIPERF

C.....initialise the VTX patrec
cc      IF (IPATRC(1).EQ.1) CALL SIPSET(TNNBX)

C.....initialise the TPC patrec
      IF (IPATRC(2).EQ.1) CALL TPCSET

C.....setup for noise hits
c      WRITE(6,*) '========================================'
      IF (TNNBX.EQ.0) THEN
c        WRITE(6,*) 'NO PAIR BACKGROUND IN TRACKING DETECTORS'
      ELSE
        WRITE(6,*) 'BACKGROUND SIMULATION ACTIVE.'
        IF (LVXDC) WRITE(6,*) 'CCD VTX INTEGRATES OVER ',TNNBX,' BX'
        IF ((TNENERGY.NE.500).AND.(TNENERGY.NE.800)) THEN
          WRITE(6,*) 'WARNING: NO BACKGROUND SETUP FOR ',
     >                 TNENERGY,' GEV'
          IF (TNENERGY.GT.(800+500)/2) THEN
            TNENERGY=800
          ELSE
            TNENERGY=500
          END IF
          WRITE(6,*) 'USING RATES FOR ',TNENERGY,' GEV'
        END IF
        IF (BFIELD.NE.30.0.AND.BFIELD.NE.40.0) THEN
          WRITE(6,*) 'WARNING: NO BACKGROUND SETUP FOR ',
     >                 BFIELD/10,' T'
          IF (BFIELD.LT.35.) THEN
             WRITE(6,*) 'USING 3T VALUES'
          ELSE
             WRITE(6,*) 'USING 4T VALUES'
          END IF
        END IF
        IF (TNENERGY.EQ.500.AND.BFIELD.LT.35.) THEN
          ! CCD (PER BX)   K. Buesser, Dec 2000
          TNVXL1=548
          TNVXL2=224
          TNVXL3= 76
          TNVXL4= 72
          TNVXL5= 26
          ! APS            crude estimates by K. Harder
          TNVTX1=TNVXL1*14.6/10.
          TNVTX2=TNVXL3*28.6/25.
          TNVTX3=((TNVXL5*28.6/25.)*1./4.+(TNSIT1*28.6/76.)*1./6.)
     >         /(1./4.+1./6.)
          TNVF2A=TNFTD1
          ! FTD            K. Buesser, Dec 2000
          TNFTD1= 71/2
          TNFTD2= 58/2
          TNFTD3= 38/2
          TNFTD4= 24/2
          TNFTD5= 26/2
          TNFTD6= 22/2
          TNFTD7= 12/2
          ! SIT            K. Buesser, Dec 2000
          TNSIT1= 22
          TNSIT2= 15
          ! FCH            K. Buesser, Dec 2000
          TNFCH=19/2
        ELSE IF (TNENERGY.EQ.500) THEN
          ! CCD (PER BX)   K. Buesser, Dec 2000
          TNVXL1=350
          TNVXL2=132
          TNVXL3= 42
          TNVXL4= 26
          TNVXL5= 17
          ! APS            K. Buesser, Sep 2001
          TNVTX1=576
          TNVTX2=98
          TNVTX3=9
          TNVF2A=32
          ! FTD            K. Buesser, Dec 2000
          TNFTD1= 34/2
          TNFTD2= 33/2
          TNFTD3= 24/2
          TNFTD4= 16/2
          TNFTD5= 16/2
          TNFTD6= 13/2
          TNFTD7= 12/2
          ! SIT            K. Buesser, Dec 2000
          TNSIT1= 23
          TNSIT2=  7
          ! FCH            K. Buesser, Dec 2000
          TNFCH=10/2
        ELSE IF (TNENERGY.EQ.800.AND.BFIELD.LT.35.) THEN
          ! CCD (PER BX)   K. Buesser, Dec 2000
          TNVXL1=616
          TNVXL2=264
          TNVXL3= 80
          TNVXL4= 31
          TNVXL5= 23
          ! APS            crude estimates by K. Harder
          TNVTX1=TNVXL1*14.6/10.
          TNVTX2=TNVXL3*28.6/25.
          TNVTX3=((TNVXL5*28.6/25.)*1./4.+(TNSIT1*28.6/76.)*1./6.)
     >         /(1./4.+1./6.)
          TNVF2A=TNFTD1
          ! FTD            K. Buesser, Dec 2000
          TNFTD1= 84/2
          TNFTD2= 68/2
          TNFTD3= 37/2
          TNFTD4= 27/2
          TNFTD5= 32/2
          TNFTD6= 16/2
          TNFTD7= 11/2
          ! SIT            K. Buesser, Dec 2000
          TNSIT1= 16
          TNSIT2=  9
          ! FCH            K. Buesser, Dec 2000
          TNFCH=9/2
        ELSE IF (TNENERGY.EQ.800) THEN
          ! CCD (PER BX)   K. Buesser, Dec 2000
          TNVXL1=422
          TNVXL2=119
          TNVXL3= 57
          TNVXL4= 23
          TNVXL5= 11
          ! APS            K. Buesser, Sep 2001
          TNVTX1=626
          TNVTX2=96
          TNVTX3=4
          TNVF2A=49
          ! FTD            K. Buesser, Dec 2000
          TNFTD1= 54/2
          TNFTD2= 35/2
          TNFTD3= 26/2
          TNFTD4= 26/2
          TNFTD5= 16/2
          TNFTD6= 14/2
          TNFTD7= 11/2
          ! SIT            K. Buesser, Dec 2000
          TNSIT1= 17
          TNSIT2=  5
          ! FCH            K. Buesser, Dec 2000
          TNFCH=11/2
        ELSE
           STOP 'TKNOISE: EXTREMELY STRANGE FFREAD SETUP...'
        END IF
      END IF
c      WRITE(6,*) '========================================'
      WRITE(6,*)

  999 RETURN
      END

