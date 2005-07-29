*  
* $Id: tpcsav.F,v 1.2 2005-07-29 15:06:17 aplin Exp $
*  
* $Log: not supported by cvs2svn $
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
      SUBROUTINE TPCSAV
C****************************************************************************
c* Version number: 0.1      Date of last modification  03.10.99             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
c*    03.10.99  -  adapted to second version of patrec interface (KH)       *
c* Simple conversion to DELPHI type output                                  *
c*------------------------------------------------------------------------- *
c*  INPUT: ALEPH format tracks in /TGFT/ common                             *
c*  OUTPUT:  DELPHI-type array                                              *
C* Main Authors:                                                            *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
      IMPLICIT NONE
#include "include/fndparms.inc"
#include "include/my_tgft.inc"
#include "include/my_tgtl.inc"
#include "include/my_tgcl.inc"
#include "include/tpco.inc"
#include "include/detnums.inc"
#include "include/tpcdbg.inc"
#include "include/alcons.inc"
#include "include/coildims.inc"

C Local variables
      INTEGER I,I2,J,ITK,NHITS,IERR,LIST1(NTPCOMX),LIST2(NTPCOMX)
      REAL    REF(6),PHIREF
      LOGICAL SKIP

C----------------------------
#include "include/bkgmode.inc"
c#include "include/gcbank.inc"
#include "include/tpcfunc2.inc"
#include "include/tpcfunc.inc"


      write(*,*) "Number of tracks found = " , NTGFT
      DO ITK=1,NTGFT

c ---  skip if number of hits is zero
C(KH)  IF (MY_JTGTN1(ITK)+MY_JTGTNR(ITK).EQ.0) THEN
       IF (MY_JTGTN1(ITK).EQ.0) THEN
         IF (TPCDBG.GT.0) THEN
           WRITE(6,*) 'TPCSAV: track ',itk,' of ',ntgft,': no hits.'
           WRITE(6,*) '         skipping it.'
         END IF
         GOTO 900
       END IF

c ---  list of hits
C(KH)  NHITS=MY_JTGTN1(ITK)+MY_JTGTNR(ITK)
       NHITS=MY_JTGTN1(ITK)
       I2=0
       DO I=1,NHITS
         SKIP=.FALSE.
         DO J=1,I-1
 
            IF (HITLNK(MY_TGCLT(MY_JTGTIO(ITK)+I)).EQ.
     >         HITLNK(MY_TGCLT(MY_JTGTIO(ITK)+J))) SKIP=.TRUE.
         END DO
         IF (.NOT.SKIP) THEN
           I2=I2+1
           LIST1(I2)=HITLNK(MY_TGCLT(MY_JTGTIO(ITK)+I))
          LIST2(I2)=ID_TPC
         END IF
       END DO

c ---  check for strange results
       IF (I2.NE.NHITS.AND.TPCDBG.GT.0) THEN
         WRITE(6,*) 'TPCSAV: only ',I2,' unique hits out of ',NHITS
       END IF

c ---  transform track fit results into DELPHI track parameters as reference
       REF(1)=ABS(MY_JTGFD0(ITK))                ! r_ref
       PHIREF=MY_JTGFP0(ITK)+PIBY2
       IF (PHIREF.GT.TWOPI) PHIREF=PHIREF-TWOPI
       REF(2)=ABS(MY_JTGFD0(ITK))*PHIREF         ! r_ref*phi_ref
       REF(3)=MY_JTGFZ0(ITK)                     ! z_ref
       REF(4)=-ATAN(MY_JTGFTL(ITK))+PIBY2        ! theta
       REF(5)=MY_JTGFP0(ITK)                     ! phi0
       REF(6)=MY_JTGFIR(ITK)/FCT/BFIELD*SIN(REF(4))  ! signed 1/p

CSJAWRT
       write(*,*) "r_ref = ",REF(1) 
       write(*,*) "r_ref*phi_ref = ",REF(2) 
       write(*,*) "z_ref = ",REF(3) 
       write(*,*) "theta = ", REF(4) 
       write(*,*) "phi0  = ",REF(5) 
       write(*,*) "momentum = " , 1/REF(6)

c ---  final fit and track storage
       CALL TKMKTE(ID_TPC,I2,LIST1,LIST2,REF,IERR)

       IF (IERR.NE.0) THEN
          IF (TPCDBG.GT.0) THEN
             WRITE(6,*) 'TPCSAV: TKMKTE FAILED.'
             WRITE(6,*) '  number of hits was ',I2
          END IF
       END IF

 900  ENDDO

      RETURN
      END


C----------------------------------------------------------------------

