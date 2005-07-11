*  
* $Id: fksref.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:55  hvogt
* Brahms V308
*
*  
************************************************************************
*                                                                      *
      SUBROUTINE FKSREF(DATS)
*                                                                      *
*   Computation of parameters of the reference trajectory              *
*   PARREF (cylinder param.) and PARPLA (plane param.)                 *
*                                                                      *
*   Input  :  DATS   :  TSR data array                                 *
*                                                                      *
************************************************************************
#include "include/fksave.inc"
#include "include/fkdebug.inc"
#include "include/fkclun.inc"
#include "include/fkpipar.inc"
#include "include/fkfild.inc"
#include "include/fkparm.inc"
#include "include/fkcom.inc"
#include "include/fkddes.inc"
*
      DIMENSION DATS(MXTS)
*
*  z , theta and 1/R
*
      PARREF(3)=DATS(19)
      PARPLA(1)=DATS(19)
      PARREF(4)=DATS(20)
      PARPLA(4)=DATS(20)
      SINTH=SIN(DATS(20))
      COSTH=COS(DATS(20))
      PINV=DATS(22)
      PARREF(6)=CONSB*PINV/SINTH
CSJAWRT
c      write(*,*) "CONSB = " , CONSB
c      write(*,*) "PINV = " , PINV
c      write(*,*) "SINTH = " , SINTH
 
*     Protect against unphysical large momentum       21-09-94
      IF(DABS(PARREF(6)).LT.1.D-10) PARREF(6)=1.D-10
      PARPLA(6)=PARREF(6)
*
*   if the reference surface is a cylinder
*
      IF(JBIT(DATS(3),2).EQ.1) THEN
        PARREF(1)=DATS(17)
        PARREF(2)=DATS(18)/DATS(17)
        PARPLA(2)=PARREF(1)*DCOS(PARREF(2))
        PARPLA(3)=PARREF(1)*DSIN(PARREF(2))
*
*   if the reference surface is a plane
*
      ELSE
        PARREF(1)=SQRT(DATS(17)**2+DATS(18)**2)
*
*   protection for TSR with x=y=0 (SD)
*
        IF((ABS(DATS(17)).LT.1E-10).AND.(ABS(DATS(18)).LT.1E-10)) THEN
          PARREF(2)=DATS(21)
        ELSE
cKMO          PARREF(2)=AMOD(ATAN2(DATS(18),DATS(17))+TWOPI,TWOPI)
          PARREF(2)=ATAN2(DATS(18),DATS(17))
          IF (PARREF(2).LT.0.) PARREF(2) = PARREF(2)+TWOPI
        ENDIF
        PARPLA(2)=DATS(17)
        PARPLA(3)=DATS(18)
*???    IF(PAR(1).LT.0.) PAR(6)=-PAR(6)
      ENDIF
*
*   beta or phi
      ANG=PARREF(2)
CKMO      PARREF(5)=AMOD(DATS(21)-ANG+TWOPI+PI,TWOPI)-PI
      PARREF(5)=DATS(21)-ANG
      IF (PARREF(5).LT.-PI) PARREF(5) = PARREF(5)+TWOPI
      IF (PARREF(5).GT.PI) PARREF(5) = PARREF(5)-TWOPI
      PARPLA(5)=DATS(21)
*
      if (IDEB.ge.2) then
        write(IPRNT,
     $       '(''      PARREF(1) = '',e20.12/'//
     $       ' ''      PARREF(2) = '',e20.12/'//
     $       ' ''      PARREF(3) = '',e20.12/'//
     $       ' ''      PARREF(4) = '',e20.12/'//
     $       ' ''      PARREF(5) = '',e20.12/'//
     $       ' ''      PARREF(6) = '',e20.12)') (parref(i),i=1,6)
        write(IPRNT,
     $       '(''      PARPLA(1) = '',e20.12/'//
     $       ' ''      PARPLA(2) = '',e20.12/'//
     $       ' ''      PARPLA(3) = '',e20.12/'//
     $       ' ''      PARPLA(4) = '',e20.12/'//
     $       ' ''      PARPLA(5) = '',e20.12/'//
     $       ' ''      PARPLA(6) = '',e20.12)') (parpla(i),i=1,6)
      end if

      END

