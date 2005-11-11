*  
* $Id: tknhit.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.2  2004/07/22 15:21:02  aplin
* added tkfunc1.inc
*
* Revision 1.1.1.1  2003/05/23 13:17:51  hvogt
* Brahms V308
*
*  
      INTEGER FUNCTION TKNHIT(CHTYPE,INDEX)
**********************************************************************
*
*     function TKNHIT
*
*     count number of hits in TE/TS/TK
*
*     Kristian Harder, December 1999
*
**********************************************************************
      IMPLICIT NONE
#include "include/fsfparams.inc"
#include "include/tkconst.inc"
#include "include/tkbank.inc"
c#include "gcbank.inc"
      CHARACTER*2 CHTYPE
      INTEGER     INDEX,I
#include "include/tkfunc1.inc"
#include "include/tkfunc.inc"

      TKNHIT=0
      IF (CHTYPE.EQ.'TE') THEN
        IF (INDEX.GE.1.AND.INDEX.LE.NTE) THEN
          TKNHIT=ITEDAT(2,INDEX)
        END IF
      ELSE IF (CHTYPE.EQ.'TS') THEN
        IF (INDEX.GE.1.AND.INDEX.LE.NTS) THEN
          DO I=ITSDAT(1,INDEX),ITSDAT(1,INDEX)+ITSDAT(2,INDEX)-1
            TKNHIT=TKNHIT+ITEDAT(2,ITSTEL(I))
          END DO
        END IF
      ELSE IF (CHTYPE.EQ.'TK') THEN
        IF (INDEX.GE.1.AND.INDEX.LE.NTK) THEN
          DO I=ITKDAT(1,INDEX),ITKDAT(1,INDEX)+ITKDAT(2,INDEX)-1
            TKNHIT=TKNHIT+ITEDAT(2,ITKTEL(I))
          END DO
        END IF
      ELSE
        STOP 'ILLEGAL ARGUMENT IN TKNHIT'
      END IF

      RETURN
      END


