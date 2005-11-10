*  
* $Id: tkglab.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.2  2004/07/22 15:21:02  aplin
* added tkfunc1.inc
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      SUBROUTINE TKGLAB(CHTYPE,INDEX,LABEL,IERR)
**********************************************************************
*
*     subroutine TKGLAB
*
*     get label of TE/TS/TK
*
*     input:  CHTYPE - 'TE','TS', or 'TK' (upcase letters!)
*             INDEX  - Tx identifier
*     output: LABEL  - Tx label
*             IERR   - 0 if successful
*
*     Kristian Harder, September 1999
*
**********************************************************************
      IMPLICIT NONE
#include "include/bkgmode.inc"
#include "include/tkconst.inc"
#include "include/fkparm.inc"
#include "include/tkbank.inc"
c#include "gcbank.inc"
      CHARACTER*2 CHTYPE
      INTEGER     INDEX,LABEL,IERR
#include "include/tkfunc1.inc"
#include "include/tkfunc.inc"

      IERR=0

      IF (CHTYPE.EQ.'TE') THEN
        IF (INDEX.LT.1.OR.INDEX.GT.NTE) THEN
          IERR=1
          GOTO 999
        END IF
        LABEL=ITEDAT(5,INDEX)
      ELSE IF (CHTYPE.EQ.'TS') THEN
        IF (INDEX.LT.1.OR.INDEX.GT.NTS) THEN
          IERR=1
          GOTO 999
        END IF
        LABEL=ITSDAT(5,INDEX)
      ELSE IF (CHTYPE.EQ.'TK') THEN
        IF (INDEX.LT.1.OR.INDEX.GT.NTK) THEN
          IERR=1
          GOTO 999
        END IF
        LABEL=ITKDAT(3,INDEX)
      ELSE
        WRITE(6,*) 'TKGLAB: called with wrong argument'
        STOP
      END IF

 999  RETURN
      END


