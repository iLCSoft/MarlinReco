*  
* $Id: tklist.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.2  2004/07/22 15:21:02  aplin
* added tkfunc1.inc
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      SUBROUTINE TKLIST(CHTYPE,INDEX,NUM,LIST,IERR)
**********************************************************************
*
*     subroutine TKLIST
*
*     get list of all TEs in TS/TK, or of all hits in TE
*
*     input:  CHTYPE - 'TE','TS', or 'TK' (upcase letters!)
*             INDEX  - Tx identifier
*     output: NUM    - number of hits (TE) / TEs (TS/TK)
*             LIST() - identifiers of hits (TE) / TEs (TS/TK)
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
      INTEGER     INDEX,NUM,LIST(*),IERR
#include "include/tkfunc1.inc"
#include "include/tkfunc.inc"
      IERR=0
      NUM=0

      IF (CHTYPE.EQ.'TE') THEN

        IF (INDEX.LT.1.OR.INDEX.GT.NTE) THEN
          IERR=1
          GOTO 999
        END IF

        NUM=ITEDAT(2,INDEX)
        CALL UCOPY(IHITTE(ITEDAT(1,INDEX)),LIST,NUM)

      ELSE IF (CHTYPE.EQ.'TS') THEN

        IF (INDEX.LT.1.OR.INDEX.GT.NTS) THEN
          IERR=1
          GOTO 999
        END IF

        NUM=ITSDAT(2,INDEX)
        CALL UCOPY(ITSTEL(ITSDAT(1,INDEX)),LIST,NUM)

      ELSE IF (CHTYPE.EQ.'TK') THEN

        IF (INDEX.LT.1.OR.INDEX.GT.NTK) THEN
          IERR=1
          GOTO 999
        END IF

        NUM=ITKDAT(2,INDEX)
        CALL UCOPY(ITKTEL(ITKDAT(1,INDEX)),LIST,NUM)

      ELSE

        WRITE(6,*) 'TKLIST: WRONG ARGUMENT'
        IERR=2
        GOTO 999

      END IF

 999  RETURN
      END


