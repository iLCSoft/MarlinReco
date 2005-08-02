*  
* $Id: tknumb.F,v 1.1 2005-08-01 08:17:47 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.2  2004/07/22 15:21:02  aplin
* added tkfunc1.inc
*
* Revision 1.1.1.1  2003/05/23 13:17:51  hvogt
* Brahms V308
*
*  
      INTEGER FUNCTION TKNUMB(CHTYPE,ID_DET)
**********************************************************************
*
*     function TKNUMB
*
*     return total number of hits/TE/TS/TK from tracking bank arrays
*
*     input:  CHTYPE - 'MCT','HIT','TE','TS', or 'TK' (upcase letters!)
*             ID_DET - only used for 'HIT'. specify subdetector ID
*                      for which number of hits is to be returned.
*                      ID_DET=0 -> all hits
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
      CHARACTER*3 CHTYPE
      INTEGER     ID_DET
#include "include/tkfunc1.inc"
#include "include/tkfunc.inc"

      IF (CHTYPE(1:2).EQ.'TE') THEN
        TKNUMB=NTE
      ELSE IF (CHTYPE(1:2).EQ.'TS') THEN
        TKNUMB=NTS
      ELSE IF (CHTYPE(1:2).EQ.'TK') THEN
        TKNUMB=NTK
      ELSE IF (CHTYPE(1:3).EQ.'MCT') THEN
         TKNUMB=TKNTRK
      ELSE IF (CHTYPE(1:3).EQ.'HIT') THEN
        IF (ID_DET.EQ.0) THEN
          TKNUMB=NTKHIT
        ELSE IF (ID_DET.GT.0.AND.ID_DET.LE.3000) THEN
          TKNUMB=IHNUMB(ID_DET)
        ELSE
          WRITE(6,*) 'TKNUMB: INVALID DETECTOR ID'
          TKNUMB=-1
        END IF
      ELSE
        WRITE(6,*) 'TKNUMB: WRONG ARGUMENT.'
        TKNUMB=-1
      END IF

 999  RETURN
      END


