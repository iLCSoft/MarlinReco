*  
* $Id: tklhte.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.2  2004/07/22 15:21:02  aplin
* added tkfunc1.inc
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      SUBROUTINE TKLHTE(INDEX,MEXCLU,NEXCLU,IEXCLU,IERR)
**********************************************************************
*
*     subroutine TKLHTE
*
*     get list of hit TEs contained in current multipoint TE
*
*     input:  INDEX  - (multipoint) TE id
*             MEXCLU - dimension of IEXCLU array
*     output: NEXCLU - number of exclusions
*             IEXCLU() ids of excluded TEs
*             IERR   - 0 if successful
*                      1 if index out of range
*                      2 if MEXCLU too small
*                      3 if TE id INDEX is no multipoint
*
*     Kristian Harder, May 2000; adapted from TKLSEX
*
**********************************************************************
      IMPLICIT NONE
#include "include/bkgmode.inc"
#include "include/tkconst.inc"
#include "include/fkparm.inc"
#include "include/tkffread.inc"
#include "include/tkbank.inc"
#include "include/detnums.inc"
C#include "include/gcbank.inc"
#include "include/tkfunc1.inc"
      INTEGER     INDEX,IERR,MEXCLU,NEXCLU,IEXCLU(*),I,TEDET
#include "include/tkfunc.inc"

      IERR=0
      NEXCLU=0

C.....check range of TE index
      IF (INDEX.LT.1.OR.INDEX.GT.NTE) THEN
        IERR=1
        GOTO 999
      END IF

C.....check if we really have a multipoint TE
      IF ((ITE(1,INDEX).NE.ID_VTX).AND.
     >    (ITE(1,INDEX).NE.ID_FTD).AND.
     >    (ITE(1,INDEX).NE.ID_VTXP)) THEN
        IERR=3
        GOTO 999
      END IF

C.....search for non-multipoints in exclusion list
      IF (ITEDAT(3,INDEX).GT.0) THEN
        DO I=1,ITEDAT(4,INDEX)
          TEDET=ITE(1,IEXCL(ITEDAT(3,INDEX)+I-1))
          IF (TEDET/100*100.NE.TEDET) THEN
            NEXCLU=NEXCLU+1
            IF (NEXCLU.GT.MEXCLU) THEN
              IERR=2
              NEXCLU=MEXCLU
              GOTO 999
            END IF
            IEXCLU(NEXCLU)=IEXCL(ITEDAT(3,INDEX)+I-1)
          END IF
        END DO
      ELSE
        NEXCLU=0
      END IF

 999  RETURN
      END


