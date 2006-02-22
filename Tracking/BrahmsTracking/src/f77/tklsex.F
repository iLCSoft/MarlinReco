*  
* $Id: tklsex.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.2  2004/07/22 15:21:02  aplin
* added tkfunc1.inc
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      SUBROUTINE TKLSEX(CHTYPE,INDEX,MEXCLU,NEXCLU,IEXCLU,IERR)
**********************************************************************
*
*     subroutine TKLSEX
*
*     get list of TE/TS excluded against given TE/TS
*
*     input:  CHTYPE - 'TE' or 'TS' (upcase letters!)
*             INDEX  - id of TE/TS
*             MEXCLU - dimension of IEXCLU array
*     output: NEXCLU - number of exclusions
*             IEXCLU() ids of excluded TEs/TSs
*             IERR   - 0 if successful
*                      1 if index out of range
*                      2 if MEXCLU too small
*
*     Kristian Harder, September 1999
*
**********************************************************************
      IMPLICIT NONE
#include "include/bkgmode.inc"
#include "include/tkconst.inc"
#include "include/fkparm.inc"
#include "include/tkffread.inc"
#include "include/tkbank.inc"
c#include "gcbank.inc"
      CHARACTER*2 CHTYPE
      INTEGER     INDEX,IERR,MEXCLU,NEXCLU,IEXCLU(*),I
#include "include/tkfunc1.inc"
#include "include/tkfunc.inc"
      IERR=0

      IF (CHTYPE.EQ.'TE') THEN

        IF (INDEX.LT.1.OR.INDEX.GT.NTE) THEN
          IERR=1
          GOTO 999
        END IF

        IF (ITEDAT(3,INDEX).GT.0) THEN
          NEXCLU=ITEDAT(4,INDEX)
          IF (NEXCLU.GT.MEXCLU) THEN
            IERR=2
            IF (TKSTDBG.GT.0) THEN
              WRITE(6,*) 'TKLSEX PROBLEM: ARGUMENT MEXCLU TOO SMALL'
              WRITE(6,*) ' MEXCLU = ',MEXCLU
              WRITE(6,*) ' NUMBER OF EXCLUSION LIST ENTRIES: ',NEXCLU
              DO I=1,NEXCLU
                WRITE(6,*) ' ENTRY ',I,': '
                WRITE(6,*) '   TE NO.: ',IEXCL(ITEDAT(3,INDEX)+I-1)
                WRITE(6,*) '   DETID : ',
     >              ITE(1,IEXCL(ITEDAT(3,INDEX)+I-1))
              END DO
              CALL TKVETO(1)
            END IF
            NEXCLU=MEXCLU
          END IF
          CALL UCOPY(IEXCL(ITEDAT(3,INDEX)),IEXCLU,NEXCLU)
        ELSE
          NEXCLU=0
        END IF

      ELSE IF (CHTYPE.EQ.'TS') THEN

        IF (INDEX.LT.1.OR.INDEX.GT.NTS) THEN
          IERR=1
          GOTO 999
        END IF

        IF (ITSDAT(3,INDEX).GT.0) THEN
          NEXCLU=ITSDAT(4,INDEX)
          IF (NEXCLU.GT.MEXCLU) THEN
            IERR=2
            NEXCLU=MEXCLU
          END IF
          CALL UCOPY(IEXCL(ITSDAT(3,INDEX)),IEXCLU,NEXCLU)
        ELSE
          NEXCLU=0
        END IF

      ELSE
        WRITE(6,*) 'TKLSEX: wrong argument'
        STOP
      END IF

 999  RETURN
      END


