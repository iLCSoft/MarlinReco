*  
* $Id: tkread.F,v 1.2 2005-07-20 17:46:38 aplin Exp $
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
* Revision 1.2  2004/07/22 15:21:02  aplin
*  added tkfunc1.inc
*
* Revision 1.1.1.1  2003/05/23 13:17:51  hvogt
* Brahms V308
*
*  
      SUBROUTINE TKREAD(CHTYPE,ID_DET,INDEX,DATA,IDATA,IERR)
**********************************************************************
*
*     subroutine TKREAD
*
*     read hit/TE/TS/TK from tracking bank arrays
*
*     input:  CHTYPE - 'MCT','HIT','TE','TS', or 'TK' (upcase letters!)
*             INDEX  - index of object to be read
*             ID_DET - only used for 'HIT' - subdetector ID to read out
*                      (if 0, linear access to full bank is used)
*     output: DATA() - track/hit/TE/TS/TK (array size 11/MXHT/MXTE/MXTS/MXTK)
*             IERR   - 0 if ok
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
      INTEGER     ID_DET,INDEX,IERR,MAXN,INDEXT,MM
      REAL        DATA(*)
      INTEGER     IDATA(*)
 
      REAL RMARINTF(MXHT)
      INTEGER IMARINTF(MXHT)


      EQUIVALENCE (RMARINTF,IMARINTF)

      CHARACTER*3 CHTYPE
#include "include/tkfunc1.inc"
#include "include/tkfunc.inc"
      IERR=0

      IF (CHTYPE(1:2).EQ.'TE') THEN
        IF (INDEX.LT.1.OR.INDEX.GT.NTE) THEN
          IERR=1
          GOTO 999
        END IF
C(KH)   CALL UCOPY(RTE(1,INDEX),DATA,MXTE)
C FIXME: SJA: zebra commented out
c        CALL UCOPY(Q(LBRTE+MXTE*(INDEX-1)+1),DATA,MXTE)
      ELSE IF (CHTYPE(1:2).EQ.'TS') THEN
        IF (INDEX.LT.1.OR.INDEX.GT.NTS) THEN
          IERR=1
          GOTO 999
        END IF
        CALL UCOPY(RTS(1,INDEX),DATA,MXTS)
      ELSE IF (CHTYPE(1:2).EQ.'TK') THEN
        IF (INDEX.LT.1.OR.INDEX.GT.NTK) THEN
          IERR=1
          GOTO 999
        END IF
        CALL UCOPY(RTK(1,INDEX),DATA,MXTK)
      ELSE IF (CHTYPE(1:3).EQ.'MCT') THEN
        IF (INDEX.LT.1.OR.INDEX.GT.TKNTRK) THEN
          IERR=1
          GOTO 999
        END IF
*        CALL UCOPY(TKMCTR(1,INDEX),DATA,11)
        DO MM=1, 11
C FIXME: SJA: need to reinclude MC tracks at some point           
ccc          DATA(MM) = TKMCTR(MM,INDEX)
        END DO
      ELSE IF (CHTYPE(1:3).EQ.'HIT') THEN

         IF (ID_DET.EQ.0) THEN
            MAXN=NTKHIT
         ELSE IF (ID_DET.GT.0.AND.ID_DET.LE.3000) THEN
            MAXN=IHNUMB(ID_DET)
         ELSE
            MAXN=-1
         END IF
         IF (INDEX.LT.1.OR.INDEX.GT.MAXN) THEN
            IERR=1
            GOTO 999
         END IF
         IF (ID_DET.GT.0) THEN
            INDEXT=INDEX+IHPOINT(ID_DET)-1
         ELSE
            INDEXT=INDEX
         END IF
*     *        CALL UCOPY(RTKHIT(1,INDEXT),DATA,MXHT)

c         DO MM=1,MXHT
c            DATA(MM) = RTKHIT(MM,INDEXT)
c         END DO
         
C SJA: DO loop above replace with individual lines due to problem
C     with integer float equivalence when replacing IQ, Q 
C     zebra write statements with marlin interface
         
         DATA(1) = RTKHIT(1,INDEXT)
c         WRITE(*,*) "1ST DATA(1) = ",DATA(1)
         DATA(2) = RTKHIT(2,INDEXT)
c         WRITE(*,*) "1ST DATA(2) = ",DATA(2)
         DATA(3) = RTKHIT(3,INDEXT)
c         WRITE(*,*) "1ST DATA(3) = ",DATA(3)
         DATA(4) = RTKHIT(4,INDEXT)
c         WRITE(*,*) "1ST DATA(4) = ",DATA(4)

c         IMARINTF(5) = INT(RTKHIT(5,INDEXT))
c         DATA(5) = RMARINTF(5)
c         DATA(5) = ITKHIT(5,INDEXT)
c         WRITE(*,*) "2ND DATA(5) = ",DATA(5)
c         IMARINTF(6) = INT(RTKHIT(6,INDEXT))
c         DATA(6) = RMARINTF(6)
c         IMARINTF(7) = INT(RTKHIT(7,INDEXT))
c         DATA(7) = RMARINTF(7)
c         IMARINTF(8) = INT(RTKHIT(8,INDEXT))
c         DATA(8) = RMARINTF(8)
c         IMARINTF(9) = INT(RTKHIT(9,INDEXT))
c         DATA(9) = RMARINTF(9)
c
         IDATA(5) = ITKHIT(5,INDEXT)
c         WRITE(*,*) "1ST DATA(5) = ",IDATA(5)
         IDATA(6) = ITKHIT(6,INDEXT)
c         WRITE(*,*) "1ST DATA(6) = ",IDATA(6)
         IDATA(7) = ITKHIT(7,INDEXT)
c         WRITE(*,*) "1ST DATA(7) = ",IDATA(7)
         IDATA(8) = ITKHIT(8,INDEXT)
c         WRITE(*,*) "1ST DATA(8) = ",IDATA(81)
         IDATA(9) = ITKHIT(9,INDEXT)
c         WRITE(*,*) "1ST DATA(9) = ",IDATA(9)

         DATA(10) = RTKHIT(10,INDEXT)
         DATA(11) = RTKHIT(11,INDEXT)


      ELSE
        WRITE(6,*) 'TKREAD: called with wrong argument'
        STOP
      END IF

 999  RETURN
      END


