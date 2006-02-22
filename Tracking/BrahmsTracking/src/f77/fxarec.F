*  
* $Id: fxarec.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
       SUBROUTINE FXAREC(NEWTK,IDNEWTK,IFXDBG,IFAULT)
*
************************************************************************
*      Lift TS banks for tracks modified in the ambiguity solver.
*      ID of list of new TS in IDLIST. (TS lifted as DEACTIVE banks)
*      J.Wickens 25/3/96
************************************************************************
*
       IMPLICIT NONE

       INTEGER IFXDBG,I,J,IDET,IERR,ID
       INTEGER IFAULT,NEWTK,IDNEWTK(*)
*
#include "include/bkgmode.inc"
#include "include/tkconst.inc"
#include "include/fxapar.inc"
*
 2001 FORMAT(1X,'FXAREC:',2I5,I7,1X,I8)
*
      IFAULT=0
      NEWTK=0
*
      DO 1 I=1,NRET
        ID=IRETLI(I)
        IF(IDTKR(ID).GT.0)GO TO 1
*
*      Build used module identifier
*
        DO J=1,NTERUS(ID)
          IDET=IMODTE(J,ID)
          CALL SBIT1(MODUSE(ID),IDET)
        ENDDO
*
*      Lift new TS bank
*      Use data from original TKR as start point for fitting
*
        CALL TKSTOR ('TK',IDNEWTK(NEWTK+1),ITKDAT(1,ID),
     >               NTERUS(ID),IDTERU(1,ID),0,IERR)
        IF (IERR.NE.0) THEN
          CALL SBIT1(IFAULT,2)
          GO TO 1
        ELSE
          NEWTK = NEWTK + 1
        ENDIF
*
        IF(IFXDBG.GT.0)THEN
          WRITE(6,2001)I,ID,NEWTK,MODUSE(ID)
        ENDIF
 1    CONTINUE
*
      RETURN
      END

