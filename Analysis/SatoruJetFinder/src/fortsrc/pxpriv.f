CDECK  ID>, PXPRIV.
      SUBROUTINE PXPRIV (NAME,NELMT,IARR)
*.*********************************************************
*. ------
*. PXPRIV
*. ------
*. SOURCE: J.W.Gary
*. Print-out of an integer vector array
*. Usage:
*.
*.      CHARACTER*6 NAME
*.      INTEGER  NELMT
*.      PARAMETER  (NELMT=1.or.more)
*.      INTEGER IARR (NELMT)
*.      INTEGER  ISIZ
*.
*.      NAME = 'IARR  '
*.      ISIZ = 1.to.NELMT
*.      CALL PXPRIV (NAME,ISIZ,IARR)
*.
*. INPUT     : NAME    The six character name of the vector
*. INPUT     : ISIZ    The number of elements to print
*. INPUT     : IARR    The integer vector array
*.
*.*********************************************************
      IMPLICIT NONE
      INTEGER  IARR (*),IARREM (10)
      INTEGER  NELMT,NROW,IR,IC,NCOL,NREM
      CHARACTER*(*)  NAME
      NROW =  NELMT / 10
      NREM = MOD (NELMT,10)
      IF (NREM.NE.0) THEN
          DO 110 IC = 1,10
              IARREM (IC) = 0
              IF (IC.LE.NREM) IARREM (IC) = IARR (NROW*10+IC)
 110      CONTINUE
          NROW = NROW + 1
      END IF
      NCOL = 10
      WRITE (6,FMT='('' Array name: '',A6)') NAME
      DO 120 IR = 1,NROW
          IF (IR.EQ.NROW.AND.NREM.NE.0) THEN
              NCOL = NREM
              WRITE (6,FMT='(I4,'':  '',10I6)')
     +          IR,(IARREM (IC),IC=1,NCOL)
          ELSE
              WRITE (6,FMT='(I4,'':  '',10I6)')
     +          IR,(IARR ((IR-1)*10+IC),IC=1,NCOL)
          END IF
 120  CONTINUE
      RETURN
      END
