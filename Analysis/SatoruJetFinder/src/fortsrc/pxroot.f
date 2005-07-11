CDECK  ID>, PXROOT.
      SUBROUTINE PXROOT (XOPER,XROOT)
*.*********************************************************
*. ------
*. PXROOT
*. ------
*. SOURCE: HERWIG (B.Webber,G.Marchesini)
*. Square root with sign retention
*. Usage     :
*.
*.      REAL  XOPER,XROOT
*.
*.      CALL PXROOT (XOPER,XROOT)
*.
*. INPUT     : XOPER   The input operand
*. OUTPUT    : XROOT   Square root of the operand
*.
*.*********************************************************
      IMPLICIT NONE
      REAL  XOPER,XROOT
      DOUBLE PRECISION  AX
      AX = XOPER
      XROOT = DSIGN (DSQRT (ABS (AX)),AX)
      RETURN
      END
