CDECK  ID>, PXNORV.
      SUBROUTINE PXNORV (ISIZ,VEC,VNOR,IERR)
*.*********************************************************
*. ------
*. PXNORV
*. ------
*. SOURCE: J.W.Gary
*. Normalize a vector of arbitrary length
*. Usage     :
*.
*.      INTEGER  NDIM
*.      PARAMETER  (NDIM=1.or.more)
*.      REAL  VEC (NDIM.or.more)
*.      REAL  VNOR
*.      INTEGER  IERR,ISIZ
*.
*.      ISIZ = 1.to.NDIM
*.      CALL PXNORV (ISIZ,VEC,VNOR,IERR)
*.
*. INPUT     : ISIZ    The length of the vector
*. INPUT     : VEC     The vector
*. OUTPUT    : VNOR    The normalized vector
*. OUTPUT    : IERR    = 0 if all is OK ;   = -1 otherwise
*.
*.*********************************************************
      IMPLICIT NONE
      INTEGER  ISIZ,IX,IERR
      REAL  VEC (*),VNOR (*)
      DOUBLE PRECISION  AX,BX
      IERR = 0
      AX = 0.0
      DO 120  IX = 1,ISIZ
          AX = AX + VEC (IX) * VEC (IX)
 120  CONTINUE
      BX = DSQRT (AX)
      IF (BX.NE.0.0) THEN
          BX = 1.0 / BX
      ELSE
          IERR = -1
          RETURN
      END IF
      DO 140 IX = 1,ISIZ
          VNOR (IX) = BX * VEC (IX)
 140  CONTINUE
      RETURN
      END
