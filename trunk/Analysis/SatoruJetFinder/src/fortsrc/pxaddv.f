CDECK  ID>, PXADDV.
      SUBROUTINE PXADDV (ISIZ,VEC1,VEC2,VECO)
*.*********************************************************
*. ------
*. PXADDV
*. ------
*. SOURCE:  J.W.Gary
*. Add two vectors of arbitrary length
*. Usage     :
*.
*.      INTEGER  ISIZ,NDIM
*.      PARAMETER  (NDIM=1.or.more)
*.      REAL  VEC1 (NDIM.or.more),
*.     +      VEC2 (NDIM.or.more),
*.     +      VECO (NDIM.or.more)
*.
*.      ISIZ = 1.to.NDIM
*.      CALL PXADDV (ISIZ,VEC1,VEC2,VECO)
*.
*. INPUT     : ISIZ    The length of the vectors
*. INPUT     : VEC1    The first vector
*. INPUT     : VEC2    The second vector
*. OUTPUT    : VECO    The vector sum of VEC1 and VEC2
*.                     (elements 1 to ISIZ)
*.
*.*********************************************************
      IMPLICIT NONE
      INTEGER  ISIZ,IX
      REAL  VEC1 (*),VEC2 (*),VECO (*)
      DOUBLE PRECISION  AX
      DO 120  IX = 1,ISIZ
          AX = VEC1 (IX) + VEC2 (IX)
          VECO (IX) = AX
 120  CONTINUE
      RETURN
      END
