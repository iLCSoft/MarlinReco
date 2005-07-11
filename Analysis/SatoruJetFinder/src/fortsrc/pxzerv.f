CDECK  ID>, PXZERV.
      SUBROUTINE PXZERV (ISZE,VEC)
*.*********************************************************
*. ------
*. PXZERV
*. ------
*. SOURCE: J.W.Gary
*. Zero a vector of arbitrary length
*. Usage     :
*.
*.      INTEGER  NDIM
*.      PARAMETER  (NDIM=1.or.more)
*.      REAL  VEC (NDIM)
*.      INTEGER  ISIZ
*.
*.      ISIZ = 1.to.NDIM
*.      CALL PXZERV (ISZE,VEC)
*.
*. INPUT     : ISIZ    The length of the vector to be zeroed
*. INPUT     : VEC     the vector to be zeroed
*.
*.*********************************************************
      IMPLICIT NONE
      INTEGER  ISZE,IX
      REAL  VEC (*)
      DO 120  IX = 1,ISZE
          VEC (IX) = 0.
 120  CONTINUE
      RETURN
      END
