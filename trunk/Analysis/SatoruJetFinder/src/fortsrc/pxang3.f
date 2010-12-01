CDECK  ID>, PXANG3.
      SUBROUTINE PXANG3 (VEC1,VEC2,COST,THET,IERR)
*.*********************************************************
*. ------
*. PXANG3
*. ------
*. SOURCE: VECSUB (V. Blobel)
*. Calculate the angle beteen two 3-vectors
*. Usage     :
*.
*.      INTEGER  IERR
*.      REAL  VEC1 (3.or.more),
*.     +      VEC2 (3.or.more)
*.      REAL  COST,THET
*.
*.      CALL PXANG3 (VEC1,VEC2,COST,THET,IERR)
*.
*. INPUT     : VEC1    The first vector
*. INPUT     : VEC2    The second vector
*. OUTPUT    : COST    Cosine of the angle between the vectors
*. OUTPUT    : THET    The angle between the vectors (radians)
*. OUTPUT    : IERR    = 0 if all is OK ;   = -1 otherwise
*.
*.*********************************************************
      IMPLICIT NONE
      DOUBLE PRECISION AX,BX,CX,DX
      REAL  VEC1 (*),VEC2 (*)
      REAL  COST,THET
      INTEGER  IX,IERR
      IERR = 0
      AX = 0D0
      BX = 0D0
      CX = 0D0
      DO 120  IX = 1,3
          AX = AX + VEC1 (IX) * VEC1 (IX)
          BX = BX + VEC2 (IX) * VEC2 (IX)
          CX = CX + VEC1 (IX) * VEC2 (IX)
 120  CONTINUE
      DX = DSQRT (AX * BX)
      IF (DX.NE.0.0) THEN
          DX = CX / DX
      ELSE
          WRITE (6,FMT='('' PXANG3: Error, DX='',E12.4)') DX
          IERR = -1
          RETURN
      END IF
      IF (DABS (DX).GT.1.D0) DX = DSIGN (1.D0,DX)
      COST = DX
      THET = DACOS (DX)
      RETURN
      END
