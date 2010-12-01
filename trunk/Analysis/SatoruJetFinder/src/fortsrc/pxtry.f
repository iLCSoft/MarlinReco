CDECK  ID>, PXTRY.
******............................................................******
       SUBROUTINE PXTRY(COSR,NTRAK,PU,PP,OAXIS,NAXIS,PNEW,NEWLIS,OK)
******............................................................******
*
      IMPLICIT NONE
      INTEGER MXTRAK
      PARAMETER (MXTRAK=200)
       INTEGER NTRAK
       REAL COSR,PU(3,MXTRAK),PP(4,MXTRAK),OAXIS(3)
       LOGICAL OK
       LOGICAL NEWLIS(MXTRAK)
       REAL NAXIS(3),PNEW(4)
*** Finds all particles in cone of size COSR about OAXIS direction.
*** Calculates 4-momentum sum of all particles in cone (PNEW) , and
*** returns this as new jet axis NAXIS (Both unit Vectors)
       INTEGER N,MU
       REAL COSVAL,NORMSQ,NORM
*
       OK = .FALSE.
       DO 100 MU=1,4
          PNEW(MU)=0.0
100    CONTINUE
       DO 110 N=1,NTRAK
          COSVAL=0.0
          DO 120 MU=1,3
             COSVAL=COSVAL+OAXIS(MU)*PU(MU,N)
120       CONTINUE
          IF (COSVAL.GE.COSR)THEN
             NEWLIS(N) = .TRUE.
             OK = .TRUE.
             DO 130 MU=1,4
                PNEW(MU) = PNEW(MU) + PP(MU,N)
130          CONTINUE
          ELSE
             NEWLIS(N)=.FALSE.
          ENDIF
110   CONTINUE
*** If there are particles in the cone, calc new jet axis
       IF (OK) THEN
          NORMSQ = 0.0
          DO 140 MU = 1,3
            NORMSQ = NORMSQ + PNEW(MU)**2
140       CONTINUE
          NORM = SQRT(NORMSQ)
          DO 150 MU=1,3
            NAXIS(MU) = PNEW(MU)/NORM
150       CONTINUE
       ENDIF
       RETURN
       END
