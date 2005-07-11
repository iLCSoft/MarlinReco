CDECK  ID>, PXUVEC.
*
******............................................................******
       SUBROUTINE PXUVEC(NTRAK,PP,PU,IERR)
******............................................................******
*
*** Routine to calculate unit vectors PU of all particles PP
      IMPLICIT NONE
      INTEGER MXTRAK
      PARAMETER (MXTRAK=200)
      INTEGER NTRAK, IERR
      REAL PP(4,MXTRAK)
      REAL PU(3,MXTRAK)
      INTEGER N,MU
      REAL MAG
       DO 100 N=1,NTRAK
          MAG=0.0
          DO 110 MU=1,3
             MAG=MAG+PP(MU,N)**2
110       CONTINUE
          MAG=SQRT(MAG)
          IF (MAG.EQ.0.0) THEN
             WRITE(*,*)' PXCONE: An input particle has zero mod(p)'
             IERR=-1
             RETURN
          ENDIF
          DO 120 MU=1,3
           PU(MU,N)=PP(MU,N)/MAG
120       CONTINUE
100    CONTINUE
       RETURN
       END
