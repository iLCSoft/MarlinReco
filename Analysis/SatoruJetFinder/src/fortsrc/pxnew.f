CDECK  ID>, PXNEW.
******............................................................******
       LOGICAL FUNCTION PXNEW(TSTLIS,JETLIS,NTRAK,NJET)
******............................................................******
      IMPLICIT NONE
*
      INTEGER MXTRAK,MXPROT
      PARAMETER (MXTRAK=200,MXPROT=100)
       INTEGER NTRAK,NJET
       LOGICAL TSTLIS(MXTRAK),JETLIS(MXPROT,MXTRAK)
*** Checks to see if TSTLIS entries correspond to a jet already found
*** and entered in JETLIS
       INTEGER N, I
       LOGICAL MATCH
*
       PXNEW = .TRUE.
       DO 100 I = 1,NJET
          MATCH = .TRUE.
          DO 110 N = 1,NTRAK
            IF(TSTLIS(N).NEQV.JETLIS(I,N)) THEN
             MATCH = .FALSE.
             GO TO 100
            ENDIF
110       CONTINUE
          IF (MATCH) THEN
           PXNEW = .FALSE.
           RETURN
          ENDIF
100    CONTINUE
       RETURN
       END
