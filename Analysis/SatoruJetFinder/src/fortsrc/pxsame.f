CDECK  ID>, PXSAME.
******............................................................******
       LOGICAL FUNCTION PXSAME(LIST1,LIST2,N)
******............................................................******
      IMPLICIT NONE
*
       LOGICAL LIST1(*),LIST2(*)
       INTEGER N
*** Returns T if the first N elements of LIST1 are the same as the
*** first N elements of LIST2.
       INTEGER I
*
       PXSAME = .TRUE.
       DO 100 I = 1,N
        IF ( LIST1(I).NEQV.LIST2(I) ) THEN
          PXSAME = .FALSE.
          RETURN
        ENDIF
100    CONTINUE
       RETURN
       END
