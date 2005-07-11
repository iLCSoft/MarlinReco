CDECK  ID>, CYJET.
      SUBROUTINE CYJET(NJET,YL,YH,IERR)
*.-----------------------------------------------------------------------
*.
*.    CYJET: Return the largest pair of y-flip values for a required n-jet
*.           configuration.
*.           INPUT:   NJET (integer) Number of jets required
*.           OUTPUT:  YL   (real   ) Lower boundary of y interval
*.                    YH   (real   ) Upper boundary of y interval
*.                    IERR (integer) Error flag, 0=OK.
*.
*.  CREATED :  11-12-1997, STAN BENTVELSEN
*.  LAST MOD:
*.-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NMXY , NMXP
      PARAMETER (NMXY = 300)
      PARAMETER (NMXP = 30)
      DOUBLE PRECISION YTRANS(NMXY)
      DOUBLE PRECISION PCMJ(NMXP,4,NMXP)
      INTEGER          NTRANS(NMXY), NJITER, NJMAX, NTRACK
      INTEGER          ICMJ(NMXP,NMXY)
      COMMON / CKCOM / YTRANS, PCMJ, NTRANS, ICMJ, NJITER, NJMAX, NTRACK
      REAL    YL,YH
      INTEGER NJET,I,IERR
      INTEGER NPRINT
      DATA    NPRINT / 0 /
      SAVE    NPRINT

      IERR = 0
      YL = -1.
      YH = -1.
      IF(NJET.GT.NJMAX) THEN
         WRITE(*,*) '#################################################'
         WRITE(*,*) '## CYJET: CAMBRIDGE JET FINDER RESOLVED FOR    ##'
         WRITE(*,'(A,I2,A)') ' ##    NUMBER OF JETS BETWEEN 1 AND '
     +        ,NJMAX,' ONLY     ##'
         WRITE(*,*) '#################################################'
         IERR = 2
         RETURN
      ENDIF
      IF(NJET.LT.1) THEN
         WRITE(*,*) 'ERROR IN NJET: ',NJET
         IERR = 3
         RETURN
      ENDIF

      DO I=1,NJITER
         IF(NTRANS(I).EQ.NJET) THEN
            YH = 1.*YTRANS(I)
            IF(I.LT.NJITER) YL = 1.*YTRANS(I+1)
            GOTO 889
         ENDIF
      ENDDO

      IF(NPRINT.LT.10) THEN
         NPRINT = NPRINT + 1
         WRITE(*,*) '#################################################'
         WRITE(*,*) '## CYJET: CAMBRIDGE JET FINDER CANNOT RESOLVE  ##'
         WRITE(*,'(A,I2,A)') ' ##         THIS CONFIGURATION TO A '
     +        ,NJET,' JET      ##'
         WRITE(*,*) '#################################################'
      ENDIF

 889  CONTINUE
      RETURN
      END
