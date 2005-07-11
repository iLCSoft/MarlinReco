CDECK  ID>, CASSO.
      SUBROUTINE CASSO(NJET,PNJ,BL,IERR)
*.-----------------------------------------------------------------------
*.
*.    CASSO: Return the jet four-momenta and jet-particle association for a
*.           required number of jets.
*.           (corresponding to largest ycut values)
*.           INPUT:   NJET (integer) Number of jets required
*.           OUTPUT:  PNJ(4,*)(real) Array of jets 4-vectors
*.                    BL(*)(integer) Particle i belongs to jet BL(i)
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
      INTEGER NJET, IERR
      REAL PNJ(4,*)
      INTEGER I,J,II
      INTEGER BL(*)
      INTEGER NPRINT
      DATA    NPRINT / 0 /
      SAVE    NPRINT

      IERR = 0
      IF(NJET.GT.NJMAX.OR.NJET.LE.0) THEN
         WRITE(*,*) '#################################################'
         WRITE(*,*) '## CASSO: CAMBRIDGE JET FINDER RESOLVED FOR    ##'
         WRITE(*,'(A,I2,A)') ' ##    NUMBER OF JETS BETWEEN 1 AND '
     +        ,NJMAX,' ONLY     ##'
         WRITE(*,*) '#################################################'
         DO I=1,NJMAX
            DO II=1,4
               PNJ(II,I) = 0
            ENDDO
         ENDDO
         IERR = 2
         RETURN
      ENDIF
      DO I=1,NJITER
         IF(NTRANS(I).EQ.NJET) THEN
            DO J=1,NJET
               DO II=1,4
                  PNJ(II,J) = SNGL(PCMJ(I,II,J))
               ENDDO
            ENDDO
            DO J=1,NTRACK
               BL(J) = ICMJ(I,J)
            ENDDO
            GOTO 889
         ENDIF
      ENDDO

      IF(NPRINT.LT.10) THEN
         NPRINT = NPRINT + 1
         WRITE(*,*) '#################################################'
         WRITE(*,*) '## CASSO: CAMBRIDGE JET FINDER CANNOT RESOLVE  ##'
         WRITE(*,'(A,I2,A)') ' ##         THIS CONFIGURATION TO A '
     +        ,NJET,' JET      ##'
         WRITE(*,*) '#################################################'
      ENDIF

 889  CONTINUE
      RETURN
      END
