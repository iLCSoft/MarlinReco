CDECK  ID>, SYMRGJ.
**************************************************************************
      SUBROUTINE SYMRGJ
     +  (NFIN,IEXAM4,N0JET,NPART,IJOUT,QJET,IJFIN,QJFIN4,IERR)
*  Jet Merge Routine
*  Authour : Satoru Yamashita
*  Creat   : 1-Oct-97
*  Mod     : 1-Feb-98 S.Yamashita buffer size increased
*
*    Inputs:
*      NFIN   ;  requested Number of Jets
*      IEXAM4 ;  Jet merge method  1-8
*      N0JET  ;  Number of original Jets (>=NFIN)
*      NPART  ;  Number of particles
*      IJOUT  ;  Original Jet assignment for "particle"s
*      QJET   ;  Original Jet 6-momenta
*    Outputs:
*      IJFIN  ;  New Jet assignment of "particle"s for the Merged Jet
*      QJFIN4 ;  New merged Jet 6-momenta
*      IERR   ;  0=O.K   Others=error
      implicit none
      INTEGER MXJET
      PARAMETER (MXJET = 30)
      INTEGER I,J,N,IMINS,JMINS
      REAL    EXAMIN
      INTEGER NFIN,IEXAM4
      INTEGER N0JET,NPART
      INTEGER IJOUT(*)
      REAL    QJET(6,*)
      INTEGER IJFIN(*)
      REAL    QJFIN4(6,*),QTEMP(6,MXJET)
      INTEGER IERR

      REAL    DURHAM(MXJET,MXJET)
      REAL    E0JADE(MXJET,MXJET)
      REAL    EJADE(MXJET,MXJET)
      REAL    ANGLE(MXJET,MXJET)
      REAL    ECOS(MXJET,MXJET)
      REAL    RMULT(MXJET,MXJET)
      REAL    PL(6,MXJET)
      INTEGER MXEXAM
      PARAMETER (MXEXAM = 8)
      REAL    EXAM(MXEXAM,1000)
      INTEGER NNJET,IDJET(MXJET),IDTEMP(MXJET)

*            EXAM(1,N) = PL(4,I)*PL(4,J)
*            EXAM(2,N) = RMULT(I,J)
*            EXAM(3,N) = DURHAM(I,J)
*            EXAM(4,N) = E0JADE(I,J)
*            EXAM(5,N) = EJADE(I,J)
*            EXAM(6,N) = ANGLE(I,J)
*            EXAM(7,N) = ECOS(I,J)
*            EXAM(8,N) = PL(6,I)*PL(6,J)
*
CC OLD  1996  S.Y
CC            EXAM(1) = THETA
CC            EXAM(2) = - SQRT(EMIN2) * COSTH
CC            EXAM(3) = DURHAM
CC            EXAM(4) = E0JADE
CC            EXAM(5) = EJADE
*
*
** INITIALYSE
      DO I=1,NFIN
         DO J=1,6
            QJFIN4(J,I) = 0.0
         ENDDO
      ENDDO
      DO I=1,NPART
         IJFIN(I) = 0
      ENDDO
      IERR = 0

      IF(N0JET.LT.NFIN) THEN
         IERR = 1
         RETURN
      ELSE IF(N0JET.EQ.NFIN) THEN  ! Just Copy
         DO I=1,NFIN
            DO J=1,6
               QJFIN4(J,I) = QJET(J,I)
            ENDDO
         ENDDO
         DO I=1,NPART
            IJFIN(I) = IJOUT(I)
         ENDDO
         RETURN
      ENDIF

      IF(IEXAM4.LT.1.OR.IEXAM4.LT.MXEXAM) THEN
         IERR = -80
         RETURN
      ENDIF

      DO I=1,N0JET
         DO J=1,6
            PL(J,I) = QJET(J,I)
         ENDDO
         IDJET(I) = I
      ENDDO

      NNJET = N0JET
 100  CONTINUE
      IF(NNJET.LE.NFIN) THEN
         GOTO 777
      ENDIF

      N=0
      DO I=1,NNJET-1
         DO J=I+1,NNJET
            N=N+1

            E0JADE(I,J) =
     +           2.*PL(4,I)*PL(4,J)*MAX(0.,(1.-
     +           (PL(1,I)*PL(1,J)+PL(2,I)*PL(2,J)+PL(3,I)*PL(3,J))/
     +           (PL(6,I)*PL(6,J))))
            DURHAM(I,J) =
     +           2.*MIN(PL(4,I)*PL(4,I),PL(4,J)*PL(4,J))*MAX(0.,(1.-
     +           (PL(1,I)*PL(1,J)+PL(2,I)*PL(2,J)+PL(3,I)*PL(3,J))/
     +           (PL(6,I)*PL(6,J))))
            EJADE(I,J) =
     +           MAX(0.,(PL(4,I)+PL(4,J))**2-(PL(1,I)+PL(1,J))**2-
     +           (PL(2,I)+PL(2,J))**2-(PL(3,I)+PL(3,J))**2)
            ANGLE(I,J) =
     +           (PL(1,I)*PL(1,J)+PL(2,I)*PL(2,J)+PL(3,I)*PL(3,J))/
     +           (PL(6,I)*PL(6,J))
            ANGLE(I,J) = MIN(ANGLE(I,J), 1.)
            ANGLE(I,J) = MAX(ANGLE(I,J),-1.)
            ANGLE(I,J) = ACOS(ANGLE(I,J))
            ECOS(I,J)  = -MIN(PL(4,I),PL(4,J)) *
     +           (PL(1,I)*PL(1,J)+PL(2,I)*PL(2,J)+PL(3,I)*PL(3,J))/
     +           (PL(6,I)*PL(6,J))
            RMULT(I,J) = PL(4,I)*PL(4,J)-PL(1,I)*PL(1,J)
     +           -PL(2,I)*PL(2,J)-PL(3,I)*PL(3,J)

            EXAM(1,N) = PL(4,I)*PL(4,J)
            EXAM(2,N) = RMULT(I,J)
            EXAM(3,N) = DURHAM(I,J)
            EXAM(4,N) = E0JADE(I,J)
            EXAM(5,N) = EJADE(I,J)
            EXAM(6,N) = ANGLE(I,J)
            EXAM(7,N) = ECOS(I,J)
            EXAM(8,N) = PL(6,I)*PL(6,J)
         ENDDO
      ENDDO

      N=0
      IMINS=0
      JMINS=0
      EXAMIN = 1.E+10
      DO I=1,NNJET-1
         DO J=I+1,NNJET
            N=N+1
            IF(EXAM(IEXAM4,N).LT.EXAMIN) THEN
               IMINS=I
               JMINS=J
               EXAMIN=EXAM(IEXAM4,N)
            ENDIF
         ENDDO
      ENDDO

*     MERGE JETS
      DO I=1,NNJET
         DO J=1,6
            QTEMP(J,I) = PL(J,I)
         ENDDO
      ENDDO
      DO I=1,N0JET
         IDTEMP(I) = IDJET(I)
      ENDDO
      DO I=1,N0JET
         IDJET(I) = 0
      ENDDO

      N = 0
      DO I=1,NNJET
         IF(I.NE.IMINS.AND.I.NE.JMINS) THEN
            N = N + 1
            DO J=1,6
               PL(J,N) = QTEMP(J,I)
            ENDDO
            DO J=1,N0JET
               IF(IDTEMP(J).EQ.I) THEN
                  IDJET(J)   = N
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      NNJET = NNJET - 1
      DO J=1,6
         PL(J,NNJET) = QTEMP(J,IMINS)+QTEMP(J,JMINS)
      ENDDO
      DO J=1,N0JET
         IF(IDTEMP(J).EQ.IMINS.OR.IDTEMP(J).EQ.JMINS) THEN
            IDJET(J)   = NNJET
         ENDIF
      ENDDO
      GOTO 100

 777  CONTINUE
      DO I=1,NFIN
         DO J=1,6
            QJFIN4(J,I) = PL(J,I)
         ENDDO
      ENDDO
      DO I=1,NPART
         IJFIN(I) = IDJET(IJOUT(I))
      ENDDO

      END
