CDECK  ID>, SYJETF2.
************************************************************************
      SUBROUTINE SYJETF2(IEXAM,N0JET,NPART,QIN,IJ,FRAC,IJOUT,QJET)
*  Particle Re-association Routine (Mode=2)
*  Jet core is always defined by the original jet
*
*  Authour : Satoru Yamashita
*  Creat   : 1-Oct-97
*  Mod     : 1-Feb-98 S.Yamashita buffer size increased
*
*    Inputs:
*      IEXAM  ;  Jet-particle Re-association method  1-6
*            EXAM(1) = THETA
*            EXAM(2) = PART(4)*PJSTMP(4,J)-PART(1)*PJSTMP(1,J)
*     +                -PART(2)*PJSTMP(2,J)-PART(3)*PJSTMP(3,J)
*            EXAM(3) = DURHAM
*            EXAM(4) = E0JADE
*            EXAM(5) = Geneva Scheme
*            EXAM(6) = EJADE   ! it makes too bad cross talk
*      N0JET  ;  Number of Jets
*      NPART  ;  Number of particles
*      QIN    ;  particles 5-momenta
*      IJ     ;  Original Jet assignment for "particle"s
*      FRAC   ;  Energy Fraction for the core
*    Outputs:
*      IJOUT  ;  New Jet assignment for "particle"s
*      QJET   ;  New Jet 6-momenta
*
      implicit none
      INTEGER IEXAM,N,I,J,K,II
      REAL    EMIN2,RINV2
      INTEGER N0JET,NPART
      REAL    QIN(5,*)
      INTEGER IJ(*)
      REAL    FRAC
      INTEGER IJOUT(*)
      REAL    QJET(6,*)
      REAL    E0JADE,EJADE,DURHAM
* ----
      INTEGER MXJET
      PARAMETER (MXJET = 30)
      INTEGER NJ(MXJET),III(500),ISORTJ(500,MXJET),
     +                 JETASS(500),ICUTJ(MXJET)
      REAL    PJSTMP(6,MXJET)
      INTEGER ISERI(500)
      REAL    PJ(5,500,MXJET)
      INTEGER ISERIJ(50,MXJET)
      REAL    Q(5,500),QSAVE(5,500)
      REAL    RSUM,RSUM0
      INTEGER NRES,IASTMP(500)
      REAL    PART(6)
      REAL    A,COSTH,THETA
      INTEGER MXEXAM
      PARAMETER (MXEXAM = 6)
      REAL    EXAM(MXEXAM)
      REAL    EXAMIN
      INTEGER IEXAMN

** INPUTS
      N=NPART
      DO I=1,N
         DO J=1,5
            Q(J,I)=QIN(J,I)
            QSAVE(J,I)=QIN(J,I)
         ENDDO
         ISERI(I)=I
      ENDDO
** INITIALYZE NUMBER OF PARTICLES IN JET
      DO I=1,MXJET
         NJ(I)= 0
      ENDDO
** PARTICLES ARRANGED IN JET
      DO I=1,N
*--- added for cone jet
         IF(IJ(I).GT.0) THEN
*--- add end
            NJ(IJ(I))=NJ(IJ(I))+1
            DO J=1,5
               PJ(J,NJ(IJ(I)),IJ(I)) = Q(J,I)
            ENDDO
            ISERIJ(NJ(IJ(I)),IJ(I)) = ISERI(I)
         ENDIF
      ENDDO

C      WRITE(6,*)NJ(1),NJ(2),NJ(3),NJ(4)
** LEADING PARTICLES SEARCH IN EACH JETS (ORIGINAL NJET) AND RE-DEFINE JETS
      DO I=1,N0JET
         DO J=1,6
            PJSTMP(J,I) = 0.
         ENDDO
         RSUM=0.
         ICUTJ(I) = 0
         DO J=1,NJ(I)
            DO K=1,5
               Q(K,J)=PJ(K,J,I)
            ENDDO
            RSUM = RSUM + Q(4,J)
         ENDDO
         CALL SYISRT(NJ(I),Q,III)
         DO J=1,NJ(I)
            ISORTJ(J,I) = III(J)
         ENDDO
         RSUM0 = 0.
         DO J=1,NJ(I)
            RSUM0 = RSUM0 + Q(4,III(J))
*       --- RE-CALCULATE JET MOMENTUM
            DO K=1,4
               PJSTMP(K,I) = PJSTMP(K,I) + Q(K,III(J))
            ENDDO
            PJSTMP(6,I)=
     +       SQRT(PJSTMP(1,I)**2+PJSTMP(2,I)**2+PJSTMP(3,I)**2)
            PJSTMP(5,I)=SQRT(MAX(0.,PJSTMP(4,I)**2-PJSTMP(6,I)**2))

            ICUTJ(I) = J
*       --- CHECK FRACTION OF ENERGY
            IF(RSUM0.GT.RSUM*FRAC) THEN
               GOTO 7
            ENDIF
         ENDDO
 7       CONTINUE
      ENDDO

C---re-calculate reference jet momentum
      DO I=1,N0JET
         DO J=1,6
            PJSTMP(J,I) = 0.
         ENDDO
         DO J=1,NJ(I)
            DO K=1,4
               PJSTMP(K,I) = PJSTMP(K,I) + PJ(K,J,I)
            ENDDO
         ENDDO
         PJSTMP(6,I)=
     +        SQRT(PJSTMP(1,I)**2+PJSTMP(2,I)**2+PJSTMP(3,I)**2)
         PJSTMP(5,I)=SQRT(MAX(0.,PJSTMP(4,I)**2-PJSTMP(6,I)**2))
      ENDDO

** JET ASSIGN CLEAR
      DO I=1,N
         JETASS(I)=0
      ENDDO
** JET ASSIGN AGAIN FOR LEADING PARTICLES
      DO I=1,N0JET
         DO J=1,ICUTJ(I)
            II=ISORTJ(J,I)
            II=ISERIJ(II,I)
            JETASS(II) = I
         ENDDO
      ENDDO
** PUT REMAINING PARTICLES INTO TEMP BUFFER
      NRES = 0
      DO I=1,N
         IF(JETASS(I).EQ.0) THEN
            NRES = NRES + 1
            IASTMP(NRES) = ISERI(I)
            DO J=1,5
               Q(J,NRES) = QSAVE(J,I)
            ENDDO
         ENDIF
      ENDDO

** SORT REMAINING PARTICLES IN ENERGY ORDERING
      IF(NRES.EQ.0) THEN
C         WRITE(*,*)' NRES 0'
         GOTO 99
      ENDIF
      CALL SYISRT(NRES,Q,III)
 99   CONTINUE
** ASSIGNMENT OF THE PARTICLES INTO JETS
      DO I=1,NRES
         II=III(I)      ! ID IN RES-BUFFER
         II=IASTMP(II)  ! SERIAL NUMBER
         DO J=1,5
            PART(J)=QSAVE(J,II)
         ENDDO
         PART(6)=SQRT(PART(1)**2+PART(2)**2+PART(3)**2)

         IEXAMN = 0
         EXAMIN = 1.E+10
         DO J=1,N0JET
            A = PJSTMP(1,J)*PART(1)+
     +          PJSTMP(2,J)*PART(2)+PJSTMP(3,J)*PART(3)
            A = A/PART(6)/PJSTMP(6,J)
            A = MIN(A, 1.)
            COSTH = MAX(A,-1.)   ! COSTH
            THETA = ACOS(COSTH)    ! THERA (ANGLE BETWEEN JET & PART)  (1)
*--
            EMIN2 = MIN(PART(4)**2,PJSTMP(4,J)**2)
            RINV2  = (PART(4)+PJSTMP(4,J))**2-(PART(1)+PJSTMP(1,J))**2
     +      -(PART(2)+PJSTMP(2,J))**2-(PART(3)+PJSTMP(3,J))**2
            EJADE = MAX(0.,RINV2)
            DURHAM= 2.*EMIN2*MAX(0.,1.-COSTH)
            E0JADE= 2.*PART(4)*PJSTMP(4,J)*MAX(0.,1.-COSTH)
            EXAM(1) = THETA
            EXAM(2) = PART(4)*PJSTMP(4,J)-PART(1)*PJSTMP(1,J)
     +                -PART(2)*PJSTMP(2,J)-PART(3)*PJSTMP(3,J)
            EXAM(3) = DURHAM
            EXAM(4) = E0JADE
            EXAM(5) = 8.*PART(4)*PJSTMP(4,J)*MAX(0.,(1.-
     +    (PART(1)*PJSTMP(1,J)+PART(2)*PJSTMP(2,J)+PART(3)*PJSTMP(3,J))/
     +    (PART(6)*PJSTMP(6,J))))/(9.*(PART(4)+PJSTMP(4,J))**2) ! geneve
            EXAM(6) = EJADE   ! it makes too bad cross talk
C            EXAM(2) = - SQRT(EMIN2) * COSTH

            IF(EXAM(IEXAM).LT.EXAMIN) THEN
               EXAMIN = EXAM(IEXAM)
               IEXAMN = J
            ENDIF

         ENDDO

*        NEW ASSIGN AND NEW JETS
         JETASS(II) = IEXAMN
C==== MODE=2 ; Do NOT recalculate Jets (always use original jets for reference)
C         DO K=1,4
C            PJSTMP(K,IEXAMN) =  PJSTMP(K,IEXAMN) + PART(K)
C         ENDDO
C         PJSTMP(6,IEXAMN) = PJSTMP(1,IEXAMN)**2+
C     +                      PJSTMP(2,IEXAMN)**2+PJSTMP(3,IEXAMN)**2
C         PJSTMP(5,IEXAMN) = PJSTMP(4,IEXAMN)**2-PJSTMP(6,IEXAMN)
C         PJSTMP(6,IEXAMN) = SQRT(MAX(0.,PJSTMP(6,IEXAMN)))
C         PJSTMP(5,IEXAMN) = SQRT(MAX(0.,PJSTMP(5,IEXAMN)))
      ENDDO

      DO I=1,NPART
         IJOUT(I) = JETASS(I)
      ENDDO
      DO I=1,N0JET
         DO J=1,6
            QJET(J,I)=PJSTMP(J,I)
         ENDDO
      ENDDO
C---re-calculate reference jet momentum
      DO I=1,N0JET
         DO J=1,6
            QJET(J,I) = 0.
         ENDDO
      ENDDO
      DO I=1,NPART
         IF(IJOUT(I).GT.0.AND.IJOUT(I).LE.N0JET) THEN
            DO J=1,4
               QJET(J,IJOUT(I)) = QJET(J,IJOUT(I)) + QIN(J,I)
            ENDDO
         ENDIF
      ENDDO
      DO I=1,N0JET
         QJET(6,I)=
     +        SQRT(QJET(1,I)**2+QJET(2,I)**2+QJET(3,I)**2)
         QJET(5,I)=SQRT(MAX(0.,QJET(4,I)**2-QJET(6,I)**2))
      ENDDO

      END
