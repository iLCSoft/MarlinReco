CDECK  ID>, PXOLAP.
******............................................................******
      SUBROUTINE PXOLAP(NJET,NTRAK,JETLIS,PJ,PP)
******............................................................******
*
*** Looks for particles assigned to more than 1 jet, and reassigns them
*** If more than a fraction OVLIM of a jet's energy is contained in
*** higher energy jets, that jet is neglected.
*** Particles assigned to the jet closest in angle (a la CDF, Snowmass).
*.
*. Modification Log.
*. 22-Apr-97: D. Chrisman - Dimension IJET with MXPROT instead of 30.
*.
      IMPLICIT NONE
      INTEGER MXTRAK, MXPROT
      PARAMETER (MXTRAK=200,MXPROT=100)
      INTEGER NJET, NTRAK
      LOGICAL JETLIS(MXPROT,MXTRAK)
      REAL PJ(4,MXPROT),PP(4,MXTRAK)
      INTEGER I,J,N,MU
      LOGICAL OVELAP
      REAL EOVER
      REAL OVLIM
      INTEGER ITERR, IJMIN, IJET(MXPROT), NJ
      REAL VEC1(3), VEC2(3), COST, THET, THMIN
      PARAMETER (OVLIM = 0.75)
*
      IF (NJET.LE.1) RETURN
*** Look for jets with large overlaps with higher energy jets.
      DO 100 I = 2,NJET
*** Find overlap energy between jets I and all higher energy jets.
       EOVER = 0.0
       DO 110 N = 1,NTRAK
         OVELAP = .FALSE.
         DO 120 J= 1,I-1
           IF (JETLIS(I,N).AND.JETLIS(J,N)) THEN
            OVELAP = .TRUE.
           ENDIF
120      CONTINUE
         IF (OVELAP) THEN
           EOVER = EOVER + PP(4,N)
         ENDIF
110     CONTINUE
*** Is the fraction of energy shared larger than OVLIM?
        IF (EOVER.GT.OVLIM*PJ(4,I)) THEN
*** De-assign all particles from Jet I
            DO 130 N = 1,NTRAK
              JETLIS(I,N) = .FALSE.
130         CONTINUE
         ENDIF
100   CONTINUE
*** Now there are no big overlaps, assign every particle in
*** more than 1 jet to the closet jet.
*** Any particles now in more than 1 jet are assigned to the CLOSET
*** jet (in angle).
      DO 140 I=1,NTRAK
         NJ=0
         DO 150 J=1, NJET
         IF(JETLIS(J,I)) THEN
            NJ=NJ+1
            IJET(NJ)=J
         ENDIF
150      CONTINUE
         IF (NJ .GT. 1) THEN
*** Particle in > 1 jet - calc angles...
            VEC1(1)=PP(1,I)
            VEC1(2)=PP(2,I)
            VEC1(3)=PP(3,I)
            THMIN=0.
            DO 160 J=1,NJ
               VEC2(1)=PJ(1,IJET(J))
               VEC2(2)=PJ(2,IJET(J))
               VEC2(3)=PJ(3,IJET(J))
               CALL PXANG3(VEC1,VEC2,COST,THET,ITERR)
               IF (J .EQ. 1) THEN
                  THMIN=THET
                  IJMIN=IJET(J)
               ELSEIF (THET .LT. THMIN) THEN
                  THMIN=THET
                  IJMIN=IJET(J)
               ENDIF
160         CONTINUE
*** Assign track to IJMIN
            DO 170 J=1,NJET
               JETLIS(J,I) = .FALSE.
170         CONTINUE
            JETLIS(IJMIN,I)=.TRUE.
         ENDIF
140   CONTINUE
*** Recompute PJ
      DO 200 I = 1,NJET
        DO 210 MU = 1,4
          PJ(MU,I) = 0.0
210     CONTINUE
        DO 220 N = 1,NTRAK
          IF( JETLIS(I,N) ) THEN
            DO 230 MU = 1,4
             PJ(MU,I) = PJ(MU,I) + PP(MU,N)
230         CONTINUE
          ENDIF
220     CONTINUE
200   CONTINUE
      RETURN
      END
