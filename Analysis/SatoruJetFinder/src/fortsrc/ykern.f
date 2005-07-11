CDECK  ID>, YKERN.
      SUBROUTINE YKERN(IMODE,NT,ITKDM,PP,IERR)
C
C  JET FINDING ROUTINE A LA THE ORIGINAL JADE (E0) SCHEME, INCLUDING
C  THE E-,P- AND P0-VARIANTS. ALSO FEATURES THE NEW DURHAM (D)
C  AND GENEVA (G) SCHEMES.
C
C  INPUTS: IMODE (JET FINDING SCHEME); PP-ARRAY CONTAINING THE
C  FOUR-MOMENTA OF THE SELECTED PARTICLES, NT SPECIFYING HOW
C  MANY LOCATIONS IN PP ARE FILLED.
C
C  OUTPUTS:  YREC(I) CONTAINS VALUE OF Y WHERE EVENT FLIPS FROM
C  (I+1)-JET TO I-JET; PJET(K,I,J) CONTAINS THE JET AXES FOUR VECTOR
C  (K=1-4) OF JET NUMBER I WHEN THE EVENT IS CLASSIFIED AS J-JET.
C  LOCATION K=7 GIVES THE NUMBER OF PARTICLES ASSIGNED TO THIS JET.
C  NOTE THAT EACH EVENT IS ALWAYS FULLY RECONSTRUCTED DOWN TO
C  1-JET CONFIGURATION.
C  PROGRAM IS SELF-CONTAINED.
C
C  IMODE=1: E0 (=JADE) 2: P  3: P0  4: E  5: DURHAM (KT)  6: GENEVA
C
C  LAST MOD  :  21-Jun-98
C
C  Modification Log.
C  08-Oct-97 D. Chrisman, Increase NYCLMX from 250 to 500.
C  21-Jun-98 D. Chrisman, Want to be able to handle more than 10 jets.
C                         Introduce parameter NJETMX. Use this
C                         in declaration of PJET and YREC in order to remove
C                         hard coded numbers. Remove variable JCHECK and
C                         replace with NJETMX.
C
      IMPLICIT NONE
C     IMPLICIT NONE
      INTEGER NT,ITKDM,NYCLMX,NJETMX,NCALL,NPRINT,I,J,K,IERR,IMODEO
      INTEGER IMODE
      PARAMETER (NYCLMX = 500, NJETMX = 20)
      INTEGER NTO,NJETO,NJJ,IM,JM,NOLD,KK,II,HISTOR(2,NYCLMX)
      REAL PL(7,NYCLMX),Y(NYCLMX,NYCLMX),PINT(10,NYCLMX)
      REAL YREC(NJETMX),EVISO
      REAL PP(ITKDM,*),JADE,D,G,E,EVIS,PVIS,PJET(10,NJETMX,NJETMX),YMINI
      COMMON /YCL/YREC,PJET,HISTOR
      COMMON /YINT/ IMODEO,NTO,NJETO,PINT
      CHARACTER*7 CM
      SAVE NCALL,NPRINT,Y
      DATA NCALL,NPRINT /0,0/
C
C  JET RESOLUTION FUNCTIONS
C
      JADE(I,J) = 2.*PL(4,I)*PL(4,J)*MAX(0.,(1.-
     + (PL(1,I)*PL(1,J)+PL(2,I)*PL(2,J)+PL(3,I)*PL(3,J))/
     + (PL(6,I)*PL(6,J))))
      D(I,J) = 2.*MIN(PL(4,I)*PL(4,I),PL(4,J)*PL(4,J))*MAX(0.,(1.-
     + (PL(1,I)*PL(1,J)+PL(2,I)*PL(2,J)+PL(3,I)*PL(3,J))/
     + (PL(6,I)*PL(6,J))))
      E(I,J) = MAX(0.,(PL(4,I)+PL(4,J))**2-(PL(1,I)+PL(1,J))**2-
     + (PL(2,I)+PL(2,J))**2-(PL(3,I)+PL(3,J))**2)
      G(I,J) = 8.*PL(4,I)*PL(4,J)*MAX(0.,(1.-
     + (PL(1,I)*PL(1,J)+PL(2,I)*PL(2,J)+PL(3,I)*PL(3,J))/
     + (PL(6,I)*PL(6,J))))/(9.*(PL(4,I)+PL(4,J))**2)
C
C INITIALIZE
C
      IERR = -1
      IF(NCALL.LE.0) THEN
        IMODEO = 0
        NTO = 0
        NJETO = 0
      ENDIF
      DO 5001 I=1,NJETMX
        YREC(I) = 0.
 5001 CONTINUE
      NCALL = NCALL + 1
C
C PRINT JET SCHEME
C
      IF(IMODE.NE.IMODEO .AND. NPRINT.LE.7) THEN
        IF(IMODE.EQ.1) THEN
          CM = 'JADE E0'
        ELSEIF(IMODE.EQ.2) THEN
          CM = 'JADE P '
        ELSEIF(IMODE.EQ.3) THEN
          CM = 'JADE P0'
        ELSEIF(IMODE.EQ.4) THEN
          CM = 'JADE E '
        ELSEIF(IMODE.EQ.5) THEN
          CM = 'DURHAM '
        ELSEIF(IMODE.EQ.6) THEN
          CM = 'GENEVA '
        ELSE
          WRITE(6,281) IMODE
 281      FORMAT(/,' ### YKERN: IMODE =',I3,' INVALID; SET TO 1 ###')
          IMODE = 1
          CM = 'JADE E0'
        ENDIF
        PRINT 789, CM
 789    FORMAT(/,8X,54('#'),/,8X,
     +  '# YCLUS JET FINDER WITH ',A7,' RECOMBINATION SCHEME #',/,
     +  8X,54('#'),/)
        NPRINT = NPRINT + 1
      ENDIF
C
C CHECK INPUT PARAMETERS
C
      IF(NT.LE.1) THEN
        PRINT 701, NT
 701    FORMAT(/,' ### YKERN: NUMBER OF INPUT PARTICLES ',I4,' < 2;',
     +  ' NO JET RECONSTRUCTION DONE. ###')
        RETURN
      ENDIF
      IF(NT.GT.NYCLMX) THEN
        PRINT 700, NT, NYCLMX,NYCLMX
  700   FORMAT(/,' #### YKERN: NUMBER OF INPUT PARTICLES ',I4,' > ',I4,
     +  '; SET TO ',I4,/,12X,' INCREASE NYCLMX IF THIS',
     +  ' WARNING OCCURS MORE OFTEN. #### ')
        NT = NYCLMX
      ENDIF
C
C COPY INPUT VECTORS INTO INTERNAL MOMENTUM ARRAY (PL)
C  (POSITION  7 OF EACH VECTOR IN PL WILL BE OVERWRITTEN BY THIS ROUTINE
C   IN ORDER TO RECORD NUMBER OF INITIAL PARTICLES BELONGING TO THIS
C   POSITION [= 1. INITIALLY])
C
      EVIS = 0.
      PVIS = 0.
      DO 5002 I=1,NT
        PL(6,I)=SQRT(PP(1,I)**2 + PP(2,I)**2 + PP(3,I)**2)
        PL(1,I)=PP(1,I)
        PL(2,I)=PP(2,I)
        PL(3,I)=PP(3,I)
        IF(IMODE.EQ.2 .OR. IMODE.EQ.3) THEN
          PL(4,I) = PL(6,I)
        ELSE
          PL(4,I) = PP(4,I)
        ENDIF
        PL(7,I) = 1.
        DO 5003 K=1,7
          PINT(K,I) = PL(K,I)
 5003   CONTINUE
        EVIS=EVIS+PL(4,I)
        PVIS=PVIS+PL(6,I)
 5002 CONTINUE
      NJJ = NT
C
      IF(EVIS.LE.0. .OR. PVIS.GT.1.001*EVIS) THEN
        WRITE(6,283) EVIS,PVIS
 283    FORMAT(' #### YKERN: INCOMPATIBLE SUMS OF ENERGIES',
     +  ' AND/OR MOMENTA:',
     +  /,12X,' SUM(E) = ',F11.4,' SUM(P) = ',F11.4,/,12X,' CHECK',
     +  ' INPUT VECTORS AND/OR ARRAY DIMENSIONS. ####')
        RETURN
      ENDIF
C
      IF(NJJ.LE.NJETMX) THEN
        DO 5004 J=1,NJJ
          DO 5005 I=1,7
            PJET(I,J,NJJ) =  PL(I,J)
 5005     CONTINUE
 5004   CONTINUE
      ENDIF
C
C CALCULATE AND STORE PAIR MASSES
C
       DO 5006 I=1,NJJ-1
         DO 5007 J=I+1,NJJ
           IF(IMODE.EQ.1) THEN
             Y(I,J) = JADE(I,J)
           ELSEIF(IMODE.EQ.5) THEN
             Y(I,J) = D(I,J)
           ELSEIF(IMODE.EQ.6) THEN
             Y(I,J) = G(I,J)
           ELSE
             Y(I,J) = E(I,J)
           ENDIF
 5007    CONTINUE
 5006  CONTINUE
C
C FIND LOCAL MINIMUM OF PAIR MASSES
C
      IM = 0
      JM = 0
 1000 CONTINUE
      YMINI = 2.*EVIS**2
      DO 5008 I=1,NJJ-1
        DO 5009 J=I+1,NJJ
          IF(Y(I,J).LT.YMINI) THEN
           YMINI = Y(I,J)
           IM = I
           JM = J
          ENDIF
 5009   CONTINUE
 5008 CONTINUE
C
      IF(IM.EQ.0 .OR. JM.EQ.0) THEN
        WRITE(6,284)
  284   FORMAT(' #### YKERN: ERROR; NO MINIMUM FOUND IN Y-ARRAY! ####')
        RETURN
      ENDIF
C
C RECOMBINE PARTICLES IM AND JM, STORE AT POSITION IM
C
      PL(1,IM) = PL(1,IM) + PL(1,JM)
      PL(2,IM) = PL(2,IM) + PL(2,JM)
      PL(3,IM) = PL(3,IM) + PL(3,JM)
      PL(6,IM) = SQRT(PL(1,IM)**2 + PL(2,IM)**2
     +         + PL(3,IM)**2)
      IF(IMODE.EQ.2) THEN
        PL(4,IM) = PL(6,IM)
      ELSEIF(IMODE.EQ.3) THEN
        EVISO = EVIS
        EVIS = EVIS-PL(4,IM)-PL(4,JM)+PL(6,IM)
        PL(4,IM) = PL(6,IM)
      ELSE
        PL(4,IM) = PL(4,IM) + PL(4,JM)
      ENDIF
C                  KEEP TRACK OF # OF PARTICLES ASSIGNED TO NEW CLUSTER:
      PL(7,IM) = PL(7,IM) + PL(7,JM)
C                  MOVE LAST PARTICLE IN LIST (NJJ) TO EMPTY SLOT (JM)
      NOLD = 0
      IF(JM.NE.NJJ) THEN
        DO 5010 KK=1,7
          PL(KK,JM) = PL(KK,NJJ)
 5010   CONTINUE
        NOLD = NJJ
      ENDIF
      HISTOR(1,NJJ) = JM
      HISTOR(2,NJJ) = IM
      NJJ = NJJ - 1
C                  REMEMBER JET AXES AND VALUE OF YIJ OF FLIP
      IF(NJJ.LE.NJETMX) THEN
        IF(IMODE.EQ.3) THEN
          YREC(NJJ) = YMINI/EVISO**2
        ELSEIF(IMODE.EQ.6) THEN
          YREC(NJJ) = YMINI
        ELSE
          YREC(NJJ) = YMINI/EVIS**2
        ENDIF
        DO 5011 I=1,NJJ
          DO 5012 K=1,7
            PJET(K,I,NJJ) = PL(K,I)
 5012     CONTINUE
 5011   CONTINUE
      ENDIF
C
C END IF 1-JET CASE REACHED
C
      IF(NJJ.LE.1) GOTO 9000
C
C NOW CALCULATE RELEVANT NEW MASS-COMBINATIONS
C
      DO 5013 II=1,NJJ
        IF(II.NE.IM) THEN
         I = MIN(II,IM)
         J = MAX(II,IM)
         IF(IMODE.EQ.1) THEN
           Y(I,J) = JADE(I,J)
         ELSEIF(IMODE.EQ.5) THEN
           Y(I,J) = D(I,J)
         ELSEIF(IMODE.EQ.6) THEN
           Y(I,J) = G(I,J)
         ELSE
           Y(I,J) = E(I,J)
         ENDIF
        ENDIF
        IF(NOLD.NE.0) THEN
          I = MIN(II,JM)
          J = MAX(II,JM)
          Y(I,J) = Y(II,NOLD)
        ENDIF
 5013 CONTINUE
C
C  BACK TO START OF LOOP
C
      GOTO 1000
C
 9000 CONTINUE
C
      IMODEO = IMODE
      NTO = NT
      IERR = 0
CCCCC
C     WRITE(6,825) (J,(HISTOR(I,J),I=1,2),J=1,NT)
C825  FORMAT(/,' HISTORY:',/,250(I3,4X,2I4,/))
CCCCC
      RETURN
      END
