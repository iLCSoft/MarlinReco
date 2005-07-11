CDECK  ID>, PXCAMJ.
      SUBROUTINE PXCAMJ(ITKDM,NT,PT,YCUT,NJ,IJ,PJ,IERR)
C
C---CAMBRIDGE JET CLUSTERING ALGORITHM
C   BASED ON YCLUS BY S BETHKE
C   REF: YU L DOKSHITZER, G D LEDER, S MORETTI, B R WEBBER
C   CAVENDISH-HEP-97/06 (JUNE 1997)
C   07/07/97 FIRST RELEASED BY BRW
C   23/08/97 COMMENTS REVISED BY BRW
C
C   INPUT:
C    ITKDM = 1st dimension of array PT, ITKDM >= 4 required
C       NT = NUMBER OF TRACKS
C   PT(,I) = 4-MOMENTUM OF TRACK I (I=1,NT)
C     YCUT = (DURHAM) JET RESOLUTION
C
C  OUTPUT:
C       NJ = NUMBER OF JETS
C    IJ(I) = J IF TRACK I BELONGS TO JET J (I=1,NT)
C  PJ(4,J) = 4-MOMENTUM OF JET J (J=1,NJ)
C
C      NB:   CLUSTERING SEQUENCE DEPENDS ON VALUE OF YCUT
C
C  Modifications:
C  23.09.97, STK: Variable 1st dim for PT, REAL call args, handle
C                 errors with error flag IERR
C  24.09.97, STK: Improved error handling, sorted declarations,
C                 added welcome message
C  29.09.97, SB : Calculate EVIS from PT array; introduce to PX library
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ITKDM,NT,NJ,IJ(*),IERR
      REAL PT(ITKDM,*),EVIS,YCUT,PJ(4,*)
      INTEGER NTRK,NV
      PARAMETER( NTRK=300, NV=NTRK*(NTRK-2)+NTRK-(NTRK-2)*(NTRK-1)/2 )
      LOGICAL IP(NTRK),LCALL
      INTEGER I,II,J,K,L,IMINI,JMINI,IAD,JJ(NTRK)
      DOUBLE PRECISION PP(5),PL(5,NTRK),V(NV),PM,VMINI,YSCA
      SAVE LCALL
      DATA LCALL / .FALSE. /
C  Welcome message:
      IF( .NOT.LCALL ) THEN
        PRINT *, ' '
        PRINT *, 'Cambridge jet finding algorithm, please refer to:'
        PRINT *, 'Yu.L. Dokshitzer, G.D. Leder, S. Moretti, B.R. Webber'
        PRINT *, 'CAVENDISH-HEP-97/06 (June 1997)'
        PRINT *, ' '
        LCALL= .TRUE.
      ENDIF
C---WARNINGS
      IERR = 0
      EVIS = 0.
      IF( NT.GT.NTRK ) THEN
        WRITE(*,'(''CAMJET: More than '',I3,'' input particles: '',I5)')
     &        NTRK,NT
        IERR= 1
        RETURN
      ELSEIF( NT.LT.2 ) THEN
        WRITE(*,'(''CAMJET: Less than 2 input particles: '',I5)') NT
        IERR= 2
        RETURN
      ENDIF
C---COPY MOMENTA INTO PL-ARRAY
      DO I=1,NT
        IP(I)=.TRUE.
        IJ(I)=I
        DO II=1,4
          PL(II,I)= DBLE(PT(II,I))
        ENDDO
        PM=PL(1,I)**2+PL(2,I)**2+PL(3,I)**2
        EVIS = EVIS + PT(4,I)
        IF (PM.GT.0D0) THEN
          PL(5,I)=1D0/SQRT(PM)
        ELSE
          PL(5,I)=1D0
        ENDIF
      ENDDO
      YSCA= DBLE(YCUT)*DBLE(EVIS)**2
C---FILL V-ARRAY: V(I,J) IS V(NT*(I-1)+J-I(I+1)/2)
      IAD = 0
      DO I=1,NT-1
        DO II=1,5
          PP(II)=PL(II,I)
        ENDDO
        DO J=I+1,NT
          IAD = IAD + 1
          V(IAD) = 2D0*(1D0-(PP(1)*PL(1,J) +PP(2)*PL(2,J)
     &                      +PP(3)*PL(3,J))*PP(5)*PL(5,J))
        ENDDO
      ENDDO
      NJ=NT
C---START MAIN LOOP.  FIRST LOOK FOR MINIMUM V
    1 VMINI = 1D10
      IMINI = 0
      IAD = 0
      DO I=1,NT-1
        IF (IP(I)) THEN
          DO J=I+1,NT
            IAD = IAD + 1
            IF (IP(J).AND.V(IAD).LT.VMINI) THEN
              VMINI = V(IAD)
              IMINI = I
              JMINI = J
            ENDIF
          ENDDO
        ELSE
          IAD=IAD+NT-I
        ENDIF
      ENDDO
C---END OF CLUSTER SEARCH FOR VMINI
      IF (IMINI.NE.0) THEN
C---NOT FINISHED YET
        IF (VMINI*MIN(PL(4,IMINI),PL(4,JMINI))**2.GE.YSCA) THEN
C---SOFT FREEZING HERE
          IF (PL(4,IMINI).LT.PL(4,JMINI)) THEN
            IP(IMINI)=.FALSE.
          ELSE
            IP(JMINI)=.FALSE.
          ENDIF
        ELSE
C---COMBINE PARTICLES IMINI AND JMINI
          DO II=1,4
            PL(II,IMINI)=PL(II,IMINI)+PL(II,JMINI)
          ENDDO
          PM=PL(1,IMINI)**2+PL(2,IMINI)**2+PL(3,IMINI)**2
          IF (PM.GT.0D0) THEN
            PL(5,IMINI)=1D0/SQRT(PM)
          ELSE
            PL(5,IMINI)=1D0
          ENDIF
C---FLAG PARTICLE JMINI AS COMBINED
          IP(JMINI)=.FALSE.
          IJ(JMINI)=IMINI
          NJ=NJ-1
C---CALCULATE RELEVANT NEW V VALUES
          DO I=1,NT
            IF (I.NE.IMINI) THEN
              IF (IJ(I).EQ.JMINI) IJ(I)=IMINI
              IF (IP(I)) THEN
                K = MIN(I,IMINI)
                L = MAX(I,IMINI)
                IAD = NT*(K-1) + L - (K*(K+1))/2
                V(IAD) = 2D0*(1D0-(PL(1,K)*PL(1,L) +PL(2,K)*PL(2,L)
     &                            +PL(3,K)*PL(3,L))*PL(5,K)*PL(5,L))
              ENDIF
            ENDIF
          ENDDO
        ENDIF
C---BACK TO START OF LOOP
        GO TO 1
      ELSE
C---FINISHED: CONSTRUCT JETS
        J=0
        DO I=1,NT
          IF (IJ(I).EQ.I) THEN
            J=J+1
            JJ(I)=J
            DO II=1,4
              PJ(II,J)= SNGL(PL(II,I))
            ENDDO
          ENDIF
        ENDDO
        DO I=1,NT
          IJ(I)=JJ(IJ(I))
        ENDDO
      ENDIF
      END
