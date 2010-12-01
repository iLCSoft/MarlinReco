CDECK  ID>, CKERN.
      SUBROUTINE CKERN(ITKDM,NT,PT,NJREQ,IERR)
*.-----------------------------------------------------------------------
*.
*.    CKERN: Resolves the clustering, and stores information
*.           in common block CKCOM.
*.           INPUT: ITKDM      (integer) First dimension of PT array
*.                                 (>=4 required)
*.                  NT         (integer) Number of tracks
*.                  PT(ITKDM,*)(real)    Four-momenta of tracks
*.                  NJREQ      (integer) Clustering limit. Note that the CPU
*.                                       consumption is roughly proportional to
*.                                       calling JADE/Durham algorithm NJREQ
*.                                       times!!
*.                                       1) If interested in ALL possible N-jet
*.                                          configurations, set NJREQ equal to
*.                                          NT. Note the remark on CPU time!
*.                                       2) If interested in first N-jet
*.                                          configuration (with largest
*.                                          ycut values), set NJREQ equal
*.                                          to number of jets N-jet.
*.           OUTPUT: IERR      (integer) Error flag, 0=OK.
*.           CALLS:  CKSORD
*.
*.
*.    Once CKERN is called, all information on the event clustering,
*.    final state jets and particle association can be accessed
*.    without additional CPU time using the utility routines.
*.
*.---CAMBRIDGE JET CLUSTERING ALGORITHM
*.   BASED ON YCLUS BY S BETHKE
*.   REF: YU L DOKSHITZER, G D LEDER, S MORETTI, B R WEBBER
*.   CAVENDISH-HEP-97/06 (JUNE 1997)
*.   07/07/97 FIRST RELEASED BY BRW
*.   23/08/97 COMMENTS REVISED BY BRW
*.   23/09/97 IMPLEMENT FOR PX LIBRARY BY STK AND SB
*.
*.
*.  CREATED :  11-12-1997, STAN BENTVELSEN
*.  LAST MOD:  11-02-1998, Z. Troscyani, Add factor epsilon
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
      INTEGER ITKDM,NT,NJ,NJREQ,IERR
      REAL PT(ITKDM,*)
      INTEGER NTRK,NV
      PARAMETER( NTRK=300, NV=NTRK*(NTRK-2)+NTRK-(NTRK-2)*(NTRK-1)/2 )
      LOGICAL IP(NTRK),LCALL
      INTEGER I,II,J,K,L,IMINI,JMINI,IAD,IJ(NTRK),JJ(NTRK)
      DOUBLE PRECISION PP(5),PL(5,NTRK),V(NV),PM,VMINI,YSCA,EVIS
C
      DOUBLE PRECISION EPSILON
      DOUBLE PRECISION YCUR(NTRK)
      INTEGER          ICUR, KIX(NTRK)
C

      SAVE LCALL
      DATA LCALL / .FALSE. /
      DATA EPSILON / 1D-12 /
      SAVE EPSILON

C  Welcome message:
      IF( .NOT.LCALL ) THEN
c        PRINT *, ' '
c        PRINT *, 'Cambridge jet finding algorithm, please refer to:'
c        PRINT *, 'Yu.L. Dokshitzer, G.D. Leder, S. Moretti, B.R. Webber'
c        PRINT *, 'CAVENDISH-HEP-97/06 (June 1997)'
        PRINT *, ' '
        PRINT *, ' #################################################'
        PRINT *, ' ###     C K E R N    V E R S I O N    104     ###'
c        PRINT *, ' ###  S. Bentvelsen and I. Meyer, EPXXX, 1998  ###'
        PRINT *, ' #################################################'
        WRITE(*,'(A,I2,A)')
     +      '  ##### CAMBRIDGE JET CLUSTERING 1 - ',NJREQ,' JETS  #####'
        PRINT *, ' #################################################'
        LCALL= .TRUE.
      ENDIF
C---WARNINGS
      IERR = 0
      EVIS = 0D0
      ICUR = 0
      NTRACK = NT
      NJREQ = MIN(NJREQ,NT)
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
C
C..   RESET ALL VALUES
C
      NJITER = 0
      DO I=1,NMXY
         YTRANS(I) = 0D0
         NTRANS(I) = 0
         DO J=1,NMXP
            ICMJ(J,I)   = 0
         ENDDO
      ENDDO

      DO I=1,NT
         EVIS = EVIS + DBLE(PT(4,I))
      ENDDO
      YSCA   = EVIS*EVIS
C---COPY MOMENTA INTO PL-ARRAY
 2    DO I=1,NT
         YCUR(I) = 0D0
         IP(I)=.TRUE.
         IJ(I)=I
         DO II=1,4
            PL(II,I)= DBLE(PT(II,I))
         ENDDO
         PM=PL(1,I)**2+PL(2,I)**2+PL(3,I)**2
         IF (PM.GT.0D0) THEN
            PL(5,I)=1D0/SQRT(PM)
         ELSE
            PL(5,I)=1D0
         ENDIF
      ENDDO
C---  FILL V-ARRAY: V(I,J) IS V(NT*(I-1)+J-I(I+1)/2)
      IAD = 0
      DO I=1,NT-1
         DO II=1,5
            PP(II)=PL(II,I)
         ENDDO
         DO J=I+1,NT
            IAD = IAD + 1
            V(IAD) = 2D0*(1D0-(PP(1)*PL(1,J) +PP(2)*PL(2,J)
     &               +PP(3)*PL(3,J))*PP(5)*PL(5,J))
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
C---  END OF CLUSTER SEARCH FOR VMINI
      IF (IMINI.NE.0) THEN
C---  NOT FINISHED YET
         ICUR       = ICUR + 1
         YCUR(ICUR) = VMINI*MIN(PL(4,IMINI),PL(4,JMINI))**2
         IF (YCUR(ICUR)+EPSILON.GE.YSCA) THEN
C---  SOFT FREEZING HERE
            IF (PL(4,IMINI).LT.PL(4,JMINI)) THEN
               IP(IMINI)=.FALSE.
            ELSE
               IP(JMINI)=.FALSE.
            ENDIF
         ELSE
C---  COMBINE PARTICLES IMINI AND JMINI
            DO II=1,4
               PL(II,IMINI)=PL(II,IMINI)+PL(II,JMINI)
            ENDDO
            PM=PL(1,IMINI)**2+PL(2,IMINI)**2+PL(3,IMINI)**2
            IF (PM.GT.0D0) THEN
               PL(5,IMINI)=1D0/SQRT(PM)
            ELSE
               PL(5,IMINI)=1D0
            ENDIF
C---  FLAG PARTICLE JMINI AS COMBINED
            IP(JMINI)=.FALSE.
            IJ(JMINI)=IMINI
            NJ=NJ-1
C---  CALCULATE RELEVANT NEW V VALUES
            DO I=1,NT
               IF (I.NE.IMINI) THEN
                  IF (IJ(I).EQ.JMINI) IJ(I)=IMINI
                  IF (IP(I)) THEN
                     K = MIN(I,IMINI)
                     L = MAX(I,IMINI)
                     IAD = NT*(K-1) + L - (K*(K+1))/2
                     V(IAD) = 2D0*(1D0-(PL(1,K)*PL(1,L) +PL(2,K)*PL(2,L)
     &                    +PL(3,K)*PL(3,L))*PL(5,K)*PL(5,L))
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
C---  BACK TO START OF LOOP
         GO TO 1
      ELSE
C
C..   ITERATION IN JET-FINDING FINISHED
C
         NJITER = NJITER + 1
         IF(NJITER.LE.NMXY) THEN
C
C..   STORE THE VALUE OF YFLIP AND THE CORRESPONDING NUMBER OF JETS
C
            YTRANS(NJITER) = YSCA
            NTRANS(NJITER) = NJ

            IF(NJITER.LE.NMXP) THEN
C
C..   STORE THE JET-DIRECTIONS AND PARTICLE ASSOCIATION
C
               J=0
               DO I=1,NT
                  IF (IJ(I).EQ.I) THEN
                     J=J+1
                     JJ(I)=J
                     DO II=1,4
                        PCMJ(NJITER,II,J)= PL(II,I)
                     ENDDO
                  ENDIF
               ENDDO
               DO I=1,NT
                  ICMJ(NJITER,I)=JJ(IJ(I))
               ENDDO
            ENDIF
            IF(NJ.LE.NJREQ) THEN
C
C..   ANOTHER ITERATION?
C
               CALL CKSORD(NT,YCUR,KIX,'I')
               DO I=1,NT
C
C..   DETERMINE THE NEXT VALUE FOR YSCA, RESET VARIABLES
C
                  IF(YCUR(KIX(I)).LT.YSCA) THEN
                     YSCA = YCUR(KIX(I))
                     ICUR = 0
                     DO J=1,NT
                        YCUR(J) = 0D0
                     ENDDO
C
C..   YUP, GO FOR THE NEXT ROUND
C
                     IF(YSCA.NE.0) GOTO 2
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDIF
C
C..   DEVIDE BY ECM**2
C
      NJMAX = -1
      DO I=1,NJITER
         YTRANS(I) = YTRANS(I)/EVIS**2
         IF(NTRANS(I).GT.NJMAX) NJMAX = NTRANS(I)
      ENDDO

      END
