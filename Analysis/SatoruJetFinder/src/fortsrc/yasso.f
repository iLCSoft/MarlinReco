CDECK  ID>, YASSO.
      SUBROUTINE YASSO(NJET,PNJ,BL,IERR)
C
C  ROUTINE TO RETURN THE ASSIGNMENT OF PARTICLES TO JET AXES, FOR
C  CLASSIFICATION AS N-JET (FOR PARTICLE K IN MOMENTUM ARRAY,
C  J=BL(K) POINTS TO JET NUMBER J IN ARRAY PNJ).
C  JETS ARE ORDERED ACCORDING TO THEIR ENERGIES: E1 >= E2 >= ...
C
C  LAST MOD  : 28-Jul-99
C
C  Modification Log.
C  28-Jul-99 D. Chrisman, ITAG and IREORD are now dimensioned with NJETMX
C                         elements.
C  08-Oct-97 D. Chrisman, Increase NYCLMX from 250 to 500.
C  21-Jun-98 D. Chrisman, Want to be able to handle more than 10 jets.
C                         Introduce parameter NJETMX. Use this
C                         in declaration of PJET, YREC and PNJ in order to
C                         remove hard coded numbers.
C
      IMPLICIT NONE
C     IMPLICIT NONE
      INTEGER NYCLMX,NJETMX,NJET,IERR,IMODEO,NJETO,NTO,ICHECK,I,J,K,N
      INTEGER I1,I2
      PARAMETER (NYCLMX = 500, NJETMX = 20)
      INTEGER HISTOR(2,NYCLMX),IMAX,BL(NYCLMX),ITAG(NJETMX),
     +         IREORD(NJETMX)
      REAL YREC(NJETMX),PINT(10,NYCLMX),PJET(10,NJETMX,NJETMX)
      REAL PNJ(10,NJETMX),EMAX
      COMMON /YCL/YREC,PJET,HISTOR
      COMMON /YINT/ IMODEO,NTO,NJETO,PINT
      IERR = -1
C
CHECK IF CALL WAS MADE TO YKERN
C
      ICHECK = 0
      IF(IMODEO.LE.0 .OR. IMODEO.GT.7) THEN
        WRITE(6,111)
 111    FORMAT(' #### YASSO: YKERN MUST BE CALLED FIRST ! ####')
        ICHECK = -1
      ENDIF
C
CHECK IF INPUT MAKES SENSE
C
      IF(NJET.LE.0 .OR. NJET.GT.NJETMX) THEN
        WRITE(6,1) NJET
 1      FORMAT(' #### YASSO: REQUEST FOR NJET=',I4,
     +  ' NOT SUPPORTED ####')
        ICHECK = -1
      ENDIF
      IF(NJET.GT.NTO) THEN
        WRITE(6,2) NJET,NTO
 2      FORMAT(' #### YASSO:',I4,' JETS OUT OF',I4,' PARTICLES NOT',
     +  ' POSSIBLE. ####')
        ICHECK = -1
      ENDIF
      IF(ICHECK.NE.0) THEN
        DO 5001 I=1,10
          DO 5002 J=1,NJETMX
            PNJ(I,J) = -1.
 5002     CONTINUE
 5001   CONTINUE
        DO 5003 I=1,NTO
          BL(I) = -1
 5003   CONTINUE
        RETURN
      ENDIF
C
      DO 5004 I=1,NTO
        BL(I) = I
 5004 CONTINUE
C
      IF(NJET.NE.NTO) THEN
      DO 5005 I=NTO,NJET+1,-1
        I1 = HISTOR(1,I)
        I2 = HISTOR(2,I)
        DO 5006 N=1,NTO
          IF(BL(N).EQ.I1) BL(N) = I2
          IF(I1.NE.I) THEN
            IF(BL(N).EQ.I) BL(N) = I1
          ENDIF
 5006   CONTINUE
 5005 CONTINUE
      ENDIF
C
C ORDER JETS ACCORDING TO THEIR ENERGY (FIRST IS LARGEST)
C AND CHANGE POINTERS ACCORDINGLY
C
      DO 5007 I=1,NJET
        ITAG(I) = 1
 5007 CONTINUE
      DO 5008 I=1,NJET
C       IF(ITAG(I).NE.0) THEN
          EMAX = 0.
          IMAX = 0
          DO 5009 J=1,NJET
            IF(ITAG(J).NE.0 .AND. EMAX.LT.PJET(4,J,NJET)) THEN
               EMAX = PJET(4,J,NJET)
               IMAX = J
            ENDIF
 5009     CONTINUE
          IF(IMAX.LE.0) THEN
            WRITE(6,9)
 9          FORMAT(' #### YASSO: JET AXIS WITH ZERO OR NEGATIVE ',
     +      'ENERGY COMPONENT DETECTED; NO ORDERING DONE. ####')
            DO 5010 N=1,NJET
              DO 5011 K=1,7
                PNJ(K,N) = PJET(K,N,NJET)
 5011         CONTINUE
 5010       CONTINUE
            RETURN
          ENDIF
          ITAG(IMAX) = 0
          IREORD(IMAX) = I
C       ENDIF
 5008 CONTINUE
C
      DO 5012 N=1,NJET
        DO 5013 K=1,7
          PNJ(K,IREORD(N)) = PJET(K,N,NJET)
 5013   CONTINUE
 5012 CONTINUE
C
      DO 5014 I=1,NTO
        BL(I) = IREORD(BL(I))
 5014 CONTINUE
C
      IERR = 0
      RETURN
      END
