CDECK  ID>, YYJET.
      SUBROUTINE YYJET(NJET,YL,YH,IERR)
C
C  ROUTINE TO RETURN THE VALUES OF YCUT BETWEEN WHICH EVENT IS
C  CLASSIFIED AS N-JET (YL < YH)
C
C  LAST MOD  :  21-Jun-98
C
C  Modification Log.
C  08-Oct-97 D. Chrisman, Increase NYCLMX from 250 to 500.
C  21-Jun-98 D. Chrisman, Want to be able to handle more than 10 jets.
C                         Introduce parameter NJETMX. Use this
C                         in declaration of PJET and YREC in order to remove
C                         hard coded numbers.
C
      IMPLICIT NONE
C     IMPLICIT NONE
      INTEGER NYCLMX,NJETMX,NJET,IERR,IMODEO,NJETO,NTO
      PARAMETER (NYCLMX = 500, NJETMX = 20)
      INTEGER HISTOR(2,NYCLMX)
      REAL YREC(NJETMX),YL,YH,PINT(10,NYCLMX),PJET(10,NJETMX,NJETMX)
      COMMON /YCL/YREC,PJET,HISTOR
      COMMON /YINT/ IMODEO,NTO,NJETO,PINT
      IERR = -1
C
CHECK IF CALL WAS MADE TO YKERN
C
      IF(IMODEO.LE.0 .OR. IMODEO.GT.7) THEN
        WRITE(6,111)
 111    FORMAT(' #### YYJET: YKERN MUST BE CALLED FIRST ! ####')
        YL = -1.
        YH = -1.
        RETURN
      ENDIF
C
CHECK IF INPUT MAKES SENSE
C
      IF(NJET.LE.0 .OR. NJET.GT.NJETMX) THEN
        WRITE(6,1) NJET
 1      FORMAT(' #### YYJET: REQUEST FOR NJET=',I12,
     +  ' NOT SUPPORTED ####')
        YL = -1.
        YH = -1.
        RETURN
      ENDIF
      IF(NJET.GT.NTO) THEN
        WRITE(6,2) NJET,NTO
 2      FORMAT(' #### YYJET:',I3,' JETS OUT OF',I3,' PARTICLES NOT',
     +  ' POSSIBLE. ####')
        YL = -1.
        YH = -1.
        RETURN
      ENDIF
C
      IF(NJET.EQ.1) THEN
        YH = 1.
      ELSE
        YH = YREC(NJET-1)
      ENDIF
      YL = YREC(NJET)
C
      IERR = 0
      RETURN
      END
