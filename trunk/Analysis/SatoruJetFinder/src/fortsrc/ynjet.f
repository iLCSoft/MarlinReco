CDECK  ID>, YNJET.
      SUBROUTINE YNJET(YCUT,NJET,IERR)
C
C  ROUTINE TO RETURN THE NUMBER OF JETS FOR A GIVEN YCUT
C
C  LAST MOD  :  21-Jun-98
C
C  Modification Log.
C  08-Oct-97 D. Chrisman, Increase NYCLMX from 250 to 500.
C  21-Jun-98 D. Chrisman, Want to be able to handle more than 10 jets.
C                         Introduce parameter NJETMX. Use this
C                         in declaration of PJET and YREC in order to remove
C                         hard coded numbers. Also replace the hard
C                         coded number in the "D0 5001" loop with the
C                         parameter NJETMX.
C
      IMPLICIT NONE
C     IMPLICIT NONE
      INTEGER NJET,IERR,NYCLMX,NJETMX,IMODEO,NJETO,NTO,I
      PARAMETER (NYCLMX = 500, NJETMX = 20)
      REAL YCUT,YREC(NJETMX),PINT(10,NYCLMX),PJET(10,NJETMX,NJETMX)
      INTEGER HISTOR(2,NYCLMX)
      COMMON /YCL/YREC,PJET,HISTOR
      COMMON /YINT/ IMODEO,NTO,NJETO,PINT
      IERR = -1
C
CHECK IF CALL WAS MADE TO YKERN
C
      IF(IMODEO.LE.0 .OR. IMODEO.GT.7) THEN
        WRITE(6,111)
 111    FORMAT(' #### YNJET: YKERN MUST BE CALLED FIRST ! ####')
        NJET = -1
        RETURN
      ENDIF
C
CHECK IF INPUT MAKES SENSE
C
      IF(YCUT.LE.0. .OR. YCUT.GT.1.) THEN
        WRITE(6,1) YCUT
 1      FORMAT(' #### YNJET: INPUT YCUT=',E12.4,' INVALID ####')
        NJET = -1
        RETURN
      ENDIF
C
      NJET = 1
      DO 5001 I=1,NJETMX
        IF(YCUT.LT.YREC(I)) NJET = I+1
 5001 CONTINUE
      IERR = 0
      RETURN
      END
