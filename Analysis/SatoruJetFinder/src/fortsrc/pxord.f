CDECK  ID>, PXORD.
******............................................................******
       SUBROUTINE PXORD(EPSLON,NJET,NTRAK,JETLIS,PJ)
******............................................................******
*
*** Routine to put jets into order and eliminate tose less than EPSLON
*.
*. Modification Log.
*. 22-Apr-97: D. Chrisman - Dimension ELIST and INDEX with MXPROT
*.                          instead of 100.
*.
      IMPLICIT NONE
      INTEGER MXTRAK,MXPROT
      PARAMETER (MXTRAK=200,MXPROT=100)
      INTEGER I, J, INDEX(MXPROT)
      REAL PTEMP(4,MXPROT), ELIST(MXPROT)
      INTEGER NJET,NTRAK
      LOGICAL JETLIS(MXPROT,MXTRAK)
      LOGICAL LOGTMP(MXPROT,MXTRAK)
      REAL EPSLON,PJ(4,MXPROT)
*** Puts jets in order of energy: 1 = highest energy etc.
*** Then Eliminate jets with energy below EPSLON
*
*** Copy input arrays.
      DO 100 I=1,NJET
         DO 110 J=1,4
            PTEMP(J,I)=PJ(J,I)
110      CONTINUE
         DO 120 J=1,NTRAK
            LOGTMP(I,J)=JETLIS(I,J)
120      CONTINUE
100   CONTINUE
      DO 150 I=1,NJET
         ELIST(I)=PJ(4,I)
150   CONTINUE
*** Sort the energies...
      CALL PXSORV(NJET,ELIST,INDEX,'I')
*** Fill PJ and JETLIS according to sort ( sort is in ascending order!!)
      DO 200 I=1, NJET
         DO 210 J=1,4
            PJ(J,I)=PTEMP(J,INDEX(NJET+1-I))
210      CONTINUE
         DO 220 J=1,NTRAK
            JETLIS(I,J)=LOGTMP(INDEX(NJET+1-I),J)
220      CONTINUE
200   CONTINUE
** Jets are now in order
*** Now eliminate jets with less than Epsilon energy
      DO 300, I=1, NJET
         IF (PJ(4,I) .LT. EPSLON) THEN
            NJET=NJET-1
            PJ(4,I)=0.
         ENDIF
300   CONTINUE
      RETURN
      END
