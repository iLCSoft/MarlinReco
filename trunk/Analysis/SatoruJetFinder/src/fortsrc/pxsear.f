CDECK  ID>, PXSEAR.
******............................................................******
      SUBROUTINE PXSEAR(COSR,NTRAK,PU,PP,VSEED,NJET,
     +                JETLIS,PJ,UNSTBL,IERR)
******............................................................******
*
      IMPLICIT NONE
      INTEGER MXTRAK, MXPROT
      PARAMETER (MXTRAK=200,MXPROT=100)
      INTEGER NTRAK, IERR
      REAL COSR,PU(3,MXTRAK),PP(4,MXTRAK),VSEED(3)
      LOGICAL UNSTBL
      LOGICAL JETLIS(MXPROT,MXTRAK)
      INTEGER NJET
      REAL  PJ(4,MXPROT)
*** Using VSEED as a trial axis , look for a stable jet.
*** Check stable jets against those already found and add to PJ.
*** Will try up to MXITER iterations to get a stable set of particles
*** in the cone.
      INTEGER MU,N,ITER
      LOGICAL PXSAME,PXNEW,OK
      LOGICAL NEWLIS(MXTRAK),OLDLIS(MXTRAK)
      REAL OAXIS(3),NAXIS(3),PNEW(4)
      INTEGER MXITER
      PARAMETER(MXITER = 30)
*
      DO 100 MU=1,3
        OAXIS(MU) = VSEED(MU)
100   CONTINUE
      DO 110 N = 1,NTRAK
        OLDLIS(N) = .FALSE.
110   CONTINUE
      DO 120 ITER = 1,MXITER
        CALL PXTRY(COSR,NTRAK,PU,PP,OAXIS,NAXIS,PNEW,NEWLIS,OK)
*** Return immediately if there were no particles in the cone.
       IF (.NOT.OK) THEN
         RETURN
       ENDIF
       IF(PXSAME(NEWLIS,OLDLIS,NTRAK)) THEN
*** We have a stable jet.
             IF (PXNEW(NEWLIS,JETLIS,NTRAK,NJET)) THEN
*** And the jet is a new one. So add it to our arrays.
*** Check arrays are big anough...
             IF (NJET .EQ. MXPROT) THEN
             WRITE (6,*) ' PXCONE:  Found more than MXPROT proto-jets'
                IERR = -1
                RETURN
             ENDIF
               NJET = NJET + 1
               DO 130 N = 1,NTRAK
                 JETLIS(NJET,N) = NEWLIS(N)
130            CONTINUE
               DO 140 MU=1,4
                 PJ(MU,NJET)=PNEW(MU)
140          CONTINUE
             ENDIF
             RETURN
       ENDIF
*** The jet was not stable, so we iterate again
       DO 150 N=1,NTRAK
         OLDLIS(N)=NEWLIS(N)
150    CONTINUE
       DO 160 MU=1,3
         OAXIS(MU)=NAXIS(MU)
160    CONTINUE
120   CONTINUE
      UNSTBL = .TRUE.
      RETURN
      END
