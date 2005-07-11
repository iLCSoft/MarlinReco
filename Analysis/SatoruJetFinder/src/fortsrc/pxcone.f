CDECK  ID>, PXCONE.
      SUBROUTINE PXCONE (NTRAK,ITKDM,PTRAK,CONER,EPSLON,MXJET,
     +                   NJET,PJET,IPASS,IJMUL,IERR)
*.*********************************************************
*. ------
*. PXCONE
*. ------
*.
*.********** Pre Release Version 26.2.93
*.
*. Driver for the Cone  Jet finding algorithm of L.A. del Pozo.
*. Based on algorithm from D.E. Soper.
*. Finds jets inside cone of half angle CONER with energy > EPSLON.
*. Output jets are ordered in energy.
*. Usage     :
*.
*.      INTEGER  ITKDM,MXTRK
*.      PARAMETER  (ITKDM=4.or.more,MXTRK=1.or.more)
*.      INTEGER  MXJET, MXTRAK, MXPROT
*.      PARAMETER  (MXJET=10,MXTRAK=200,MXPROT=100)
*.      INTEGER  IPASS (MXTRAK),IJMUL (MXJET)
*.      INTEGER  NTRAK,NJET,IERR
*.      REAL  PTRAK (ITKDM,MXTRK),PJET (5,MXJET)
*.      REAL  CONER, EPSLON
*.      NTRAK = 1.to.MXTRAK
*.      CONER   = 0.7   (suggested value)
*.      EPSLON  = 7.0   (suggested value)
*.      CALL PXCONE (NTRAK,ITKDM,PTRAK,CONER,EPSLON,MXJET,
*.     +             NJET,PJET,IPASS,IJMUL,IERR)
*.
*. INPUT     :  NTRAK     Number of particles
*. INPUT     :  ITKDM     First dimension of PTRAK array
*. INPUT     :  PTRAK     Array of particle 4-momenta (Px,Py,Pz,E)
*. INPUT     :  CONER     Cone size (half angle) in radians
*. INPUT     :  EPSLON    Minimum Jet energy (GeV)
*. INPUT     :  MXJET     Maximum possible number of jets
*. OUTPUT    :  NJET      Number of jets found
*. OUTPUT    :  PJET      5-vectors of jets
*. OUTPUT    :  IPASS(k)  Particle k belongs to jet number IPASS(k)
*.                        IPASS = -1 if not assosciated to a jet
*. OUTPUT    :  IJMUL(i)  Jet i contains IJMUL(i) particles
*. OUTPUT    :  IERR      = 0 if all is OK ;   = -1 otherwise
*.
*. CALLS     : PXSEAR, PXSAME, PXNEW, PXTRY, PXORD, PXUVEC, PXOLAP
*. CALLED    : User
*.
*. AUTHOR    :  L.A. del Pozo
*. CREATED   :  26-Feb-93
*. LAST MOD  :  08-Oct-97
*.
*. Modification Log.
*. 08-Oct-97: D. Chrisman - Call PXADDV with the correct number
*.                          of arguments.
*. 26-Jun-96: D. Chrisman - Save ROLD and EPSOLD
*. 3-Mar-93: L A del Pozo - Check Cone size is sensible.
*. 2-Mar-93: L A del Pozo - Fix Bugs in PXOLAP
*. 1-Mar-93: L A del Pozo - Remove Cern library routine calls
*. 1-Mar-93: L A del Pozo - Add Print out of welcome and R and Epsilon
*.
*.*********************************************************
      IMPLICIT NONE
*** External Arrays
      INTEGER  ITKDM,MXJET,NTRAK,NJET,IERR
      INTEGER  IPASS (*),IJMUL (MXJET)
      REAL  PTRAK (ITKDM,*),PJET (5,*), CONER, EPSLON
*** Internal Arrays
      INTEGER MXPROT, MXTRAK
      PARAMETER (MXPROT=100, MXTRAK=200)
      REAL PP(4,MXTRAK), PU(3,MXTRAK), PJ(4,MXPROT)
      LOGICAL JETLIS(MXPROT,MXTRAK)
*** Used in the routine.
      REAL COSR,COS2R, VSEED(3), VEC1(3), VEC2(3)
      LOGICAL UNSTBL
      INTEGER I,J,N,MU,N1,N2, ITERR
      REAL ROLD, EPSOLD
      SAVE ROLD, EPSOLD
      INTEGER NCALL, NPRINT
      SAVE NCALL,NPRINT
      DATA NCALL,NPRINT /0,0/
      IERR=0
*
*** INITIALIZE
      IF(NCALL.LE.0)  THEN
         ROLD = 0.
         EPSOLD = 0.
      ENDIF
      NCALL = NCALL + 1
*
*** Print welcome and Jetfinder parameters
      IF((CONER.NE.ROLD .OR. EPSLON.NE.EPSOLD) .AND. NPRINT.LE.10) THEN
         WRITE (6,*)
         WRITE (6,*) ' *********** PXCONE: Cone Jet-finder ***********'
         WRITE(6,1000)'   Cone Size R = ',CONER,' Radians'
         WRITE(6,1001)'   Min Jet energy Epsilon = ',EPSLON,' GeV'
         WRITE (6,*) ' ***********************************************'
         WRITE (6,*)
1000     FORMAT(A18,F5.2,A10)
1001     FORMAT(A29,F5.2,A5)
         NPRINT = NPRINT + 1
         ROLD=CONER
         EPSOLD=EPSLON
      ENDIF
*** Check input Value of R is sensible.
      IF (CONER .GT. 1.5708) THEN
         WRITE (6,*) ' PXCONE: CONER > 1.57 rad (90 degrees)'
         IERR=-1
         RETURN
      ENDIF
*
*** Copy calling array PTRAK  to internal array PP(4,NTRAK)
*
      IF (NTRAK .GT. MXTRAK) THEN
         WRITE (6,*) ' PXCONE: Ntrak too large'
         IERR=-1
         RETURN
      ENDIF
      DO  100 I=1, NTRAK
         DO  101 J=1,4
            PP(J,I)=PTRAK(J,I)
101      CONTINUE
100   CONTINUE
*
*** Zero output variables
*
      NJET=0
      DO 102 I = 1, NTRAK
         DO 103 J = 1, MXPROT
           JETLIS(J,I) = .FALSE.
103      CONTINUE
102   CONTINUE
      CALL PXZERV(4*MXPROT,PJ)
      CALL PXZERV(MXJET,IJMUL)
*
      COSR = COS(CONER)
      COS2R = COS(2*CONER)
      UNSTBL = .FALSE.
      CALL PXUVEC(NTRAK,PP,PU,IERR)
      IF (IERR .NE. 0) RETURN
*** Look for jets using particle diretions as seed axes
*
      DO 110 N = 1,NTRAK
        DO 120 MU = 1,3
          VSEED(MU) = PU(MU,N)
120     CONTINUE
        CALL PXSEAR(COSR,NTRAK,PU,PP,VSEED,
     +                   NJET,JETLIS,PJ,UNSTBL,IERR)
         IF (IERR .NE. 0) RETURN
110   CONTINUE
*** Now look between all pairs of jets as seed axes.
      DO 140 N1 = 1,NJET-1
         VEC1(1)=PJ(1,N1)
         VEC1(2)=PJ(2,N1)
         VEC1(3)=PJ(3,N1)
         CALL PXNORV(3,VEC1,VEC1,ITERR)
         DO 150 N2 = N1+1,NJET
            VEC2(1)=PJ(1,N2)
            VEC2(2)=PJ(2,N2)
            VEC2(3)=PJ(3,N2)
            CALL PXNORV(3,VEC2,VEC2,ITERR)
            CALL PXADDV(3,VEC1,VEC2,VSEED)
            CALL PXNORV(3,VSEED,VSEED,ITERR)
            CALL PXSEAR(COSR,NTRAK,PU,PP,VSEED,NJET,
     +      JETLIS,PJ,UNSTBL,IERR)
            IF (IERR .NE. 0) RETURN
150      CONTINUE
140   CONTINUE
      IF (UNSTBL) THEN
        IERR=-1
        WRITE (6,*) ' PXCONE: Too many iterations to find a proto-jet'
        RETURN
      ENDIF
*** Now put the jet list into order by jet energy, eliminating jets
*** with energy less than EPSLON.
       CALL PXORD(EPSLON,NJET,NTRAK,JETLIS,PJ)
*
*** Take care of jet overlaps
       CALL PXOLAP(NJET,NTRAK,JETLIS,PJ,PP)
*
*** Order jets again as some have been eliminated, or lost energy.
       CALL PXORD(EPSLON,NJET,NTRAK,JETLIS,PJ)
*
*** All done!, Copy output into output arrays
      IF (NJET .GT. MXJET) THEN
         WRITE (6,*) ' PXCONE:  Found more than MXJET jets'
         IERR=-1
         GOTO 99
      ENDIF
      DO 300 I=1, NJET
         DO 310 J=1,4
            PJET(J,I)=PJ(J,I)
310      CONTINUE
300   CONTINUE
      DO 320 I=1, NTRAK
         IPASS(I)=-1
         DO 330 J=1, NJET
            IF (JETLIS(J,I)) THEN
               IJMUL(J)=IJMUL(J)+1
               IPASS(I)=J
            ENDIF
330      CONTINUE
320   CONTINUE
99    RETURN
      END
