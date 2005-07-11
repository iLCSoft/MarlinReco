CDECK  ID>, PXLCL5.
      SUBROUTINE PXLCL5 (NTRAK,ITKDM,PTRAK,IMODE,XMIN,MNJET,
     +                   MXJET,NJET,PJET,IPASS,IJMUL,IERR)
*.*********************************************************
*. ------
*. PXLCL5
*. ------
*. Jet-finding routine using the Jetset algorithm.
*. The implementation here is without a common block, however.
*. Thus this routine may be used regardless of whether the
*. Jetset6.3 or Jetset7.1 library might be linked.  It is
*. not necessary to link to Jetset, however.
*. Usage     :
*.
*.      INTEGER  ITKDM,MXTRK
*.      PARAMETER  (ITKDM=5.or.more,MXTRK=1.or.more)
*.      INTEGER  MXJET
*.      PARAMETER  (MXJET=10)
*.      INTEGER  IPASS (MXTRK),IJMUL (MXJET)
*.      INTEGER  NTRAK,NJET,IERR,IMODE,MNJET
*.      REAL  PTRAK (ITKDM,MXTRK),PJET (5,MXJET)
*.      REAL  XMIN
*.
*.      NTRAK = 1.to.MXTRAK
*.      IMODE = 1.to.4 (Jetset7.1 default = 1)
*.      XMIN  = 2.5    (Jetset7.1 default = 2.5 for IMODE = 1,2,3;
*.                                       = 0.05 for IMODE = 4)
*.      MNJET = 1.to.MXJET  (= 1 for most purposes)
*.      CALL PXLCL5 (NTRAK,ITKDM,PTRAK,IMODE,XMIN,MNJET,
*.     +             MXJET,NJET,PJET,IPASS,IJMUL,IERR)
*.
*. INPUT     : NTRAK     Total number of particles
*. INPUT     : ITKDM     First dimension of PTRAK array
*. INPUT     : PTRAK     Particle momentum array: Px,Py,Pz,E,M
*. INPUT     : IMODE     Jetfinder mode (= MSTU(46) in Jetset7.1)
*. INPUT     : XMIN      Jet resolution parameter
*.                           equivalent to PARU(44) for IMODE = 1,2,3
*.                           equivalent to PARU(45) for IMODE = 4
*. INPUT     : MNJET     The minimum number of jets to reconstruct
*. INPUT     : MXJET     The maximum number of jets permitted
*. OUTPUT    : NJET      Number of jets found
*. OUTPUT    : PJET      5-momenta of reconstructed jets
*. OUTPUT    : IPASS(k)  Particle k belongs to jet number IPASS(k)
*. OUTPUT    : IJMUL(i)  Jet i contains IJMUL(i) particles
*. OUTPUT    : IERR      = 0 if all is OK ;   = -1 otherwise
*.
*. CALLS     : PXLUC5,PXMAS4,PXPRNT,PXPRIV
*. CALLED    : By User
*.
*. AUTHOR    :  J.W.Gary
*. CREATED   :  19-Jun-88
*. LAST MOD  :  15-Feb-89
*.
*. Modification Log.
*. 15-Feb-89  Integrate with PXLUC5  J.W.Gary
*. 12-May-97 D. Chrisman, remove declaration of unused variables
*.             DMIN, MINTR and TGEN.
*.
*.*********************************************************
      IMPLICIT NONE
      INTEGER  IOLUN,NRLUDM
      PARAMETER  (IOLUN=6,NRLUDM=2000)
      INTEGER  NTRAK,MXJET,IX1,IX2,ITKDM,NJET,IERR,NLUND,
     +         MSTU46,MSTU47,IMODE,MNJET
      INTEGER  IPASS (*),IJMUL (*)
      REAL  PTRAK (ITKDM,*),PLUND (NRLUDM,5),PJET (5,*)
      REAL  XMIN,PARU44,PARU45
      LOGICAL  LPRT
      DATA  LPRT / .FALSE. /

      IERR = 0
*  select jetfinder mode
*  ------ --------- ----
      MSTU46 = IMODE
      MSTU47 = MNJET
      IF (MSTU46.LT.1.OR.MSTU46.GT.4) THEN
          WRITE (6,FMT='(''PXLCL5: Error, MSTU46 ='',I12)') MSTU46
          GO TO 995
      END IF
      PARU44 = XMIN
      PARU45 = XMIN
      IF (NTRAK.LE.1) THEN
          WRITE (IOLUN,FMT='('' PXLCL5: Error, NTRAK ='',I4)')
     +           NTRAK
          GO TO 995
      END IF
*  Pack Jetset arrays
*  ---- ------ ------
      NLUND = NTRAK
      DO 110  IX1 = 1,NTRAK
          DO 100  IX2 = 1,5
              PLUND (IX1,IX2) = PTRAK (IX2,IX1)
 100      CONTINUE
 110  CONTINUE
*  Call Jetset routine for cluster-finding
*  ---- ------ ------- --- ------- -------
      CALL PXLUC5 (NLUND,NRLUDM,PLUND,MSTU46,MSTU47,
     +             PARU44,PARU45,NJET,IJMUL,IPASS)
      IF (NJET.LT.1) THEN
          WRITE (IOLUN,FMT='('' PXLCL5: ERROR, NJET='',
     +           I6)') NJET
          IERR = -1
          RETURN
      ELSE IF (NJET.GT.MXJET) THEN
          WRITE (IOLUN,FMT='
     +       ('' PXLCL5: ERROR, NJET='',I6,
     +        '' exceeds MXJET'',I6)') NJET,MXJET
          GO TO 995
      END IF
*  Copy jet 4-momenta to output buffer
*  ---- --- - ------- -- ------ ------
      DO 210  IX2 = 1,NJET
          DO 200  IX1 = 1,4
              PJET (IX1,IX2) = PLUND (NLUND+IX2,IX1)
 200      CONTINUE
          CALL PXMAS4 (PJET (1,IX2),PJET (5,IX2))
 210  CONTINUE
      IF (LPRT) THEN
          WRITE (IOLUN,FMT='('' PXLCL5:'')')
          CALL PXPRNT (1,NJET,5,5,PJET,'E')
          CALL PXPRIV ('IJMUL',NJET,IJMUL)
          CALL PXPRIV ('IPASS',NTRAK,IPASS)
      END IF

      RETURN
 995  IERR = -1
      RETURN
      END
