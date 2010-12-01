CDECK  ID>, SYJKRN.
**************************************************************************
      SUBROUTINE SYJKRN(MD,NJETRQ,THRESH,IMODE,
     +  YCUTS,IEXAM,FRAC,IEXAM4,NPAR,IDIM,PPAR,MXJET,
     +  NJET,IASS,JDIM,PJETS,Y34,Y45,IERR)
* Kernel Routine for SY-Jetfindings Scheme or Traditional Jetfindings
* Authour   : Satoru Yamashita
* Creat     : 5-SEP-98
*           :
* Inputs:
*   MD      ; mode (see below)
*   NJETRQ  : fixed number of jets ; =0 variable -- use YCUT etc...
*   THRESH  : threshold (GeV) for the primary jet finding
*   IMODE   : Kernel Jet finding mode
*      1- 6 : YKERN
*       IF(IMODE.EQ.1) THEN
*          CM = 'JADE E0'
*        ELSEIF(IMODE.EQ.2) THEN
*          CM = 'JADE P '
*        ELSEIF(IMODE.EQ.3) THEN
*          CM = 'JADE P0'
*        ELSEIF(IMODE.EQ.4) THEN
*          CM = 'JADE E '
*        ELSEIF(IMODE.EQ.5) THEN
*          CM = 'DURHAM '
*        ELSEIF(IMODE.EQ.6) THEN
*          CM = 'GENEVA '
*        ELSE
*          WRITE(6,281) IMODE
* 281      FORMAT(/,' ### YKERN: IMODE =',I3,' INVALID; SET TO 1 ###')
*          IMODE = 1
*          CM = 'JADE E0'
*        ENDIF
*
*
*
*         7 : cambridge
*      9-11 : LUCLUS
*        12 : cone
*   YCUTS   : (1) ycut or xcut or Cone Radius (2) EPS
*   IEXAM   : Particle re-association method 1-6
*            EXAM(1) = THETA
*            EXAM(2) = PART(4)*PJSTMP(4,J)-PART(1)*PJSTMP(1,J)
*     +                -PART(2)*PJSTMP(2,J)-PART(3)*PJSTMP(3,J)
*            EXAM(3) = DURHAM
*            EXAM(4) = E0JADE
*            EXAM(5) = Geneva Scheme
*            EXAM(6) = EJADE   ! it makes too bad cross talk
*
*   FRAC    : core rate
*   IEXAM4  : Jet merge method 1-8
*            EXAM(1,N) = PL(4,I)*PL(4,J)
*            EXAM(2,N) = PL(4,I)*PL(4,J)-PL(1,I)*PL(1,J)
*     +           -PL(2,I)*PL(2,J)-PL(3,I)*PL(3,J)
*            EXAM(3,N) = DURHAM(I,J)
*            EXAM(4,N) = E0JADE(I,J)
*            EXAM(5,N) = EJADE(I,J)
*            EXAM(6,N) = ANGLE(I,J)
*            EXAM(7,N) = ECOS(I,J)
*            EXAM(8,N) = PL(6,I)*PL(6,J)
*   NPAR    : Number of "particles"
*   IDIM    : first dimension of PPAR array
*   PPAR    : Particle Momenta
*   MXJET   : Maximum Number of Jets (buffer size for PJETS = MXJETS*JDIM)
*   JDIM    : first dimension of PJETS array
*
* Outputs:
*   NJET    : Number of Jets
*   IASS    : Jet assosiation for "particle"s
*   PJETS   : Jet momenta
*   Y34     : YHI (sometimes not filled)
*   Y45     : YLO (sometimes not filled)
*   IERR    : Error flag  0==All O.K.
*     -9999: abnormal error
*      -999: boundary error
*       -99: input mismatch
*      Others : error from various internal calls
* ==========================================================================
*   MD
* There are 11 variations
* MD   Name      : Final NJ    : 1st process  : 2nd process : 3rd process
* 0a. Traditional: variable    :  fixed Ycuts :     -       :   -
* 0b. Traditional: fixed Njets :  fixed Njets :     -       :   -
* 0c. Just Merge : fixed Njets :  fixed Ycuts : do merge    :   -
* 1a. LMODE1     : variable    :  fixed Ycuts : do SYJETF1  :   -
* 1b. LMODE1     : fixed final :  fixed Ycuts : do SYJETF1  : do merge
* 1c. LMODE1     : fixed final :  fixed Njets : do SYJETF1  :   -
* 1d. LMODE1     : fixed final :  fixed Ycuts : do merge    : do SYJETF1
* 2a. LMODE2     : variable    :  fixed Ycuts : do SYJETF2  :   -
* 2b. LMODE2     : fixed final :  fixed Ycuts : do SYJETF2  : do merge
* 2c. LMODE2     : fixed final :  fixed Njets : do SYJETF2  :   -
* 2d. LMODE2     : fixed final :  fixed Ycuts : do merge    : do SYJETF2
* ============================================================================
*  MD & NJETRQ &  IMODE  & YCUTS & IEXAM &  Frac   &  IEXAM4 &     ISYMOD
*  0a &   -    &   1-12  &  > 0  &   -   &   -     &     -   &        1
*  0b &  > 0   &   1-12  &   -   &   -   &   -     &     -   &        2
*  0c &  > 0   &   1-12  &  > 0  &   -   &   -     &    1-8  &        3
*  1a &   -    &   1-12  &  > 0  &  1-6  & 0.0-1.0 &     -   &       11
*  1b &  > 0   &   1-12  &  > 0  &  1-6  & 0.0-1.0 &    1-8  &       12
*  1c &  > 0   &   1-12  &   -   &  1-6  & 0.0-1.0 &     -   &       13
*  1d &  > 0   &   1-12  &  > 0  &  1-6  & 0.0-1.0 &    1-8  &       14
*  2a &   -    &   1-12  &  > 0  &  1-6  & 0.0-1.0 &     -   &       21
*  2b &  > 0   &   1-12  &  > 0  &  1-6  & 0.0-1.0 &    1-8  &       22
*  2c &  > 0   &   1-12  &   -   &  1-6  & 0.0-1.0 &     -   &       23
*  2d &  > 0   &   1-12  &  > 0  &  1-6  & 0.0-1.0 &    1-8  &       24
* ============================================================================
****************************************************************************
      implicit none
      CHARACTER*(*) MD
      INTEGER IDIM,JDIM,MXJET
      INTEGER NPAR
      REAL    PPAR(IDIM,*)
      REAL    PJETS(JDIM,*)
      INTEGER NJETRQ,IMODE,IEXAM,IEXAM4,NJET,IASS(*),IERR
      REAL    THRESH,YCUTS(*),FRAC,Y34,Y45

      INTEGER MXWRK
      PARAMETER (MXWRK = 300)
      INTEGER NPAR0
      REAL    PPAR0(5,MXWRK)
      INTEGER    IUSE(MXWRK)
      REAL    YCUT,R0,EPS0

      INTEGER I,J
      INTEGER IJ(MXWRK),IJN(MXWRK),IJJ(MXWRK)
      REAL    PJ6(6,50),PJN6(6,50),PTMP6(6,50)
      REAL    PORD(50),RMAX
      INTEGER IIII(50),II
      REAL    YDHI,YDLO
      INTEGER IMODE0
      INTEGER ISYMOD

*-----Initialysation for outputs
      IERR = -9999
      NJET = 0
      DO I=1,NPAR
         IASS(I) = -1
      ENDDO
      CALL VZERO(PJETS,MXJET*JDIM)
      Y34 = 0.0
      Y45 = 0.0
      ISYMOD = 0

*-----Check boundary
      IF(MD(1:2).EQ.'0a'.OR.MD(1:2).EQ.'0A') ISYMOD = 1
      IF(MD(1:2).EQ.'0b'.OR.MD(1:2).EQ.'0B') ISYMOD = 2
      IF(MD(1:2).EQ.'0c'.OR.MD(1:2).EQ.'0C') ISYMOD = 3
      IF(MD(1:2).EQ.'1a'.OR.MD(1:2).EQ.'1A') ISYMOD = 11
      IF(MD(1:2).EQ.'1b'.OR.MD(1:2).EQ.'1B') ISYMOD = 12
      IF(MD(1:2).EQ.'1c'.OR.MD(1:2).EQ.'1C') ISYMOD = 13
      IF(MD(1:2).EQ.'1d'.OR.MD(1:2).EQ.'1D') ISYMOD = 14
      IF(MD(1:2).EQ.'2a'.OR.MD(1:2).EQ.'2A') ISYMOD = 21
      IF(MD(1:2).EQ.'2b'.OR.MD(1:2).EQ.'2B') ISYMOD = 22
      IF(MD(1:2).EQ.'2c'.OR.MD(1:2).EQ.'2C') ISYMOD = 23
      IF(MD(1:2).EQ.'2d'.OR.MD(1:2).EQ.'2D') ISYMOD = 24
*      write(*,*)'symode selected: ',ISYMOD
*
      IF(ISYMOD.EQ.0) THEN
         WRITE(*,'(A,A)')'unknown MODE detected:'//MD(1:2)
*         CALL REPORT('SYJKRN',1,'E')
         IERR  = -999
         RETURN
      ENDIF
      IF(IDIM.LT.4) THEN
         WRITE(*,'(A,I5)')'too small IDIM',IDIM
*         CALL REPORT('SYJKRN',2,'E')
         IERR  = -999
         RETURN
      ENDIF
      IF(NPAR.GT.MXWRK) THEN
         WRITE(*,'(A,I5)')'too many particles',NPAR
*         CALL REPORT('SYJKRN',3,'E')
         IERR  = -999
         RETURN
      ENDIF
      IF(IMODE.LE.0.OR.IMODE.GT.12) THEN
         WRITE(*,'(A,I5)')'unknown Jet find mode',IMODE
*         CALL REPORT('SYJKRN',4,'E')
         IERR = -999
         RETURN
      ENDIF

*-----check inputs
      IERR = -99
      IF(ISYMOD.EQ.1) THEN
         IF(YCUTS(1).LE.0.0) THEN
            WRITE(*,'(A,I6,E12.3)')'ycuts small?',ISYMOD,YCUTS(1)
*            CALL REPORT('SYJKRN',5,'E')
            RETURN
         ENDIF
      ELSEIF(ISYMOD.EQ.2) THEN
         IF(NJETRQ.LE.0) THEN
            WRITE(*,'(A,I6,I7)')'No NjetRQ',ISYMOD,NJETRQ
*            CALL REPORT('SYJKRN',6,'E')
            RETURN
         ENDIF
      ELSEIF(ISYMOD.EQ.3) THEN
         IF(NJETRQ.LE.0) THEN
            WRITE(*,'(A,I6,I7)')'No NjetRQ',ISYMOD,NJETRQ
*            CALL REPORT('SYJKRN',7,'E')
            RETURN
         ENDIF
         IF(YCUTS(1).LE.0.0) THEN
            WRITE(*,'(A,I6,E12.3)')'ycuts small?',ISYMOD,YCUTS(1)
*            CALL REPORT('SYJKRN',8,'E')
            RETURN
         ENDIF
         IF(IEXAM4.LE.0.OR.IEXAM4.GT.8) THEN
            WRITE(*,'(A,I6,I7)')'No merge way',ISYMOD,IEXAM4
*            CALL REPORT('SYJKRN',9,'E')
            RETURN
         ENDIF
      ELSEIF(ISYMOD.EQ.11.OR.ISYMOD.EQ.21) THEN
         IF(YCUTS(1).LE.0.0) THEN
            WRITE(*,'(A,I6,E12.3)')'ycuts small?',ISYMOD,YCUTS(1)
*            CALL REPORT('SYJKRN',10,'E')
            RETURN
         ENDIF
         IF(IEXAM.LE.0.OR.IEXAM.GT.6) THEN
            WRITE(*,'(A,I6,I7)')'No re-assign way',ISYMOD,IEXAM
*            CALL REPORT('SYJKRN',11,'E')
            RETURN
         ENDIF
      ELSEIF(ISYMOD.EQ.12.OR.ISYMOD.EQ.22) THEN
         IF(NJETRQ.LE.0) THEN
            WRITE(*,'(A,I6,I7)')'No NjetRQ',ISYMOD,NJETRQ
*            CALL REPORT('SYJKRN',12,'E')
            RETURN
         ENDIF
         IF(YCUTS(1).LE.0.0) THEN
            WRITE(*,'(A,I6,E12.3)')'ycuts small?',ISYMOD,YCUTS(1)
*            CALL REPORT('SYJKRN',13,'E')
            RETURN
         ENDIF
         IF(IEXAM4.LE.0.OR.IEXAM4.GT.8) THEN
            WRITE(*,'(A,I6,I7)')'No merge way',ISYMOD,IEXAM4
*            CALL REPORT('SYJKRN',14,'E')
            RETURN
         ENDIF
         IF(IEXAM.LE.0.OR.IEXAM.GT.6) THEN
            WRITE(*,'(A,I6,I7)')'No re-assign way',ISYMOD,IEXAM
*            CALL REPORT('SYJKRN',15,'E')
            RETURN
         ENDIF
      ELSEIF(ISYMOD.EQ.13.OR.ISYMOD.EQ.23) THEN
         IF(NJETRQ.LE.0) THEN
            WRITE(*,'(A,I6,I7)')'No NjetRQ',ISYMOD,NJETRQ
*            CALL REPORT('SYJKRN',16,'E')
            RETURN
         ENDIF
         IF(IEXAM.LE.0.OR.IEXAM.GT.6) THEN
            WRITE(*,'(A,I6,I7)')'No re-assign way',ISYMOD,IEXAM
*            CALL REPORT('SYJKRN',17,'E')
            RETURN
         ENDIF
      ELSEIF(ISYMOD.EQ.14.OR.ISYMOD.EQ.24) THEN
         IF(NJETRQ.LE.0) THEN
            WRITE(*,'(A,I6,I7)')'No NjetRQ',ISYMOD,NJETRQ
*            CALL REPORT('SYJKRN',18,'E')
            RETURN
         ENDIF
         IF(YCUTS(1).LE.0.0) THEN
            WRITE(*,'(A,I6,E12.3)')'ycuts small?',ISYMOD,YCUTS(1)
*            CALL REPORT('SYJKRN',19,'E')
            RETURN
         ENDIF
         IF(IEXAM4.LE.0.OR.IEXAM4.GT.8) THEN
            WRITE(*,'(A,I6,I7)')'No merge way',ISYMOD,IEXAM4
*            CALL REPORT('SYJKRN',20,'E')
            RETURN
         ENDIF
         IF(IEXAM.LE.0.OR.IEXAM.GT.6) THEN
            WRITE(*,'(A,I6,I7)')'No re-assign way',ISYMOD,IEXAM
*            CALL REPORT('SYJKRN',21,'E')
            RETURN
         ENDIF
      ENDIF

*-----Preparation
      NPAR0 = NPAR
      DO I=1,NPAR0
         CALL UCOPY(PPAR(1,I),PPAR0(1,I),MIN(5,IDIM))
         IF(IDIM.LT.5) THEN
            PPAR0(5,I)=PPAR0(4,I)**2-
     +           PPAR0(1,I)**2-PPAR0(2,I)**2-PPAR0(3,I)**2
            PPAR0(5,I)=SQRT(MAX(0.0,PPAR0(5,I)))
         ENDIF
      ENDDO
*      IF(MOD(IMODE0,100).EQ.12) THEN
      IF (IMODE.EQ.12) THEN    ! bugfix JS
         YCUT = 0.0
         R0   = YCUTS(1)
         EPS0 = YCUTS(2)
      ELSE
         YCUT = YCUTS(1)
         R0   = 0.0
         EPS0 = 0.0
      ENDIF

*-----Check threshold
      DO I=1,NPAR0
         IF(PPAR0(4,I).GT.THRESH) THEN
            IUSE(I) = 1
         ELSE
            IUSE(I) = 0
         ENDIF
      ENDDO

*###### S T A R T #######
*-----Traditional Jet findings variable Njet
      IF(ISYMOD.EQ.1) THEN
         IMODE0 = IMODE
         CALL SYJETF0(IMODE0,YCUT,EPS0,R0,0,NPAR0,PPAR0,IUSE,
     +        NJET,IJ,PTMP6,YDHI,YDLO,IERR)
         IF(IERR.NE.0) GOTO 999
         CALL SYJCP6(NPAR0,PPAR0,NJET,IJ,6,PJ6)
*-----Traditional Jet findings fixed Njet
      ELSEIF(ISYMOD.EQ.2) THEN
         IMODE0 = IMODE+100
         CALL SYJETF0(IMODE0,YCUT,EPS0,R0,NJETRQ,NPAR0,PPAR0,IUSE,
     +        NJET,IJ,PTMP6,Y34,Y45,IERR)
         IF(IERR.NE.0) GOTO 999
         CALL SYJCP6(NPAR0,PPAR0,NJET,IJ,6,PJ6)
*-----Traditional Jet findings variable Njet and then merge if Njet>Njetrq
      ELSEIF(ISYMOD.EQ.3) THEN
         IMODE0 = IMODE
         CALL SYJETF0(IMODE0,YCUT,EPS0,R0,0,NPAR0,PPAR0,IUSE,
     +        NJET,IJN,PTMP6,YDHI,YDLO,IERR)
         IF(IERR.NE.0) NJET = 0
         IF(NJET.LT.NJETRQ) THEN
            IMODE0 = IMODE + 100
            CALL SYJETF0(IMODE0,YCUT,EPS0,R0,NJETRQ,NPAR0,PPAR0,IUSE,
     +           NJET,IJN,PTMP6,Y34,Y45,IERR)
            IF(IERR.NE.0) GOTO 999
         ENDIF
         CALL SYJCP6(NPAR0,PPAR0,NJET,IJN,6,PJN6)
         CALL SYMRGJ(NJETRQ,IEXAM4,NJET,NPAR0,IJN,PJN6,IJ,PJ6,IERR)  ! merge
         IF(IERR.NE.0) GOTO 999
         NJET = NJETRQ
*-----Re-assignment Jet findings variable Njet ; no merge
      ELSEIF(ISYMOD.EQ.11.OR.ISYMOD.EQ.21) THEN
         IMODE0 = IMODE
         CALL SYJETF0(IMODE0,YCUT,EPS0,R0,0,NPAR0,PPAR0,IUSE,
     +        NJET,IJN,PTMP6,YDHI,YDLO,IERR)
         IF(IERR.NE.0) GOTO 999
         IF(ISYMOD.LT.20) THEN
*          3 times iteration
            CALL SYJETF1(IEXAM,NJET,NPAR0,PPAR0,IJN,FRAC,IJ ,PJ6)
            CALL SYJETF1(IEXAM,NJET,NPAR0,PPAR0,IJ ,FRAC,IJN,PJ6)
            CALL SYJETF1(IEXAM,NJET,NPAR0,PPAR0,IJN,FRAC,IJ ,PJ6)
         ELSE  ! 21
*          3 times iteration
            CALL SYJETF2(IEXAM,NJET,NPAR0,PPAR0,IJN,FRAC,IJ ,PJ6)
            CALL SYJETF2(IEXAM,NJET,NPAR0,PPAR0,IJ ,FRAC,IJN,PJ6)
            CALL SYJETF2(IEXAM,NJET,NPAR0,PPAR0,IJN,FRAC,IJ ,PJ6)
         ENDIF
*-----Re-assignment Jet findings variable Njet ;  merge if njet>njetrq
      ELSEIF(ISYMOD.EQ.12.OR.ISYMOD.EQ.22) THEN
         IMODE0 = IMODE
         CALL SYJETF0(IMODE0,YCUT,EPS0,R0,0,NPAR0,PPAR0,IUSE,
     +        NJET,IJN,PTMP6,YDHI,YDLO,IERR)
         IF(IERR.NE.0) NJET = 0
         IF(NJET.LT.NJETRQ) THEN
            IMODE0 = IMODE + 100
            CALL SYJETF0(IMODE0,YCUT,EPS0,R0,NJETRQ,NPAR0,PPAR0,IUSE,
     +           NJET,IJN,PTMP6,Y34,Y45,IERR)
            IF(IERR.NE.0) GOTO 999
         ENDIF
         IF(ISYMOD.LT.20) THEN
*          3 times iteration
            CALL SYJETF1(IEXAM,NJET,NPAR0,PPAR0,IJN,FRAC,IJJ,PJN6)
            CALL SYJETF1(IEXAM,NJET,NPAR0,PPAR0,IJJ,FRAC,IJN,PJN6)
            CALL SYJETF1(IEXAM,NJET,NPAR0,PPAR0,IJN,FRAC,IJJ,PJN6)
         ELSE  ! 21
*          3 times iteration
            CALL SYJETF2(IEXAM,NJET,NPAR0,PPAR0,IJN,FRAC,IJJ,PJN6)
            CALL SYJETF2(IEXAM,NJET,NPAR0,PPAR0,IJJ,FRAC,IJN,PJN6)
            CALL SYJETF2(IEXAM,NJET,NPAR0,PPAR0,IJN,FRAC,IJJ,PJN6)
         ENDIF
*        Merge
         CALL SYMRGJ(NJETRQ,IEXAM4,NJET,NPAR0,IJJ,PJN6,IJ,PJ6,IERR)
         IF(IERR.NE.0) GOTO 999
         NJET = NJETRQ
*-----Re-assignment Jet findings variable Njet ; no merge
      ELSEIF(ISYMOD.EQ.13.OR.ISYMOD.EQ.23) THEN
         IMODE0 = IMODE + 100
         CALL SYJETF0(IMODE0,YCUT,EPS0,R0,NJETRQ,NPAR0,PPAR0,IUSE,
     +        NJET,IJN,PTMP6,Y34,Y45,IERR)
         IF(IERR.NE.0) GOTO 999
         IF(ISYMOD.LT.20) THEN
*          3 times iteration
            CALL SYJETF1(IEXAM,NJET,NPAR0,PPAR0,IJN,FRAC,IJ ,PJ6)
            CALL SYJETF1(IEXAM,NJET,NPAR0,PPAR0,IJ ,FRAC,IJN,PJ6)
            CALL SYJETF1(IEXAM,NJET,NPAR0,PPAR0,IJN,FRAC,IJ ,PJ6)
         ELSE  ! 21
*          3 times iteration
            CALL SYJETF2(IEXAM,NJET,NPAR0,PPAR0,IJN,FRAC,IJ ,PJ6)
            CALL SYJETF2(IEXAM,NJET,NPAR0,PPAR0,IJ ,FRAC,IJN,PJ6)
            CALL SYJETF2(IEXAM,NJET,NPAR0,PPAR0,IJN,FRAC,IJ ,PJ6)
         ENDIF
*-----Re-assignment Jet findings variable Njet ;  merge if njet>njetrq
      ELSEIF(ISYMOD.EQ.14.OR.ISYMOD.EQ.24) THEN
         IMODE0 = IMODE
         CALL SYJETF0(IMODE0,YCUT,EPS0,R0,0,NPAR0,PPAR0,IUSE,
     +        NJET,IJN,PTMP6,YDHI,YDLO,IERR)
         IF(IERR.NE.0) NJET = 0
         IF(NJET.LT.NJETRQ) THEN
            IMODE0 = IMODE + 100
            CALL SYJETF0(IMODE0,YCUT,EPS0,R0,NJETRQ,NPAR0,PPAR0,IUSE,
     +           NJET,IJN,PTMP6,Y34,Y45,IERR)
            IF(IERR.NE.0) GOTO 999
         ENDIF
         CALL SYJCP6(NPAR0,PPAR0,NJET,IJN,6,PJN6)
         CALL SYMRGJ(NJETRQ,IEXAM4,NJET,NPAR0,IJN,PJN6,IJJ,PJN6,IERR)
         IF(IERR.NE.0) GOTO 999
         NJET = NJETRQ
         IF(ISYMOD.LT.20) THEN
*          3 times iteration
            CALL SYJETF1(IEXAM,NJET,NPAR0,PPAR0,IJJ,FRAC,IJ ,PJ6)
            CALL SYJETF1(IEXAM,NJET,NPAR0,PPAR0,IJ ,FRAC,IJJ,PJ6)
            CALL SYJETF1(IEXAM,NJET,NPAR0,PPAR0,IJJ,FRAC,IJ ,PJ6)
         ELSE  ! 21
*          3 times iteration
            CALL SYJETF2(IEXAM,NJET,NPAR0,PPAR0,IJJ,FRAC,IJ ,PJ6)
            CALL SYJETF2(IEXAM,NJET,NPAR0,PPAR0,IJ ,FRAC,IJJ,PJ6)
            CALL SYJETF2(IEXAM,NJET,NPAR0,PPAR0,IJJ,FRAC,IJ ,PJ6)
         ENDIF
      ELSE
         WRITE(*,'(A,I10)')'unknown ISYMOD',ISYMOD
*         CALL REPORT('SYJKRN',22,'E')
         IERR = -9999
         GOTO 999
      ENDIF
*=========== Energy Ordering ================
      DO I=1,NPAR0
         IJJ(I)=IJ(I)
      ENDDO
      DO I=1,NJET
         PORD(I) = PJ6(4,I)
      ENDDO
      DO I=1,NJET
         RMAX = -1.E+09
         II=0
         DO J=1,NJET
            IF(PORD(J).GT.RMAX) THEN
               RMAX = PORD(J)
               II = J
            ENDIF
         ENDDO
         IIII(I)  = II
         PORD(II) = -1.E+09
      ENDDO
      DO I=1,NPAR0
         IASS(I) = 0
         DO J=1,NJET
            IF(IIII(J).EQ.IJJ(I)) THEN
               IASS(I) = J
            ENDIF
         ENDDO
      ENDDO
      DO I=1,MIN(NJET,MXJET)
         DO J=1,MIN(6,JDIM)
            PJETS(J,I)=PJ6(J,IIII(I))
         ENDDO
      ENDDO

      IERR = 0
*---- ALL DONE : GOOD
*      WRITE(*,'(A,I10,I10)')'complete',NJET,ISYMOD
*      CALL REPORT('SYJKRN',23,'I')

      RETURN

 999  CONTINUE
*----- error, so all initialyze again.
      WRITE(*,'(A,I10,I10,I6)')'Error Found',IERR,NJET,ISYMOD
*      CALL REPORT('SYJKRN',24,'E')
      NJET = 0
      DO I=1,NPAR
         IASS(I)=0
      ENDDO
      CALL VZERO(PJETS,MXJET*JDIM)
      RETURN
      END

*********************************************************************
      SUBROUTINE SYJCP6(NPAR,PPAR,NJ,IASS,JDIM,PJ6)
* Just Re-calculate Jet Momenta
*  Author ; Satoru Yamashita
*  Creat  ; 10-May-98
*********************************************************************
      implicit none
      INTEGER NPAR,IASS(*),JDIM,NJ
      REAL    PPAR(5,*),PJ6(JDIM,*),RR,QQ
      INTEGER I
      CALL VZERO(PJ6,NJ*JDIM)
      DO I=1,NPAR
         IF(IASS(I).GT.0.AND.IASS(I).LE.NJ) THEN
            PJ6(1,IASS(I)) = PJ6(1,IASS(I)) + PPAR(1,I)
            PJ6(2,IASS(I)) = PJ6(2,IASS(I)) + PPAR(2,I)
            PJ6(3,IASS(I)) = PJ6(3,IASS(I)) + PPAR(3,I)
            PJ6(4,IASS(I)) = PJ6(4,IASS(I)) + PPAR(4,I)
         ENDIF
      ENDDO
      IF(JDIM.GE.5) THEN
         DO I=1,NJ
            RR = PJ6(1,I)**2+PJ6(2,I)**2+PJ6(3,I)**2
            QQ = PJ6(4,I)**2-RR
            PJ6(5,I)=SQRT(MAX(QQ,0.))
            IF(JDIM.GE.6) PJ6(6,I) = SQRT(RR)
         ENDDO
      ENDIF
      END
