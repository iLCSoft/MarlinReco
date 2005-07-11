CDECK  ID>, SYJETF0.
      SUBROUTINE SYJETF0(IMODE,YCUT,EPS0,R0,NJETRQ,NPAR0,PPAR0,IUSE,
     +                   NJET,IPASS0,PJET6,Y34,Y45,IERR)
*  Traditional Jet Finders
*  Authour : Satoru Yamashita
*  Creat   : 1-Oct-97
*
*     IMODE : 1,2,3,4,5,6,7,8,9,10,11,12   : fixed YCUT, or EPS-R cut
*     IMODE : 101,102,103,104,105,,,,112   : fixed NJET
*
      implicit none
      INTEGER IMODE
      REAL    YCUT,RIN,EPSIN
      INTEGER NJETRQ,NPAR0,IPASS0(*),IUSE(*),JMODE
      INTEGER NPAR,IERR,I,NJET,IPASS(300),K,IJMUL(10),J
      REAL PPAR(5,300),PPAR0(5,*)
      REAL PJJ(10,10),Y34,Y45,PJET(5,10),PJET6(6,*)
      REAL P4(4,10)
      REAL EPS0,R0
      REAL ROUT,EPSOUT,YDUM
      INTEGER NJETOUT,NJDUM

*====== initialyzation
      NJET = 0
      DO I=1,NPAR0
         IPASS0(I)=-1
      ENDDO
      CALL VZERO(PJET,5*10)
      CALL VZERO(IPASS0,NPAR0)
      Y34 = 0.
      Y45 = 0.
      IERR = 0
      NPAR = 0
      DO I=1,NPAR0
         IF(IUSE(I).GT.0) THEN
            NPAR = NPAR + 1
            DO K=1,5
               PPAR(K,NPAR) = PPAR0(K,I)
            ENDDO
         ENDIF
      ENDDO
      IF(NPAR.LE.2) THEN
         WRITE(*,'(A,I3,I10)')'NPAR LT 3',NPAR,NPAR0
*         CALL REPORT('JETFIN',1,'E')
         IERR = -10
         GOTO 999
      ENDIF
      IF(IMODE.GT.100.AND.NPAR.LT.NJETRQ) THEN
         WRITE(*,'(A,I3,I10)')'NPAR LT NJRQ',NPAR,NJETRQ
*         CALL REPORT('JETFIN',2,'E')
         IERR = -5
         GOTO 999
      ENDIF


*======  normal YKERN with fixed YCUT
      IF(IMODE.GE.1.AND.IMODE.LE.12) THEN

         IF(IMODE.GE.1.AND.IMODE.LE.6) THEN
            JMODE=IMODE
            CALL YKERN(JMODE,NPAR,5,PPAR,IERR)
            IF(IERR.NE.0) THEN
               WRITE(*,'(A,I3,I10)')'YKERN ERROR 1',I,IERR
*               CALL REPORT('JETFIN',3,'E')
               GOTO 999
            ENDIF
            CALL YNJET(YCUT,NJET,IERR) ! using input YCUT
            IF(IERR.NE.0) THEN
               WRITE(*,'(A,I3,I10)')'YNJET ERROR 2',I,IERR
*               CALL REPORT('JETFIN',4,'E')
               GOTO 999
            ENDIF
            IF(NJET.GT.10) THEN
               WRITE(*,'(A,I3,I10)')
     +              'NJET EXCEED 10 USE 10',NJET,IMODE
*               CALL REPORT('JETFIN',5,'W')
               NJET = 10
            ENDIF
            CALL YASSO(NJET,PJJ,IPASS,IERR)
            IF(IERR.NE.0) THEN
               WRITE(*,'(A,I3,I10)')'YASSO ERROR 1',I,IERR
*               CALL REPORT('JETFIN',6,'E')
               GOTO 999
            ENDIF
            DO K=1,NJET
               CALL UCOPY(PJJ(1,K),PJET(1,K),4)
            ENDDO

         ELSEIF(IMODE.EQ.7) THEN ! Cambridge
            CALL PXCAMJ(5,NPAR,PPAR,YCUT,NJET,IPASS,P4,IERR)
            IF(IERR.NE.0) THEN
               WRITE(*,'(A,I3,I10)')'PXCAMJ ERROR ',NJET,IERR
*               CALL REPORT('JETFIN',7,'E')
               GOTO 999
            ENDIF
            DO K=1,NJET
               CALL UCOPY(P4(1,K),PJET(1,K),4)
               PJET(5,K) = SQRT(MAX(0.,
     +           PJET(4,K)**2-PJET(1,K)**2-PJET(2,K)**2-PJET(3,K)**2))
            ENDDO

         ELSEIF(IMODE.GE.8.AND.IMODE.LE.11) THEN ! LUCLUS default
C     ! YCUT==XMIN
            JMODE = IMODE - 7
            CALL PXLCL5(NPAR,5,PPAR,JMODE,YCUT,1,10,
     +           NJET,PJET,IPASS,IJMUL,IERR)
            IF(IERR.NE.0) THEN
               WRITE(*,'(A,I3,I10)')'PXLCL5 ERROR ',NJET,IERR
*               CALL REPORT('JETFIN',8,'E')
               GOTO 999
            ENDIF
            DO K=1,NJET
               PJET(5,K) = SQRT(MAX(0.,
     +           PJET(4,K)**2-PJET(1,K)**2-PJET(2,K)**2-PJET(3,K)**2))
            ENDDO

         ELSEIF(IMODE.EQ.12) THEN ! CONE
            RIN   = R0
            EPSIN = EPS0
            CALL PXCONE(NPAR,5,PPAR,RIN,EPSIN,10,
     +           NJET,PJET,IPASS,IJMUL,IERR)
            IF(IERR.NE.0) THEN
               WRITE(*,'(A,I3,I10)')'PXCONE ERROR ',NJET,IERR
*               CALL REPORT('JETFIN',9,'E')
               GOTO 999
            ENDIF
            DO K=1,NJET
               PJET(5,K) = SQRT(MAX(0.,
     +           PJET(4,K)**2-PJET(1,K)**2-PJET(2,K)**2-PJET(3,K)**2))
            ENDDO

         ENDIF

      ELSEIF(IMODE.GE.101.AND.IMODE.LE.112) THEN
         NJET = NJETRQ

         IF(IMODE.GE.101.AND.IMODE.LE.106) THEN
            JMODE=IMODE-100
            CALL YKERN(JMODE,NPAR,5,PPAR,IERR) ! DURHAM
            IF(IERR.NE.0) THEN
               WRITE(*,'(A,I5,I10)')'YKERN ERROR 2',IMODE,IERR
*               CALL REPORT('JETFIN',10,'E')
               GOTO 999
            ENDIF
            CALL YYJET(NJET,Y45,Y34,IERR)
            IF(IERR.NE.0) THEN
               WRITE(*,'(A,I3,I10)')'YYJET ERROR 1',I,IERR
*               CALL REPORT('JETFIN',11,'E')
               GOTO 999
            ENDIF

*           temp
C            write(*,*)'NJET=',NJET,NPAR0,NPAR,IMODE,JMODE,PPAR(4,NPAR)

            CALL YASSO(NJET,PJJ,IPASS,IERR)
C            write(*,*)'over YASSO'

            IF(IERR.NE.0) THEN
               WRITE(*,'(A,I3,I10)')'YASSO ERROR 1',I,IERR
*               CALL REPORT('JETFIN',12,'E')
               GOTO 999
            ENDIF
            DO K=1,NJET
               CALL UCOPY(PJJ(1,K),PJET(1,K),4)
               PJET(5,K) = SQRT(MAX(0.,
     +           PJET(4,K)**2-PJET(1,K)**2-PJET(2,K)**2-PJET(3,K)**2))
            ENDDO

         ELSEIF(IMODE.EQ.107) THEN ! CAMBRIDGE WITH FIXED NJET
            CALL CKERN(5,NPAR,PPAR,NJET,IERR)
            IF(IERR.NE.0) THEN
               WRITE(*,'(A,I5,I10)')'CKERN ERROR',IMODE,IERR
*               CALL REPORT('JETFIN',13,'E')
               GOTO 999
            ENDIF
            CALL CYJET(NJET,Y45,Y34,IERR)
            IF(IERR.NE.0) THEN
               WRITE(*,'(A,I3,I10)')'CYJET ERROR',I,IERR
*               CALL REPORT('JETFIN',14,'E')
               GOTO 999
            ENDIF
            CALL CASSO(NJET,P4,IPASS,IERR)
            IF(IERR.NE.0) THEN
               WRITE(*,'(A,I3,I10)')'CASSO ERROR 1',I,IERR
*               CALL REPORT('JETFIN',15,'E')
               GOTO 999
            ENDIF
            DO K=1,NJET
               CALL UCOPY(P4(1,K),PJET(1,K),4)
               PJET(5,K) = SQRT(MAX(0.,
     +           PJET(4,K)**2-PJET(1,K)**2-PJET(2,K)**2-PJET(3,K)**2))
            ENDDO

         ELSEIF(IMODE.GE.108.AND.IMODE.LE.111) THEN ! LUCLUS NOT YET
            JMODE = IMODE - 107
            YDUM = 50.0
            IF(IMODE.EQ.111) YDUM = 0.99
            NJDUM = NJET
            CALL PXLCL5(NPAR,5,PPAR,JMODE,YDUM,NJDUM,10,
     +           NJET,PJET,IPASS,IJMUL,IERR)
            IF(IERR.NE.0) THEN
               WRITE(*,'(A,I3,I10)')'PXLCL5 ERROR ',NJET,IERR
*               CALL REPORT('JETFIN',8,'E')
               GOTO 999
            ENDIF
            DO K=1,NJET
               PJET(5,K) = SQRT(MAX(0.,
     +          PJET(4,K)**2-PJET(1,K)**2-PJET(2,K)**2-PJET(3,K)**2))
            ENDDO

         ELSE IF(IMODE.EQ.112) THEN
            EPSOUT = 10.0       ! EPS = 10 GEV
            CALL SYCONE(NJET,NPAR,5,PPAR,ROUT,EPSOUT,
     +           10,NJETOUT,PJET,IPASS,IJMUL,IERR)
            IF(IERR.NE.0) THEN
               WRITE(*,'(A,I3,I10)')'SYCONE ERROR',NJET,IERR
*               CALL REPORT('JETFIN',17,'E')
               GOTO 999
            ENDIF
            IF(NJET.NE.NJETOUT) THEN
               WRITE(*,'(A,I3,I10)')
     +              'SYCONE NJ NOT MATCH',NJET,NJETOUT
*               CALL REPORT('JETFIN',18,'E')
               IERR = 1
               GOTO 999
            ENDIF
            DO J=1,NPAR
               IF(IPASS(J).LE.0) THEN
                  IPASS(J)=0
               ENDIF
            ENDDO
            Y34 = ROUT
            Y45 = EPSOUT
         ENDIF
      ELSE
         WRITE(*,'(A,I5,I10)')'NO JET ALGORITHM',IMODE,IERR
*         CALL REPORT('JETFIN',19,'E')
         IERR = -999
         GOTO 999
      ENDIF

      DO I=1,NPAR
         IF(IPASS(I).GT.NJET.OR.IPASS(I).LT.1) THEN
            IPASS(I) = 0
         ENDIF
      ENDDO
      K = 0
      DO I=1,NPAR0
         IF(IUSE(I).GT.0) THEN
            K = K + 1
            IPASS0(I) = IPASS(K)
         ENDIF
      ENDDO
      DO K=1,NJET
         PJET(5,K) = SQRT(MAX(0.,
     +        PJET(4,K)**2-PJET(1,K)**2-PJET(2,K)**2-PJET(3,K)**2))
      ENDDO
      DO K=1,NJET
         CALL UCOPY(PJET(1,K),PJET6(1,K),5)
         PJET6(6,K) = SQRT(PJET(1,K)**2+PJET(2,K)**2+PJET(3,K)**2)
      ENDDO

      RETURN

 999  CONTINUE
      write(*,*)' BAD OK1',NJET,NPAR0,IMODE,IERR,NJETRQ
      NJET = 0
      CALL VZERO(PJET6,6*10)
      CALL VZERO(IPASS0,NPAR0)
      Y34 = 0.
      Y45 = 0.
C      write(*,*)' BAD OK2'

      RETURN

      END
