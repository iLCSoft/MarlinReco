*  
* $Id: trempt.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:35  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      SUBROUTINE TREMPT(PBEST,LOUT,LIN,NIN)
C****************************************************************************
c* Version number: 0.0      Date of last modification  26.08.99             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
C* Main modifiers from ALEPH to BRAHMS code:                                *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
C-----*************************************--------------------------
C! REMOVE BAD TPC POINTS TO GET A TRACK WITH BETTER CHI**2
C!
C!    AUTHOR        :- MARTIN POPPE 87/9/25
C!
C!
C!    INPUT  : LIN = ARRAY OF 'TPCO' COORDINATE NUMBERS OF A TRACK
C!             NIN = NUMBER OF INPUT PARTICLES
C!    OUTPUT : LOUT= ARRAY OF COORDINATE NUMBERS WITH 1 POINT
C!                   REMOVED SUCH THAT CHI**2 IS MINIMISED
C!             PBEST=FIT PROBABILITY OF BEST CASE
C?
C?    DESCRIPTION
C?    ===========
C?
C?    EACH OF THE POINTS IN "LIN" IS THROWN AWAY ONCE
C?    AT A TIME. THE REMAINING SET OF POINTS IS FITTED.
C?    LOUT IS FILLED WITH THOSE (NIN-1) POINTS FOR WHICH
C?    CHI**2 IS MINIMAL. THE INPUT ARGUMENTS REMAIN UNCHANGED
C------------------------------------------------------------------
C My stuff:

      INTEGER MPT
      PARAMETER(MPT=500)
      DOUBLE PRECISION XF(MPT),YF(MPT),WF(MPT)
      REAL RF(MPT),PF(MPT),ZF(MPT),WZF(MPT)

#include "include/fndparms.inc"
#include "include/tpco.inc"
C---------------------------------------------------------------------
      LOGICAL LFLIP
c GAB:
C      DIMENSION ERR(15),VAL(5),LIN(*),LOUT(*),LEVEN(12),LODD(12)
      DIMENSION ERR(15),VAL(5),LIN(*),LOUT(*),LEVEN(MPT),LODD(MPT)

      SAVE

      IF(NIN.GT.10) GO TO 100
C
C     ONLY FOR SHORT TRACKS, ALL POINTS ARE INVESTIGATED
C
      NM1 = NIN-1
      IBEST = 1

c GAB:
c      CALL TFTHEL
c     1  ('TPCO',NM1,LIN(2),1,VAL,ERR,CBEST)
           DO I=1,NM1
            XF(I)=DBLE(TPCO_XV(LIN(I+1)))
            YF(I)=DBLE(TPCO_YV(LIN(I+1)))
            RF(I)=TPCO_RV(LIN(I+1))
            PF(I)=TPCO_PH(LIN(I+1))
            WF(I)=1.D0/(DBLE(TPCO_SR(LIN(I+1)))+1.D-9)
            ZF(I)=TPCO_ZV(LIN(I+1))
            WZF(I)=1./(TPCO_SZ(LIN(I+1))+1.E-6)
           ENDDO

           CALL TFITHL(NM1,XF,YF,RF,PF,WF,ZF,WZF,1,
     1                 VAL,ERR,CHI2,CHISZ)

           CBEST=CHI2+CHISZ

C GAB:
      DO 10 I=2,NIN
        IDUMY = LIN(I)
        LIN(I)= LIN(1)
        LIN(1)= IDUMY
C        CALL TFTHEL
C     1  ('TPCO',NM1,LIN(2),1,VAL,ERR,CHI2)
           DO J=2,NM1
            XF(J)=DBLE(TPCO_XV(LIN(J)))
            YF(J)=DBLE(TPCO_YV(LIN(J)))
            RF(J)=TPCO_RV(LIN(J))
            PF(J)=TPCO_PH(LIN(J))
            WF(J)=1.D0/(DBLE(TPCO_SR(LIN(J)))+1.D-9)
            ZF(J)=TPCO_ZV(LIN(J))
            WZF(J)=1./(TPCO_SZ(LIN(J))+1.E-6)
           ENDDO
           CALL TFITHL(NM1,XF,YF,RF,PF,WF,ZF,WZF,1,
     1                 VAL,ERR,CHI2,CHISZ)
           CHI2=CHI2+CHISZ

        IF(CHI2.LT.CBEST)THEN
             IBEST = I
             CBEST = CHI2
        END IF
   10 CONTINUE
C
C-----> RESTORE LIN IN ITS ORIGINAL FORM
C
      IDUMY = LIN(1)
      DO 20 I=2,NIN
   20   LIN(I-1)=LIN(I)
      LIN(NIN)=IDUMY
C
C-----> LIN -> LOUT
C
      IF(IBEST.NE.1) THEN
        DO 30 I=1,IBEST-1
   30   LOUT(I) = LIN(I)
      END IF
      IF(IBEST.NE.NIN) THEN
        DO 40 I=IBEST,NM1
   40   LOUT(I) = LIN(I+1)
      END IF
c MVL      PBEST = PROB(CBEST,(2*NM1-5))
      IF(CBEST.GT.0.0.and.CBEST.lT.300.0.and.2*NM1-5.GT.0)THEN
        PBEST=max(PROB(CBEST,(2*NM1-5)),1.e-25)
      ELSE
        PBEST=1.e-25
      ENDIF
c MVL
      RETURN
  100 CONTINUE
C
C-----> LONG TRACKS ARE SPLIT IN EVEN AND ODD PARTS
C
      LFLIP = .TRUE.
      DO 110 I=1,NIN
      IF(LFLIP) THEN
        K=I/2 + 1
        LODD(K)= LIN(I)
      ELSE
        K = I/2
        LEVEN(K)=LIN(I)
      END IF
      LFLIP = .NOT.LFLIP
  110 CONTINUE
      NEVEN = NIN/2
      NODD  = NIN-NEVEN
C GAB:
C      CALL TFTHEL
C     1  ('TPCO',NEVEN,LEVEN,1,VAL,ERR,CEVEN)
           DO I=1,NEVEN
            XF(I)=DBLE(TPCO_XV(LEVEN(I)))
            YF(I)=DBLE(TPCO_YV(LEVEN(I)))
            RF(I)=TPCO_RV(LEVEN(I))
            PF(I)=TPCO_PH(LEVEN(I))
            WF(I)=1.D0/(DBLE(TPCO_SR(LEVEN(I)))+1.D-9)
            ZF(I)=TPCO_ZV(LEVEN(I))
            WZF(I)=1./(TPCO_SZ(LEVEN(I))+1.E-6)
           ENDDO

           CALL TFITHL(NEVEN,XF,YF,RF,PF,WF,ZF,WZF,1,
     1                 VAL,ERR,CHI2,CHISZ)
           CEVEN=CHI2+CHISZ

C GAB:
C      CALL TFTHEL
C     1  ('TPCO',NODD,LODD,1,VAL,ERR,CODD )
           DO I=1,NODD
            XF(I)=DBLE(TPCO_XV(LODD(I)))
            YF(I)=DBLE(TPCO_YV(LODD(I)))
            RF(I)=TPCO_RV(LODD(I))
            PF(I)=TPCO_PH(LODD(I))
            WF(I)=1.D0/(DBLE(TPCO_SR(LODD(I)))+1.D-9)
            ZF(I)=TPCO_ZV(LODD(I))
            WZF(I)=1./(TPCO_SZ(LODD(I))+1.E-6)
           ENDDO

           CALL TFITHL(NODD,XF,YF,RF,PF,WF,ZF,WZF,1,
     1                 VAL,ERR,CHI2,CHISZ)
           CODD=CHI2+CHISZ

      C0  = 0.0
      CEVEN = CEVEN/FLOAT(NEVEN)
      CODD  = CODD /FLOAT(NODD)
      IF(CEVEN.GT.CODD)THEN
         DO 140 I =1,NEVEN
           LODD(NODD+1) = LEVEN(I)
C GAB:
C           CALL TFTHEL
C     1     ('TPCO',(NODD+1),LODD,1,VAL,ERR,CC)
           DO J=1,NODD+1
            XF(J)=DBLE(TPCO_XV(LODD(J)))
            YF(J)=DBLE(TPCO_YV(LODD(J)))
            RF(J)=TPCO_RV(LODD(J))
            PF(J)=TPCO_PH(LODD(J))
            WF(J)=1.D0/(DBLE(TPCO_SR(LODD(J)))+1.D-9)
            ZF(J)=TPCO_ZV(LODD(J))
            WZF(J)=1./(TPCO_SZ(LODD(J))+1.E-6)
           ENDDO

           CALL TFITHL(NODD+1,XF,YF,RF,PF,WF,ZF,WZF,1,
     1                 VAL,ERR,CHI2,CHISZ)
           CC=CHI2+CHISZ

           IF(CC.GT.C0) THEN
                 KK = I
                 C0 = CC
           END IF
  140   CONTINUE
        KSHIT = 2*KK
      ELSE
         DO 150 I =1,NODD
           LEVEN(NEVEN+1) = LODD(I)
C GAB:
C           CALL TFTHEL
C     1     ('TPCO',(NEVEN+1),LEVEN,1,VAL,ERR,CC)
           DO J=1,NEVEN+1
            XF(J)=DBLE(TPCO_XV(LEVEN(J)))
            YF(J)=DBLE(TPCO_YV(LEVEN(J)))
            RF(J)=TPCO_RV(LEVEN(J))
            PF(J)=TPCO_PH(LEVEN(J))
            WF(J)=1.D0/(DBLE(TPCO_SR(LEVEN(J)))+1.D-9)
            ZF(J)=TPCO_ZV(LEVEN(J))
            WZF(J)=1./(TPCO_SZ(LEVEN(J))+1.E-6)
           ENDDO

           CALL TFITHL(NEVEN+1,XF,YF,RF,PF,WF,ZF,WZF,1,
     1                 VAL,ERR,CHI2,CHISZ)
           C0=CHI2+CHISZ
           IF(CC.GT.C0) THEN
                 KK = I
                 C0 = CC
           END IF
        KSHIT = 2*KK+1
  150   CONTINUE
      END IF
      DO 160 I=1,KSHIT
  160   LOUT(I) = LIN(I)
      NM1 = NIN-1
      IF(KSHIT.LT.NIN)THEN
        DO 170 I=KSHIT,NM1
  170     LOUT(I) = LIN(I+1)
      END IF
C GAB:
C        CALL TFTHEL
C     1  ('TPCO',NM1,LOUT(1),1,VAL,ERR,CBEST)
           DO I=1,NM1
            XF(I)=DBLE(TPCO_XV(LOUT(I)))
            YF(I)=DBLE(TPCO_YV(LOUT(I)))
            RF(I)=TPCO_RV(LOUT(I))
            PF(I)=TPCO_PH(LOUT(I))
            WF(I)=1.D0/(DBLE(TPCO_SR(LOUT(I)))+1.D-9)
            ZF(I)=TPCO_ZV(LOUT(I))
            WZF(I)=1./(TPCO_SZ(LOUT(I))+1.E-6)
           ENDDO
           CALL TFITHL(NM1,XF,YF,RF,PF,WF,ZF,WZF,1,
     1                 VAL,ERR,CHI2,CHISZ)
           CBEST=CHI2+CHISZ

c MVL      PBEST = PROB(CBEST,(2*NM1-5))
      IF(CBEST.GT.0.0.and.CBEST.lT.300.0.and.2*NM1-5.GT.0)THEN
        PBEST=max(PROB(CBEST,(2*NM1-5)),1.e-25)
      ELSE
        PBEST=1.e-25
      ENDIF
c MVL
      RETURN
      END

c------------------------------------------------------------------------------
