*  
* $Id: tkink.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      SUBROUTINE TKINK(DCHI2,KINK,LIN,NIN)
C****************************************************************************
c* Version number: 0.0      Date of last modification  26.08.99             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
C* Main modifiers from ALEPH to BRAHMS code:                                *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
C************************************--------------------------
C!    FIND KINK IN A TPC TRACK
C!
C!    AUTHOR        :- GEORG STIMPFL 87/9/25
C!
C!    DO NOT CALL THIS SUBROUTINE FOR LESS THAN 10 POINTS
C!    ===================================================
C!
C!    INPUT  : LIN = ARRAY OF 'TPCO' COORDINATE NUMBERS OF A TRACK
C!             NIN = NUMBER OF INPUT PARTICLES
C!    OUTPUT : KINK= FIRST COORDINATE AFTER KINK
C!                   REMOVED SUCH THAT CHI**2 IS MINIMISED
C!             DCHI2=IMPROVEMENT OF CHI**2 BY PUTTING THE KINK INTO
C!                   THE RIGHT PLACE
C?
C?    DESCRIPTION
C?    ===========
C?
C?    THE KINK POINT IS CALCULATED AS THAT ONE FOR WHICH
C?    THE SUM OF THE CHI**2 VALUES FOR BOTH ENDS OF THE TRACK
C?    IS MINIMAL. DCHI2 IS DEFINED AS
C?    DCHI2 = CHI**2(BEST KINK POINT) - CHI**2(WORST KINK POINT)
C------------------------------------------------------------------
C
C GAB:
      INTEGER MPT
      PARAMETER(MPT=500)
      DOUBLE PRECISION XF(MPT),YF(MPT),WF(MPT)
      REAL RF(MPT),PF(MPT),ZF(MPT),WZF(MPT)

      DIMENSION ERR(15),VAL(5),LIN(*)
#include "include/fndparms.inc"
#include "include/tpco.inc"

      SAVE

      SC1= 0.0
      SC2= 1.0E30
      IFIRST = 4
      ILAST  = NIN-4
      DO 10 I=IFIRST,ILAST
        NP1 = I
        NP2 = NIN-NP1
C        CALL TFTHEL
C     1  ('TPCO',NP1,LIN(1),1,VAL,ERR,C1)
C GAB:
        DO J=1,NP1
         XF(J)=DBLE(TPCO_XV(LIN(J)))
         YF(J)=DBLE(TPCO_YV(LIN(J)))
         RF(J)=TPCO_RV(LIN(J))
         PF(J)=TPCO_PH(LIN(J))
         WF(J)=1.D0/(DBLE(TPCO_SR(LIN(J)))+1.D-9)
         ZF(J)=TPCO_ZV(LIN(J))
         WZF(J)=1./(TPCO_SZ(LIN(J))+1.E-6)
        ENDDO

        CALL TFITHL(NP1,XF,YF,RF,PF,WF,ZF,WZF,1,
     1                 VAL,ERR,CHI2,CHISZ)

         C1=CHI2+CHISZ

C GAB:
C        CALL TFTHEL
C     1  ('TPCO',NP2,LIN(I),1,VAL,ERR,C2)
        DO J=1,NP2
         XF(J)=DBLE(TPCO_XV(LIN(I+J-1)))
         YF(J)=DBLE(TPCO_YV(LIN(I+J-1)))
         RF(J)=TPCO_RV(LIN(I+J-1))
         PF(J)=TPCO_PH(LIN(I+J-1))
         WF(J)=1.D0/(DBLE(TPCO_SR(LIN(I+J-1)))+1.D-9)
         ZF(J)=TPCO_ZV(LIN(I+J-1))
         WZF(J)=1./(TPCO_SZ(LIN(I+J-1))+1.E-6)
        ENDDO

        CALL TFITHL(NP2,XF,YF,RF,PF,WF,ZF,WZF,1,
     1                 VAL,ERR,CHI2,CHISZ)

        C2=CHI2+CHISZ
        SC = C1+C2
        IF(SC.GT.SC1) SC1=SC
        IF(SC.LE.SC2) THEN
          SC2 = SC
          KINK= I + 1
        END IF
   10 CONTINUE
      DCHI2 = SC1 - SC2
      RETURN
      END
C------------------------------------------------------------------
