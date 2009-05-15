*  
* $Id: tfittk.F,v 1.4 2009-05-15 15:12:29 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.3  2008/06/26 09:35:56  aplin
*
* A clean up has been made of the arrays that hold the TE's/Hits for
* track fitting. An effort has been made to ensure that both the arrays
* are big enough to hold the required number of hits, and that the arrays
* are all of the same dimensions between calling routines.
*
* All verbose output, except for the case where the tracking is dropped
* for a complete event, has now been hidden behind debug IF
* Statements. The level of verbosity can be controlled via the variables
* TPCDBG, TKSTDBG, FSFDBG and IDEB  in tkinit.F
*
* Revision 1.2  2005/08/08 07:09:13  aplin
* Made f77 tracking code use GEAR to define the geomtery of the TPC. LTPDRO now defines the maximum number of rows is used to define the size of arrays, this is limited to 224 due the use of 7 '32 bit' bit registers in trkfnd.F increased, though at present it is not likely that anybody would want more. The number of TPC padrows is defined at run time by NRTPC which should of course not exceed LTPDRO, although this is checked and the programe exits with a verbose error message. A wrapper function gettpcgeom is used to pass the GEAR TPC parameters from C++ to f77. MarlinUtil/include/marlin_tpcgeom.h have MarlinUtil/src/marlin_tpcgeom.cc consequently been removed as they are no longer needed.
*
* Revision 1.1.1.1  2005/07/11 15:15:48  aplin
* Initial version
*
* Revision 1.2  2005/06/20 18:49:46  aplin
* change tpc geometery so that it is read from  same header file as cpp i.e. marlin_tpcgeom.h
*
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
      SUBROUTINE TFITTK(IER)
C****************************************************************************
c* Version number: 0.0      Date of last modification  26.08.99             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
C* Main modifiers from ALEPH to BRAHMS code:                                *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
C-----**********************----------------------------------------
      IMPLICIT NONE
C! Fit TPC track candidates and produce TPC tracks
C!
C!    AUTHOR        :- MARTIN POPPE 87/9/25
C!    Modified      :  I.Tomalin    12/3/94
C!                     Protection put in requiring chains to have at
C!                     least 3 coords. (Necessary because TOVRLP
C!                     sometimes leaves chains with less than this).
C!                  :  D.Casper     09/5/97
C!                     Protect against crazy fits by calling sanity check
C!
C!    INPUT  BANKS  : 'TCAL' , 'TTCC' , 'TPCO'
C!    OUTPUT BANKS  : 'T1FT' , 'T1TL' , 'T1CL'
C!
C?    DESCRIPTION
C?    ===========
C?
C?    THE HELIX PARAMETERS OF ALL TRACK CANDIDATES ARE FITTED.
C?    TRACKS WITH A BAD CHI**2 PROBABILITY ARE MASSAGED:
C?
C?    FIRST, IT IS TRIED TO IMPROVE THE FIT BY EXCLUDING
C?    THE POINT WITH THE LARGEST CONTRIBUTION TO CHI**2.
C?
C?    IF THIS METHOD DOES NOT LEAD TO AN ACCEPTABLE TRACK,
C?    A KINK SEARCH "A LA STIMPFL" IS MADE: THE CRITERIUM FOR
C?    A KINK IS THAT THE DIFFERENCE IN CHI**2 BETWEEN THE "BEST"
C?    KINK POINT AND THE "WORST" KINK POINT HYPOTHESIS IS
C?    LARGER THAN A CUT VALUE (CALLED DCM IN THIS PROGRAM).
C?    IF A KINK HAS BEEN FOUND, THE TRACK IS SPLIT INTO TWO PIECES.
C?    IF THE FIT PROBABILITY OF ONE OF THESE PIECES IS TOO LOW,
C?    THE POINT(S) WITH THE LARGEST CONTRIBUTION TO CHI**2
C?    ARE REMOVED.
C?
C?    FOR THE REMAINING TRACK CANDIDATES, POINTS
C?    OF A TRACK CANDIDATE ARE THROWN AWAY "ONE BY ONE" TO OBTAIN
C?    A BEST FIT.
C?
C?    IF NEITHER MASSAGE HELPS, THE TRACK CANDIDATE (OR A PIECE
C?    THEREOF IN CASE OF A KINK) IS THROWN AWAY.
C?
C-------------------------------------------------------------------
#include "include/padrow.inc"
#include "include/fndparms.inc"
#include "include/tpco.inc"
#include "include/tpgeom.inc"
c#include "marlin_tpcgeom.h"
#include "include/my_tcal.inc"
#include "include/my_ttcc.inc"
#include "include/my_t1ft.inc"
#include "include/my_t1tl.inc"
#include "include/my_t1cl.inc"
#include "include/my_tctc.inc"
c#include "gcflag.inc"
C(KH) +SEQ,TFTKWB.
#include "include/strcons.inc"


c----------------------------------------------------------------------
      REAL VV0(5),EE0(15),CHI2
      INTEGER MLIST
      PARAMETER (MLIST=LTPDRO)
      INTEGER LIST(MLIST)
C------------------------------------------------------------------------
      LOGICAL FIRST
      REAL    DCHI2(LTPDRO)
      INTEGER IBUFF(LTPDRO)
C(SJA) DCHI2 and IBUFF should be of the same dimention as the max number of tpcrows
C(SJA)      REAL    DCHI2(500)
C(SJA)      INTEGER IBUFF(500)

      DATA FIRST/.TRUE./

      INTEGER IWORST(MAXCND)
      INTEGER I,II,IER,NTKOK,NPTOK,ICAND,JT1FT,NPT,KSTCC,KSCAL,KSGTL
      INTEGER KSGCL,NLONG,IPT,LOST,NGET,IGET,IPUT,KINK,IACC1,IACC2
      INTEGER NPT1,NPT2,NCOR,NTK,ITK,JT1TL,I1,I2,I3,ICO
      INTEGER KKTTCC,N1CORD,KKT1CL,IERR

      REAL PPP,DCF,PBEST

      REAL PROB
      EXTERNAL PROB


cC-----> OBTAIN LIST AND NUMBER OF POINTS FOR TRACK CANDIDATE 'I'
c GAB:
C      KKTTCC(II) = KTTCC+LMHLEN+ITABL(KTCAL,II,JTCAIO)
      KKTTCC(II) = MY_JTCAIO(II)

C GAB:
c      N1CORD(II) = ITABL(KTCAL,II,JTCANF)
      N1CORD(II) = MY_JTCANF(II)
cC-----> OBTAIN LIST AND NUMBER OF POINTS FOR TRACK CANDIDATE 'I'
c GAB:
c      KKT1CL(II) = KT1CL+LMHLEN+ITABL(KT1TL,II,JT1TIO)
       KKT1CL(II) = MY_JT1TIO(II)
cC-----> POINTER BEFORE COORDINATE LIST IN 'T1CL'
cC
cC---------------------- end of statement functions ------------------
cC
c      IF(FIRST) THEN
C        IWORST=0
c        NT1FT = NAMIND('T1FT')
c        NT1TL = NAMIND('T1TL')
c        NT1CL = NAMIND('T1CL')
c        NTTRC = NAMIND('TTRC')
c        NTCAL = NAMIND('TCAL')
c        NTTCC = NAMIND('TTCC')
c        NTPCO = NAMIND('TPCO')
c        FIRST = .FALSE.
c      END IF


      IER = 0

       NT1FT=0
       NT1TL=0
cC
cC---> Check for necessary input banks
cC
c      KTTRC = IW(NTTRC)
c      KTCAL = IW(NTCAL)
c      KTTCC = IW(NTTCC)
c      KTPCO = IW(NTPCO)
c      IF(KTTRC.LE.0)                              GO TO 998
c      IF(KTCAL.LE.0)                              GO TO 998
c      IF(KTTCC.LE.0)                              GO TO 998
c      IF(KTPCO.LE.0)                              GO TO 998
cC
cC---> Make sure that output banks are not already present
cC
c      IF (IW(NT1FT).NE.0) CALL BDROP(IW,'T1FT')
c      IF (IW(NT1TL).NE.0) CALL BDROP(IW,'T1TL')
c      IF (IW(NT1CL).NE.0) CALL BDROP(IW,'T1CL')
cC
C---> Create workbank for array IWORST
C--->   IWORST = 0       -> Track candidate has been accepted
C       IWORST = IPT > 0 -> Ipt is the worst point on the track
C                             (track not yet accepted)
c      NCAND = LROWS(KTCAL)

C GAB:
      IF(NCAND.GT.MAXCND)THEN
       WRITE(*,*)'TFITTK: MAXCND Too small.  Tracks may be lost.'
       CALL TKVETO(1)
      ENDIF

      IF(NCAND.LE.0)                                GO TO 1000
c      IW(1)=1
c      CALL WBANK(IW,IWORST,NCAND,*995)
cC
cC---> Creation of output banks
cC---> Create banks of size such that all particles are kinks
cC       (Twice as large as no kink, -> always large enough)
cC
c      LT1FT = LT1FTA*2*NCAND + LMHLEN
c      LT1TL = LT1TLA*2*NCAND + LMHLEN
c      LT1CL = IW(KTTCC)
c      CALL AUBOS('T1FT',0,LT1FT,KT1FT,IGARB)
c      IF (IGARB.EQ.2) GO TO 999
c      IF (IGARB.EQ.1) THEN
c        KTTRC = IW(NTTRC)
c        KTCAL = IW(NTCAL)
c        KTTCC = IW(NTTCC)
c        KTPCO = IW(NTPCO)
c      END IF
c      CALL AUBOS('T1TL',0,LT1TL,KT1TL,IGARB)
c      IF (IGARB.EQ.2) GO TO 999
c      IF (IGARB.EQ.1) THEN
c        KT1FT = IW(NT1FT)
c        KTTRC = IW(NTTRC)
c        KTCAL = IW(NTCAL)
c        KTTCC = IW(NTTCC)
c        KTPCO = IW(NTPCO)
c      END IF
c      CALL AUBOS('T1CL',1,LT1CL,KT1CL,IGARB)
c      IF (IGARB.EQ.2) GO TO 999
c      IF (IGARB.EQ.1) THEN
c         KT1FT = IW(NT1FT)
c         KT1TL = IW(NT1TL)
c         KTTRC = IW(NTTRC)
c         KTCAL = IW(NTCAL)
c         KTTCC = IW(NTTCC)
c         KTPCO = IW(NTPCO)
c      END IF
cC
cC---> Define number of columns
cC
c      IW(KT1FT+LMHCOL) = LT1FTA
c      IW(KT1TL+LMHCOL) = LT1TLA
c      IW(KT1CL+LMHCOL) = 1
cC
C*******************************************************************C
C                                                                   C
C       SIMPLE TRACKS = TRACKS WITH DECENT CHI**2 BEFORE FUDGEING   C
C                                                                   C
C*******************************************************************C
C
C GAB:
      NTKOK = 0
      NPTOK = 0

c gab:
      DO 100 ICAND = 1,MIN(NCAND,MAXCND)
c        JT1FT = KROW(KT1FT,NTKOK+1)
c gab:
        JT1FT=NTKOK+1
        NT1FT=JT1FT

        NPT   = N1CORD(ICAND)
C Protection against chains with less than three coords.
        IF (NPT.LT.3) THEN
c          IW(IWORST+ICAND)=0
          IWORST(ICAND)=0
          GOTO 100
        END IF
        KSTCC = KKTTCC(ICAND)
c(kh)
        IF (NPT.GT.MLIST) THEN
           WRITE(6,*) 'TFITTK: MLIST TOO SMALL. NEED ',NPT
           CALL TKVETO(1)
           IER=10
           RETURN
        END IF
c GAB:
        DO I=1,NPT
           LIST(I)= MY_JTCCIT(MY_JTCAIO(ICAND)+I)
        ENDDO


        
c     GAB:
c     CALL TMSHEL('TPCO',NPT,IW(KSTCC+1),2,RW(JT1FT+JT1FIR),
c     &                           RW(JT1FT+JT1FEM),RW(JT1FT+JT1FCD))

        CALL MY_TMSHEL(NPT,LIST,2,VV0,EE0,CHI2)
        
        MY_JT1FIR(JT1FT)=VV0(1)
        MY_JT1FTL(JT1FT)=VV0(2)

        DO I=1,15
           MY_JT1FEM(JT1FT,I)=EE0(I)
        ENDDO

        MY_JT1FCD(JT1FT)=CHI2
        
c GAB:
c     IW(JT1FT+JT1FDF) = 2*NPT-5
c     PPP = PROB(RW(JT1FT+JT1FCD),(2*NPT-5))
        
        MY_JT1FDF(JT1FT)= 2*NPT-5

C SJA:Protection again zero values which have no meaning        
        IF (MY_JT1FCD(JT1FT).GT.0.AND.(2*NPT-5).GT.0.AND.NPT.GT.0) THEN
           PPP = PROB(MY_JT1FCD(JT1FT),(2*NPT-5))
        ELSE
           WRITE(*,*) 
     +          "tfittk: CHI2 or NDF sent to CERNLIB PROB() are zero"
           PPP = 1e-25
        ENDIF
        
C---- >  Protect against hopeless fits
C     GAB:
C     CALL TFCHEK(RW(JT1FT+JT1FIR),.TRUE.,IERR)

        CALL TFCHEK(VV0,.TRUE.,IERR)
        
        IF (PPP.GT.PMIN .AND. IERR.EQ.0) THEN
           NTKOK = NTKOK+1

C     GAB:
c     IW(IWORST+ICAND)=0
c     IW(JT1FT+JT1FTC)=ICAND
c     IW(JT1FT+JT1FTG)=0
           IWORST(ICAND)=0
           MY_JT1FTC(JT1FT)=ICAND
           MY_JT1FTG(JT1FT)=0

C GAB
c     C----->   Candidate -> 'T1TL'
c     KSCAL = KROW(KTCAL,ICAND)
c     KSGTL = KROW(KT1TL,NTKOK)
           KSCAL = ICAND
           KSGTL = NTKOK
           NT1TL=KSGTL

c     IW(KSGTL+JT1TIO) = NPTOK
c     IW(KSGTL+JT1TN1) = IW(KSCAL+JTCANF)
c     IW(KSGTL+JT1TNR) = IW(KSCAL+JTCANR)
           MY_JT1TIO(KSGTL) = NPTOK
           MY_JT1TN1(KSGTL) = MY_JTCANF(KSCAL)
           MY_JT1TNR(KSGTL) = MY_JTCANR(KSCAL)
C----->   Candidate -> 'T1CL'
           KSGCL = KKT1CL(NTKOK)

           
C     GAB:
c     NLONG = IW(KSGTL+JT1TN1)+IW(KSGTL+JT1TNR)
c     DO 80 I=1,NLONG
c     IW(KSGCL+I) = IW(KSTCC+I)
c     80     CONTINUE
           NLONG = MY_JT1TN1(KSGTL)+MY_JT1TNR(KSGTL)
           DO 80 I=1,NLONG
              MY_T1CLT(KSGCL+I)=MY_JTCCIT(MY_JTCAIO(ICAND)+I)
 80        CONTINUE
           NPTOK = NPTOK+NLONG

        ELSE
C-----> Else, determine which was the worst point
           CHI2  = 0.0
           

           CALL TMSHEC(DCHI2)
           
           DO 90 IPT =1,NPT
              IF(DCHI2(IPT).GT.CHI2) THEN
                 CHI2 = DCHI2(IPT)
C GAB:
C     IW(IWORST+ICAND)=IPT
                 IWORST(ICAND)=IPT
              END IF
 90        CONTINUE
        END IF
 100  CONTINUE

C
C*******************************************************************C
C                                                                   C
C       IMPROVE TRACKS BY REMOVING POINTS OF LARGE CHI**2           C
C                                                                   C
C*******************************************************************C
C
      DO 200 ICAND =1,NCAND
c GAB:
c        IF(IW(IWORST+ICAND).EQ.0) GO TO 200
        IF(IWORST(ICAND).EQ.0) GO TO 200

        NPT   = N1CORD(ICAND)
        IF(NPT.LE.MINPT)       GO TO 200
        KSTCC = KKTTCC(ICAND)
C GAB:
C        JT1FT = KROW(KT1FT,NTKOK+1)
C---->  Fill buffer with points
        JT1FT=NTKOK+1
        NT1FT=JT1FT

        DO 110 IPT=1,NPT
C GAB:
C          IBUFF(IPT) = IW(KSTCC+IPT)
          IBUFF(IPT)=MY_JTCCIT(MY_JTCAIO(ICAND)+IPT)
  110   CONTINUE
        LOST = 0
C

  111   CONTINUE
C
C---->  Remove worst point
C GAB:
C        IF (IW(IWORST+ICAND).LT.NPT) THEN
c          DO 112 IPT=IW(IWORST+ICAND)+1,NPT
        IF (IWORST(ICAND).LT.NPT) THEN
          DO 112 IPT=IWORST(ICAND)+1,NPT
            IBUFF(IPT-1)=IBUFF(IPT)
  112     CONTINUE
        END IF
        NPT = NPT-1
        LOST= LOST+1
C---- >  Now fit again
C GAB:
C     CALL TMSHEL('TPCO',NPT,IBUFF,2,RW(JT1FT+JT1FIR),
C     &                     RW(JT1FT+JT1FEM),RW(JT1FT+JT1FCD))
        CALL MY_TMSHEL(NPT,IBUFF,2,VV0,EE0,CHI2)
        
        MY_JT1FIR(JT1FT)=VV0(1)
        MY_JT1FTL(JT1FT)=VV0(2)
        DO I=1,15
           MY_JT1FEM(JT1FT,I)=EE0(I)
        ENDDO
        MY_JT1FCD(JT1FT)=CHI2
C     GAB:
c     IW(JT1FT+JT1FDF) = 2*NPT-5
c     PPP = PROB(RW(JT1FT+JT1FCD),(2*NPT-5))

        MY_JT1FDF(JT1FT) = 2*NPT-5

C SJA:Commented out the hard coded CHI2 condition inserted below, and reinstated the original value of PPP
        IF (MY_JT1FCD(JT1FT).GT.0.AND.(2*NPT-5).GT.0.AND.NPT.GT.0) THEN
           PPP = PROB(MY_JT1FCD(JT1FT),(2*NPT-5))
        ELSE
           WRITE(*,*) 
     +          "tfittk: CHI2 or NDF sent to CERNLIB PROB() are zero"
           PPP = 1e-25
        ENDIF

c     MVL        PPP = PROB(MY_JT1FCD(JT1FT),(2*NPT-5))
c        IF(MY_JT1FCD(JT1FT).GT.0.0.and.MY_JT1FCD(JT1FT).lT.300.0
c     >       .and.2*NPT-5.GT.0)THEN
c           PPP=max(PROB(MY_JT1FCD(JT1FT),(2*NPT-5)),1e-25)
c        ELSE
c           PPP=1e-25
c        ENDIF
c     MVL
C---- >  Take out the garbage
C     GAB:
c     CALL TFCHEK(RW(JT1FT+JT1FIR),.TRUE.,IERR)
        CALL TFCHEK(VV0,.TRUE.,IERR)
 
        IF (PPP.GT.PMIN .AND. IERR.EQ.0) THEN
          NTKOK = NTKOK+1

C GAB:
c          IW(IWORST+ICAND)=0
c          IW(JT1FT+JT1FTC)=ICAND
c          IW(JT1FT+JT1FTG)=0
          IWORST(ICAND)=0
          MY_JT1FTC(JT1FT)=ICAND
          MY_JT1FTG(JT1FT)=0

cC----->   Candidate -> 'T1TL'
c          KSCAL = KROW(KTCAL,ICAND)
c          KSGTL = KROW(KT1TL,NTKOK)

C GAB:
           KSCAL=ICAND
           KSGTL=NTKOK
           NT1TL=NTKOK
C GAB:
c          IW(KSGTL+JT1TIO) = NPTOK
c          IW(KSGTL+JT1TN1) = NPT
c          IW(KSGTL+JT1TNR) = IW(KSCAL+JTCANR)
          MY_JT1TIO(KSGTL) = NPTOK
          MY_JT1TN1(KSGTL) = NPT
          MY_JT1TNR(KSGTL) = MY_JTCANR(KSCAL)

C----->   Candidate -> 'T1CL'
          KSGCL = KKT1CL(NTKOK)
          DO 130 IPT=1,NPT
C GAB:
c            IW(KSGCL+IPT)=IBUFF(IPT)
            MY_T1CLT(KSGCL+IPT)=IBUFF(IPT)

  130     CONTINUE
C GAB
c          NGET  = IW(KSGTL+JT1TNR)
          NGET  = MY_JT1TNR(KSGTL)

          IF (NGET.GT.0) THEN
C GAB:
c            IGET  = KSTCC+IW(KSCAL+JTCANF)
c            IPUT  = KSGCL+NPT
            IGET  = MY_JTCANF(KSCAL)
            IPUT  = KSGCL+NPT

            DO 131 IPT=1,NGET
C GAB
c              IW(IPUT+IPT)=IW(IGET+IPT)
              MY_T1CLT(IPUT+IPT)=MY_JTCCIT(MY_JTCAIO(ICAND)+IPT)
  131       CONTINUE
          END IF
          NPTOK = NPTOK+NPT+NGET
        ELSE
C-----> Else, determine again, which was the worst point
          CHI2  = 0.0
          CALL TMSHEC(DCHI2)
          DO 140 IPT =1,NPT
          IF(DCHI2(IPT).GT.CHI2) THEN
            CHI2 = DCHI2(IPT)
C GAB:
C            IW(IWORST+ICAND)=IPT
            IWORST(ICAND)=IPT
          END IF
  140     CONTINUE
        END IF

c GAB:
c        IF(IW(IWORST+ICAND).NE.0
c     1                 .AND.LOST.LT.LOMAX
c     2                        .AND.NPT.GT.MINPT)       GO TO 111

        IF(IWORST(ICAND).NE.0
     1                 .AND.LOST.LT.LOMAX
     2                        .AND.NPT.GT.MINPT) THEN
         GO TO 111
        ENDIF
C          Jump back to 111 if another point is to be removed
C-----------------------------------------------------------------+
  200 CONTINUE

C
C*******************************************************************C
C                                                                   C
C              SEARCH FOR KINKS                                     C
C                                                                   C
C*******************************************************************C
C
C GAB:
       DO 300 ICAND =1,NCAND
C GAB
C        IF (IW(IWORST+ICAND).EQ.0)                 GO TO 300
        IF (IWORST(ICAND).EQ.0)                 GO TO 300
        NPT = N1CORD(ICAND)
        IF (NPT.LE.MKINK)                          GO TO 300
        KSTCC = KKTTCC(ICAND)
C GAB:
C        JT1FT = KROW(KT1FT,NTKOK+1)
C        CALL TKINK(DCF,KINK,IW(KSTCC+1),NPT)
         JT1FT=NTKOK+1
         NT1FT=JT1FT

         DO I=1,NPT
          LIST(I)= MY_JTCCIT(MY_JTCAIO(ICAND)+I)
         ENDDO

        CALL TKINK(DCF,KINK,LIST,NPT)
        IF (DCF.LT.DCM)                            GO TO 300
        IACC1 = 0
        IACC2 = 0
        NPT1 = KINK-1
        NPT2 = NPT-NPT1

C GAB: Check that the daughter tracks are viable:
        IF((NPT1.LT.MINPT).AND.(NPT2.LT.MINPT))GO TO 300

C-----> Kink is accepted, process 1st half of candidate
C
c GAB:
        IF(NPT1.GE.MINPT)THEN
         DO 210 IPT=1,NPT1
C GAB:
C          IBUFF(IPT) = IW(KSTCC+IPT)
           IBUFF(IPT) = LIST(IPT)
  210    CONTINUE
         LOST = 0

C
C GAB:
c         CALL TMSHEL('TPCO',NPT1,IBUFF,2,RW(JT1FT+JT1FIR),
c     &                    RW(JT1FT+JT1FEM),RW(JT1FT+JT1FCD))

         CALL MY_TMSHEL(NPT1,IBUFF,2,VV0,EE0,CHI2)
         MY_JT1FIR(JT1FT)=VV0(1)
         MY_JT1FTL(JT1FT)=VV0(2)
         DO I=1,15
          MY_JT1FEM(JT1FT,I)=EE0(I)
         ENDDO
         MY_JT1FCD(JT1FT)=CHI2

C-----> See if track is garbage
C GAB
C        CALL TFCHEK(RW(JT1FT+JT1FIR),.TRUE.,IERR)
         CALL TFCHEK(VV0,.TRUE.,IERR)

         IF (IERR.NE.0) THEN
C GAB:
C            IW(ICAND+IWORST) = 0
             IWORST(ICAND) = 0
             GOTO 240
         ENDIF
C GAB:
c        IW(JT1FT+JT1FDF) = 2*NPT1-5
c        PPP = PROB(RW(JT1FT+JT1FCD),(2*NPT1-5))

         MY_JT1FDF(JT1FT) = 2*NPT1-5
C SJA:Commented out the hard coded CHI2 condition inserted below, and reinstated the original value of PPP
         IF(MY_JT1FCD(JT1FT).GT.0.AND.(2*NPT1-5).GT.0.AND.NPT1.GT.0)THEN
            PPP = PROB(MY_JT1FCD(JT1FT),(2*NPT1-5))
         ELSE
            WRITE(*,*) 
     +           "tfittk: CHI2 or NDF sent to CERNLIB PROB() are zero"
            PPP = 1e-25
         ENDIF

c MVL         PPP = PROB(MY_JT1FCD(JT1FT),(2*NPT1-5))
C      IF(MY_JT1FCD(JT1FT).GT.0.0.and.MY_JT1FCD(JT1FT).lT.300.0
C     > .and.2*NPT1-5.GT.0)THEN
C        PPP=max(PROB(MY_JT1FCD(JT1FT),(2*NPT1-5)),1e-25)
C      ELSE
C        PPP=1e-25
C      ENDIF
c MVL

C GAB:
         IF(PPP.GE.PMIN)THEN
          IACC1 = 1
          NTKOK = NTKOK+1

C GAB:
c         IW(IWORST+ICAND)=0
c         IW(JT1FT+JT1FTC)=ICAND
c         IW(JT1FT+JT1FTG)=0
          IWORST(ICAND)=0
          MY_JT1FTC(JT1FT)=ICAND
          MY_JT1FTG(JT1FT)=0

C--->   Candidate -> 'T1TL'
c         KSCAL = KROW(KTCAL,ICAND)
c         KSGTL = KROW(KT1TL,NTKOK)
C GAB:
          KSCAL=ICAND
          KSGTL=NTKOK
          NT1TL=NTKOK
C GAB:
c        IW(KSGTL+JT1TIO) = NPTOK
c        IW(KSGTL+JT1TN1) = NPT1
c        IW(KSGTL+JT1TNR) = 0
          MY_JT1TIO(KSGTL) = NPTOK
          MY_JT1TN1(KSGTL) = NPT1
          MY_JT1TNR(KSGTL) = 0
C--->   Candidate -> 'T1CL'
          KSGCL = KKT1CL(NTKOK)
          DO 230 IPT=1,NPT1
C           IW(KSGCL+IPT)=IBUFF(IPT)
            MY_T1CLT(KSGCL+IPT)=IBUFF(IPT)
  230     CONTINUE
          NPTOK = NPTOK+NPT1
C GAB:
         ENDIF
        ENDIF
C
C-----> Now process 2nd half of candidate
C
  240   CONTINUE

C GAB:
        IF(NPT2.GT.MINPT)THEN

C GAB:
c         JT1FT = KROW(KT1FT,NTKOK+1)
         JT1FT=NTKOK+1
         NT1FT=JT1FT

         DO 250 IPT=1,NPT2
C GAB:
C          IBUFF(IPT) = IW(KSTCC+IPT+KINK-1)
           IBUFF(IPT) =MY_JTCCIT(MY_JTCAIO(ICAND)+IPT+KINK-1)
  250    CONTINUE
         LOST = 0
C GAB:
C         CALL TMSHEL('TPCO',NPT2,IBUFF,2,RW(JT1FT+JT1FIR),
C     &                     RW(JT1FT+JT1FEM),RW(JT1FT+JT1FCD))
         CALL MY_TMSHEL(NPT2,IBUFF,2,VV0,EE0,CHI2)
         MY_JT1FIR(JT1FT)=VV0(1)
         MY_JT1FTL(JT1FT)=VV0(2)
         DO I=1,15
          MY_JT1FEM(JT1FT,I)=EE0(I)
         ENDDO
         MY_JT1FCD(JT1FT)=CHI2

cC-----> Check for garbage
c         CALL TFCHEK(RW(JT1FT+JT1FIR),.TRUE.,IERR)
         CALL TFCHEK(VV0,.TRUE.,IERR)

         IF (IERR.NE.0) THEN
C GAB:
C             IW(IWORST+ICAND) = 0
             IWORST(ICAND) = 0
             GOTO 300
         ENDIF
C GAB
c        IW(JT1FT+JT1FDF) = 2*NPT2-5
          MY_JT1FDF(JT1FT) = 2*NPT2-5

C GAB:
C SJA:Commented out the hard coded CHI2 condition inserted below, and reinstated the original value of PPP
         IF(MY_JT1FCD(JT1FT).GT.0.AND.(2*NPT2-5).GT.0.AND.NPT2.GT.0)THEN
             PPP = PROB(MY_JT1FCD(JT1FT),(2*NPT2-5))
          ELSE
            WRITE(*,*) 
     +           "tfittk: CHI2 or NDF sent to CERNLIB PROB() are zero"
             PPP = 1e-25
          ENDIF

C MVL          PPP = PROB(MY_JT1FCD(JT1FT),(2*NPT2-5))
C      IF(MY_JT1FCD(JT1FT).GT.0.0.and.MY_JT1FCD(JT1FT).lT.300.0
C     > .and.2*NPT2-5.GT.0)THEN
C        PPP=max(PROB(MY_JT1FCD(JT1FT),(2*NPT2-5)),1e-25)
C      ELSE
C        PPP=1e-25
C      ENDIF
c MVL
         IF(PPP.GT.PMIN)THEN

          IACC2 = 1
          NTKOK = NTKOK+1

C GAB:
c         IW(IWORST+ICAND)=0
c         IW(JT1FT+JT1FTC)=ICAND
c         IW(JT1FT+JT1FTG)=0
          IWORST(ICAND)=0
          MY_JT1FTC(JT1FT)=ICAND
          MY_JT1FTG(JT1FT)=0

cC-----> Candidate -> 'T1TL'
c         KSCAL = KROW(KTCAL,ICAND)
c         KSGTL = KROW(KT1TL,NTKOK)
C GAB:
          KSCAL=ICAND
          KSGTL=NTKOK
          NT1TL=NTKOK
C GAB:
c         IW(KSGTL+JT1TIO) = NPTOK
c         IW(KSGTL+JT1TN1) = NPT2
c         IW(KSGTL+JT1TNR) = IW(KSCAL+JTCANR)
          MY_JT1TIO(KSGTL) = NPTOK
          MY_JT1TN1(KSGTL) = NPT2
          MY_JT1TNR(KSGTL) = MY_JTCANR(ICAND)


C--->   Candidate -> 'T1CL'
          KSGCL = KKT1CL(NTKOK)
          DO 260 IPT=1,NPT2
C GAB:
c           IW(KSGCL+IPT)=IBUFF(IPT)
            MY_T1CLT(KSGCL+IPT)=IBUFF(IPT)
  260     CONTINUE
C
C GAB:
c         NGET  = IW(KSGTL+JT1TNR)
          NGET  = MY_JT1TNR(KSGTL)

          IF (NGET.GT.0) THEN
C GAB:
C           IGET  = KSTCC+IW(KSCAL+JTCANF)
            IGET  = MY_JTCANF(KSCAL)
            IPUT  = KSGCL+NPT2
            DO 261 IPT=1,NGET
c              IW(IPUT+IPT)=IW(IGET+IPT)
              MY_T1CLT(IPUT+IPT)= MY_JTCCIT(MY_JTCAIO(ICAND)+IGET+IPT)
  261       CONTINUE
          END IF
          NPTOK = NPTOK+NPT2+NGET
C GAB:
         ENDIF
        ENDIF
C GAB:
c        IF ((IACC1+IACC2).GT.0) IW(IWORST+ICAND)=0
        IF ((IACC1+IACC2).GT.0) IWORST(ICAND)=0
  300  CONTINUE

C
C*******************************************************************C
C                                                                   C
C       IMPROVE TRACKS BY REMOVING POINTS OF LARGE CHI**2           C
C              (2-ND STRATEGY = FIND BEST CHI**2 )                  C
C                                                                   C
C*******************************************************************C
C
      DO 400 ICAND =1,NCAND
C GAB:
c        IF(IW(IWORST+ICAND).EQ.0) GO TO 400
        IF(IWORST(ICAND).EQ.0) GO TO 400
        NPT   = N1CORD(ICAND)
        IF(NPT.LE.MINPT)       GO TO 400
        KSTCC = KKTTCC(ICAND)
C---->  Fill buffer with points
        DO 310 IPT=1,NPT
C GAB:
c          IBUFF(IPT) = IW(KSTCC+IPT)
          IBUFF(IPT) = MY_JTCCIT(MY_JTCAIO(ICAND)+IPT)
  310   CONTINUE
        LOST = 0
C
C------ Jump back to 311 if another point is to be removed
C
  311   CONTINUE
C
c GAB: Check that removing the point doesn't put us below the limit
        IF(NPT.EQ.MINPT)GO TO 400

C---->  Remove worst point
        CALL TREMPT(PBEST,IBUFF,IBUFF,NPT)
        NPT = NPT-1
        LOST= LOST+1
        IF (PBEST.GE.PMIN) THEN
C          NTKOK = NTKOK+1

C GAB:
C          JT1FT = KROW(KT1FT,NTKOK+1)
          JT1FT = NTKOK+1
          NT1FT=JT1FT

C GAB:
c          CALL TMSHEL('TPCO',NPT,IBUFF,2,RW(JT1FT+JT1FIR),
c     &                     RW(JT1FT+JT1FEM),RW(JT1FT+JT1FCD))

         CALL MY_TMSHEL(NPT,IBUFF,2,VV0,EE0,CHI2)
         MY_JT1FIR(JT1FT)=VV0(1)
         MY_JT1FTL(JT1FT)=VV0(2)
         DO I=1,15
          MY_JT1FEM(JT1FT,I)=EE0(I)
         ENDDO
         MY_JT1FCD(JT1FT)=CHI2
C GAB:
cC----->   Check for garbage
c          CALL TFCHEK(RW(JT1FT+JT1FIR),.TRUE.,IERR)
          CALL TFCHEK(VV0,.TRUE.,IERR)

          IF (IERR.NE.0) THEN
C GAB:
c            IW(IWORST+ICAND) = 0
            IWORST(ICAND) = 0
            GOTO 400
          ENDIF
          NTKOK = NTKOK+1

C GAB:
c          IW(JT1FT+JT1FDF) = 2*NPT-5
c          IW(IWORST+ICAND)=0
c          IW(JT1FT+JT1FTC)=ICAND
c          IW(JT1FT+JT1FTG)=0
          MY_JT1FDF(JT1FT) = 2*NPT-5
          IWORST(ICAND)=0
          MY_JT1FTC(JT1FT)=ICAND
          MY_JT1FTG(JT1FT)=0

C----->   Candidate -> 'T1TL'
C GAB:
c          KSCAL = KROW(KTCAL,ICAND)
c          KSGTL = KROW(KT1TL,NTKOK)
           KSCAL=ICAND
           KSGTL=NTKOK
           NT1TL=NTKOK
C GAB:
c          IW(KSGTL+JT1TIO) = NPTOK
c          IW(KSGTL+JT1TN1) = NPT
c          IW(KSGTL+JT1TNR) = IW(KSCAL+JTCANR)
          MY_JT1TIO(KSGTL) = NPTOK
          MY_JT1TN1(KSGTL) = NPT
          MY_JT1TNR(KSGTL) = MY_JTCANR(KSCAL)

C----->   Candidate -> 'T1CL'
          KSGCL = KKT1CL(NTKOK)
          DO 330 IPT=1,NPT
C GAB:
c            IW(KSGCL+IPT)=IBUFF(IPT)
           MY_T1CLT(KSGCL+IPT)=IBUFF(IPT)

  330     CONTINUE

C GAB:
C          NGET  = IW(KSGTL+JT1TNR)
          NGET  = MY_JT1TNR(KSGTL)
          IF (NGET.GT.0) THEN
C GAB:
C            IGET  = KSTCC+IW(KSCAL+JTCANF)
            IGET  = KSTCC+MY_JTCANF(KSCAL)
            IPUT  = KSGCL+NPT
            DO 331 IPT=1,NGET
C GAB
c              IW(IPUT+IPT)=IW(IGET+IPT)
           MY_T1CLT(IPUT+IPT)=MY_JTCCIT(MY_JTCAIO(ICAND)+IGET+IPT)

  331       CONTINUE
          END IF
          NPTOK = NPTOK+NPT+NGET
        END IF
C GAB:
c        IF (IW(IWORST+ICAND).NE.0
c     1                .AND.LOST.LT.LOMAX
c     2                      .AND.NPT.GT.MINPT)       GO TO 311
        IF (IWORST(ICAND).NE.0
     1                .AND.LOST.LT.LOMAX
     2                      .AND.NPT.GT.MINPT)       GO TO 311
cC----------- Jump back to 111 if another point is to be removed
  400 CONTINUE
C
C*******************************************************************C
C                                                                   C
C               KEEP REMAINING TRACKS AS THEY ARE                   C
C GAB:  PROVIDED THEY PASS THE N_HIT AND MIN_PROB REQUIREMENTS      C
C                                                                   C
C*******************************************************************C
C

      DO 500 ICAND =1,NCAND
c GAB:
c        IF (IW(IWORST+ICAND).EQ.0) GO TO 500
c        JT1FT = KROW(KT1FT,NTKOK+1)
        IF (IWORST(ICAND).EQ.0) GO TO 500
        JT1FT = NTKOK+1

c GAB:
        IF(NPT.GE.MINPT)THEN

         NPT   = N1CORD(ICAND)
C GAB:
c         KSTCC = KKTTCC(ICAND)
c         CALL TMSHEL('TPCO',NPT,IW(KSTCC+1),2,RW(JT1FT+JT1FIR),
c     &                         RW(JT1FT+JT1FEM),RW(JT1FT+JT1FCD))
         DO I=1,NPT
          LIST(I)= MY_JTCCIT(MY_JTCAIO(ICAND)+I)
         ENDDO
         CALL MY_TMSHEL(NPT,LIST,2,VV0,EE0,CHI2)
         MY_JT1FIR(JT1FT)=VV0(1)
         MY_JT1FTL(JT1FT)=VV0(2)
         DO I=1,15
          MY_JT1FEM(JT1FT,I)=EE0(I)
         ENDDO
         MY_JT1FCD(JT1FT)=CHI2

C--->   Garbage check
C GAB:
C        CALL TFCHEK(RW(JT1FT+JT1FIR),.TRUE.,IERR)
         CALL TFCHEK(VV0,.TRUE.,IERR)

         IF (IERR.NE.0) THEN
C GAB:
C            IW(IWORST+ICAND) = 0
             IWORST(ICAND) = 0
             GOTO 500
         ENDIF
C GAB:

C SJA:Commented out the hard coded CHI2 condition inserted below, and reinstated the original value of PPP
        IF ( CHI2.GT.0.AND.(2*NPT-5).GT.0.AND.NPT.GT.0) THEN
           PPP = PROB(CHI2,(2*NPT-5))
        ELSE
           WRITE(*,*) 
     +          "tfittk: CHI2 or NDF sent to CERNLIB PROB() are zero"
           PPP = 1e-25
        ENDIF

c MVL         PPP=PROB(CHI2,2*NPT-5)
C      IF(CHI2.GT.0.0.and.CHI2.lT.300.0.and.2*NPT-5.GT.0)THEN
C        PPP=max(PROB(CHI2,(2*NPT-5)),1.e-25)
C      ELSE
C        PPP=1.e-25
C      ENDIF
c MVL
         IF(PPP.GE.PMIN)THEN
C GAB:
c         IW(JT1FT+JT1FDF) = 2*NPT-5
          MY_JT1FDF(JT1FT) = 2*NPT-5

C
C--->      IWORST = 0       -> Track candidate has been accepted
C                           (track not yet accepted)
           NTKOK = NTKOK+1

C GAB:
c         IW(IWORST+ICAND)=0
c         IW(JT1FT+JT1FTC)=ICAND
c         IW(JT1FT+JT1FTG)=0
          IWORST(ICAND)=0
          MY_JT1FTC(JT1FT)=ICAND
          MY_JT1FTG(JT1FT)=0
C---->  Candidate -> 'T1TL'
C GAB:
c         KSCAL = KROW(KTCAL,ICAND)
c         KSGTL = KROW(KT1TL,NTKOK)
          KSCAL = ICAND
          KSGTL = NTKOK
          NT1TL=NTKOK
C GAB:
c          IW(KSGTL+JT1TIO) = NPTOK
c          IW(KSGTL+JT1TN1) = IW(KSCAL+JTCANF)
c          IW(KSGTL+JT1TNR) = IW(KSCAL+JTCANR)
          MY_JT1TIO(KSGTL) = NPTOK
          MY_JT1TN1(KSGTL) = MY_JTCANF(KSCAL)
          MY_JT1TNR(KSGTL) = MY_JTCANR(KSCAL)

cC-----> Candidate -> 'T1CL'
          KSGCL = KKT1CL(NTKOK)
C GAB:
C          NLONG = IW(KSGTL+JT1TN1)+IW(KSGTL+JT1TNR)
          NLONG = MY_JT1TN1(KSGTL)+MY_JT1TNR(KSGTL)
          DO 480 I=1,NLONG
C GAB:
c            IW(KSGCL+I) = IW(KSTCC+I)
            MY_T1CLT(KSGCL+I)= MY_JTCCIT(MY_JTCAIO(ICAND)+I)
  480     CONTINUE
          NPTOK = NPTOK+NLONG
         ENDIF
        ENDIF
  500 CONTINUE

cC
cC---> Fill remaining words in banks and set the lengths
cC
c      IW(KT1FT+LMHROW)=NTKOK
c      IW(KT1TL+LMHROW)=NTKOK
c      IW(KT1CL+LMHROW)=NPTOK

C GAB:
       NT1FT=NTKOK
       NT1CL=NPTOK


c      CALL AUBPRS('T1FTT1TLT1CL')
c      CALL BLIST(IW,'T+','T1FTT1TLT1CL')
C

C*******************************************************************C
C                                                                   C
C              RELABEL THE 'TPCO' BANK                              C
C                                                                   C
C*******************************************************************C
C
c GAB:
c      NCOR = LROWS(KTPCO)
       NCOR=NTPCO
      DO 800 I=1,NCOR
C        IW(KROW(KTPCO,I)+JTPCTN)=0
        TPCO_TN(I)=0
  800 CONTINUE
C GAB:
C      NTK = LROWS(KT1FT)
      NTK=NT1FT

      DO 820 ITK = 1,NTK
C GAB:
c        JT1TL = KROW(KT1TL,ITK)
         JT1TL=ITK
C GAB:
c        I1 = IW(JT1TL+JT1TIO)
c        I2 = IW(JT1TL+JT1TN1)
c        I3 = IW(JT1TL+JT1TNR)
        I1 = MY_JT1TIO(JT1TL)
        I2 = MY_JT1TN1(JT1TL)
        I3 = MY_JT1TNR(JT1TL)

        DO 810 IPT =1,I2
C GAB:
C          ICO = IW((KT1CL+LMHLEN+I1)+IPT)
C          IW(KROW(KTPCO,ICO)+JTPCTN) = ITK
          ICO = MY_T1CLT(I1+IPT)
          TPCO_TN(ICO) = ITK
  810   CONTINUE
        DO 813 IPT=1,I3
C GAB:
c          ICO = IW((KT1CL+LMHLEN+I1+I2)+IPT)
c          IW(KROW(KTPCO,ICO)+JTPCTN) = -ITK
          ICO = MY_T1CLT(I1+I2+IPT)
c(kh) bug workaround
c(kh)     TPCO_TN(ICO) = -ITK
          IF (ICO.GT.0) THEN
             TPCO_TN(ICO) = -ITK
          ELSE
             WRITE(6,*) 'TFITTK: internal error. ICO=',ICO
             CALL TKVETO(1)
          END IF
  813   CONTINUE
  820 CONTINUE
C
C     Call the subroutine that converts the parameter definition
C     from TASSO convention to standard convention.
C     This is done by actually refitting the track yet again and
C     overwriting the T1FT bank which was filled above.
C
      CALL UCTGFT
C
      GO TO 1000
c  995 CONTINUE
c      IER=1
c      CALL RERROR('TFITTK',IER,'No room in BOS for workbank')
c      GO TO 1000
c  998 CONTINUE
c      IER = 2
c      CALL RERROR('TFITTK',IER,'Input banks not found')
c      GO TO 1000
c  999 CONTINUE
c      IER = 3
c      CALL RERROR('TFITTK',IER,'Not enough space in IW')
c      CALL BDROP(IW,'T1FTT1TLT1CL')
cC
 1000 CONTINUE
c      CALL WDROP(IW,IWORST)

      RETURN

      END

c-------------------------------------------------------------------
