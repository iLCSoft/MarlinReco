*  
* $Id: tflnks.F,v 1.6 2008-06-26 09:35:56 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.5  2008/01/23 17:22:22  aplin
*
*      The TPC Tracking code has been altered to cope with more than 224 pad rows
*      the max number of pad rows can be changed in padrow.inc by changing
*
*      The maximum number of TPC Padrows that the tracking code
*      can accept must be a multiple of 32
*
*      Note this is not the number of TPC Padrows in the dectector
*      but the maximum allowed because of the predifined array sizes
*      and a set of 32 bit registers
*
*      Currently it is set to 512
*
*      By increasing the value of N32BITREG in padrow.inc you can increase the max
*      number of pad rows
*
* Revision 1.4  2006/02/24 16:36:55  aplin
* increased NMXHTRW again by a factor of 2 in thitsrt.F
* increased the minimum number of hits to form a chain in trkfnd.F from MNLNK=4 to MNLNK=7
* set the minimum radius/pad ring for inward search in first track finding to MINMRG=120 in unmerge.F
* reduced the maximum number of links NMAXLNK to 25000
*
* Revision 1.3  2006/02/22 12:19:21  owendt
* merged new non-cvs versions of reco code (O.Wendt/F.Gaede)
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
      SUBROUTINE TFLNKS(IER)
      IMPLICIT NONE
C****************************************************************************
c* Version number: 0.1      Date of last modification   6.11.99             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
c* Addional code to join links of loopers                                   *
C* Main modifiers from ALEPH to BRAHMS code:                                *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
C-----*****************---------------------------------------------
C! FIND TPC LINKS BETWEEN CHAINS
C!
C!    AUTHOR        :- MARTIN POPPE 87/5/25
C!    Modified      :  P. Colrain 95/10/26
C!                     D. Casper  98/2/04  Correct typo in "long-range" code
C!                                         Select best link, not first...
C!
C!    OUTPUT:  IER = ERROR CODE:  =0 IF OK, =1 IF INCOMPLETE
C!
C!    INPUT  BANKS  : 'TCHA' , 'TCTC'
C!    OUTPUT BANK   : 'TLNK'
C!
C!
C!    DESCRIPTION
C!    ===========
C!
C!    Creates a table of links between chains defining
C!    which chains belong to the same particle. This table is
C!    stored in the bank 'TLNK' which has as many rows as there
C!    are chains. Four types of links  are possible: "up", "down"
C!    "top" , and "bottom". These links are stored in columns 1 to 4.
C!    For example, the number 16 in row 3 of column 10 means that
C!    there is a "top" link between chains 10 an 16, i.e. this is
C!    a spiraling particle.
C-------------------------------------------------------------------C

C
#include "include/tkffread.inc"
#include "include/padrow.inc"
#include "include/fndparms.inc"
#include "include/alcons.inc"
#include "include/tpco.inc"
#include "include/tchajj.inc"
#include "include/tpgeom.inc"
c#include "marlin_tpcgeom.h"
#include "include/tpgeop.inc"
#include "include/my_ttrc.inc"
#include "include/my_tcha.inc"
#include "include/my_tlnk.inc"
#include "include/my_tctc.inc"
c#include "include/gcflag.inc"
#include "include/tpcdims.inc"
#include "include/patrec.inc"
C GAB: for fitting chains:
#include "include/tchain.inc"
#include "include/my_jndb.inc"

      INTEGER IER
      INTEGER*4 ITPTES,IOUT1,IOUT2,IOUT3,IOUT4

      INTEGER MPT
c(KH) PARAMETER(MPT=200)
      PARAMETER(MPT=LTPDRO)
      DOUBLE PRECISION XF(MPT),YF(MPT),WF(MPT)
      REAL RF(MPT),PF(MPT),ZF(MPT),WZF(MPT)

      REAL CHISZ
      INTEGER ITSTL(N32BITREG),ITSTU(N32BITREG)
      DOUBLE PRECISION FLOAT2,FLOAT3,FLOAT4,FLOAT5,FLOAT6,FLOAT7
C-----------------------------------------------------------------------
      INTEGER MCAND
      REAL EPS
C changed MCAND=900 to MCAND=1800 because of released energy cut in Mokka 5.04
C      PARAMETER(MCAND=900,EPS=0.0000001)
      PARAMETER(MCAND=1800,EPS=0.0000001)
C(KH) INTEGER LILLY(LTPDRO)

C changed LILLY(500) to LILLY(2000) because of released energy cut in Mokka 5.04
C(SJA) set back to LTPDRO to ensure control over all max row based arrays
      INTEGER LILLY(LTPDRO)
      REAL VAL(5),ERR(15)
c GAB:
c      INTEGER IBUF1(MCAND),IBUF2(MCAND),IPATL(MCAND)
      INTEGER IBUF1(MCAND),IBUF2(MCAND)
c GAB:
c      INTEGER I1PAT(MCAND),I2PAT(MCAND)
c(kh) INTEGER*4 I1PAT(MCAND,4),I2PAT(MCAND,4)
      INTEGER*4 I1PAT(MCAND,N32BITREG),I2PAT(MCAND,N32BITREG)
      LOGICAL LB2(MCAND)
      LOGICAL FIRST,LFUND
      DATA FIRST/.TRUE./

      INTEGER I,J,JJ,K1,K2,IF15,IL06,NCNDU,NCNDL,LSTPT,LB1PT,LSTCO,LB1CO
      INTEGER ICNDU,JCH1,JCH2,NHITS,ICNDL,ITEST1,ITEST2,ITEST3,ITEST4
      INTEGER ITEST5,ITEST6,ITEST7
      INTEGER ITEST(N32BITREG)
      INTEGER ITOP,IBOT,LOOPL,NQL,NQU,IPT,NLILL,IDUM1,IDUM2,IDUM3,NCNDR
      INTEGER ICNDR,LOOPR,ISEL
      REAL CHMAX,CHTMX,CHPMX,OCUT,DROW,CHFRE,R1,R2,Z1,Z2,ZEXT,RMATK
      REAL DR1,DR2,FUDGE,DPHI,CHI2,CHX,RMINI,RLRGE,FUDG1,RMTRK
      DOUBLE PRECISION FACTOR1,FACTOR2,FTMP(N32BITREG)

c      INTEGER LINK_BEST
c      REAL CHI2_BEST

      REAL OME,ALP,PH0,DD0,ZZ0,CHF,RADTK,RTOP,EOME,EALP,EPH0,EDD0
      REAL EOMES,EALPS,EPH0S,EDD0S,EZZ0S,EZZ0

      SAVE CHMAX,CHTMX,CHPMX,IOUT1,IOUT2,IOUT3,IOUT4,IF15,IL06,
     & FLOAT2,FLOAT3,FLOAT4,FLOAT5,FLOAT6,FLOAT7,OCUT

C GAB: Top/bottom chain linking variables
       INTEGER IH1BEG,IH1END,IH2BEG,IH2END,K,IHIT,ICH,ILAST,IFRST
       INTEGER KMIN,ISWAP
       REAL D1,D2,RMIN
C
C-------------------- Statement Function Definitions ----------------
C
#include "include/bkgmode.inc"
c#include "include/gcbank.inc"
#include "include/tpcfunc2.inc"

#include "include/tpcfunc.inc"


C-----> track parameters from 'TCHA'  follow now:
      OME(JJ)  = MY_JTCHIR(JJ)
C-----> inv. rad of track from offset in 'TCHA'
      ALP(JJ)  = MY_JTCHTL(JJ)
C-----> TAN(LAM) of track from offset in 'TCHA'
      PH0(JJ)  = MY_JTCHP0(JJ)
C-----> PHI(0)   of track from offset in 'TCHA'
      DD0(JJ)  = MY_JTCHD0(JJ)
C-----> IMPACT   of track from offset in 'TCHA'
      ZZ0(JJ)  = MY_JTCHZ0(JJ)
C-----> Z(0)     of track from offset in 'TCHA'
      CHF(JJ)  = MY_JTCHC1(JJ)
C-----> chi**2/freedom of track from offset in 'TCHA'
C-----> derived track quantities  follow
      RADTK(JJ)= 1.0/ABS(MY_JTCHIR(JJ)+0.01*EPS)
      RTOP(JJ) = 2.0*RADTK(JJ)-MY_JTCHD0(JJ)
C-----> inverse errors squared on track quantities follow now:
      EOME(JJ)= MY_JTCHER(JJ,1)
C                           ... FOR 1/R
      EALP(JJ)= MY_JTCHER(JJ,3)
C                           ... FOR TAN(LAM)
      EPH0(JJ)= MY_JTCHER(JJ,6)
C                           ... FOR PHI(0)
      EDD0(JJ)= MY_JTCHER(JJ,10)
C                           ... FOR D(0)
      EZZ0(JJ)= MY_JTCHER(JJ,15)
C                           ... FOR Z(0)

C-----> error weights scaled with mean chi**2/freedom
      EOMES(JJ)= MY_JTCHER(JJ,1)/(MY_JTCHC2(JJ)+1.5)
C                           ... FOR 1/R
      EALPS(JJ)= MY_JTCHER(JJ,3)/(MY_JTCHC2(JJ)+1.5)
C                           ... FOR TAN(LAM)
      EPH0S(JJ)= MY_JTCHER(JJ,6)/(MY_JTCHC2(JJ)+1.5)
C                           ... FOR PHI(0)
      EDD0S(JJ)= MY_JTCHER(JJ,10)/(MY_JTCHC2(JJ)+1.5)
C                           ... FOR D(0)
      EZZ0S(JJ)= MY_JTCHER(JJ,15)/(MY_JTCHC2(JJ)+1.5)
C                           ... FOR Z(0)

cC
cC------------------------------------------------------------------
cC
      IER=0
      IF (FIRST) THEN
        CHMAX = 100.0
        CHTMX = 500.0
        CHPMX = 1000.0

c(kh)   IOUT1 = 2**(LTPDRO-96)
c(kh)   IOUT2 = 2**(LTPDRO-96-1)
c(kh)   IOUT3 = 2**(LTPDRO-96-2)
c(kh)   IOUT4 = 2**(LTPDRO-96-3)
        IOUT1 = 2**(NRTPC-(NRTPC/32)*32)
        IOUT2 = 2**(NRTPC-(NRTPC/32)*32-1)
        IOUT3 = 2**(NRTPC-(NRTPC/32)*32-2)
        IOUT4 = 2**(NRTPC-(NRTPC/32)*32-3)
        IF15  = 2**15 - 1
c(kh)   IL06  = (2**7-1) * (2**(LTPDRO-96-6))
c(SJA)        IL06  = (2**7-1) * (2**(NRTPC-(NRTPC/32)*32-6))
        IL06  = (2**N32BITREG-1) * 
     &       (2**(NRTPC-(NRTPC/32)*32-(N32BITREG-1)))



        FLOAT2=2.D0**32
        FLOAT3=2.D0**64
        FLOAT4=2.D0**96
        FLOAT5=2.D0**128
        FLOAT6=2.D0**160
        FLOAT7=2.D0**192

        OCUT= 1.0/(5.0*TPCACRO)

        FIRST = .FALSE.
      END IF
      DROW= TPDRST(1)

C
C-->  Get steering constants from the TTRC bank
C
        CHFRE = MY_JTTRCC
        CHMAX = MY_JTTRSU**2  * 10.0
        CHTMX = MY_JTTRST**2  * 10.0
        CHPMX = MY_JTTRSI**2  * 2.0
C
C-------------------------------------------------------------------C
C              SELECT CANDIDATES FOR UP/DOWN LINKS                  C
C              ***********************************                  C
C-------------------------------------------------------------------C
C---> Find candidates for upper chains in up-down links
C     criterium: at most 2 hits in the inner 5 pad rows
C
      NCNDU = 0
c GAB:
       DO I=1,MY_NCHAIN
        MY_JTLNUP(I)=0
        MY_JTLNDN(I)=0
        MY_JTLNTP(I)=0
        MY_JTLNBM(I)=0
       ENDDO

C GAB:
      DO 10 I=1,MY_NCHAIN
        IF (NCNDU.EQ.MCAND) THEN
         IF(TKSTDBG.GT.0) WRITE(*,*)' TFLNKS: MCAND too small'
         CALL TKVETO(1)
         GO TO 11
        ENDIF

C GAB:
         ITPTES = MY_JTCHPI(I,1)
        NHITS = IAND(ITPTES,32)/32 +IAND(ITPTES,16)/16 +
     $          IAND(ITPTES,8)/8   +IAND(ITPTES,4)/4   +IAND(ITPTES,2)/2
        IF (NHITS.GT.2) GO TO 10
        NCNDU = NCNDU+1

C GAB:
        IF(NCNDU.GT.MCAND)THEN
         IF(TKSTDBG.GT.0) WRITE(*,*)'TFLNKS:MCAND too small'
         CALL TKVETO(1)
         IER=1
         RETURN
        ENDIF

        IBUF1(NCNDU) = I
C GAB:
        DO J=1,N32BITREG
         I1PAT(NCNDU,J) = MY_JTCHPI(I,J)
        ENDDO
   10 CONTINUE
   11 CONTINUE
C
C---> IBUF1 now contains the chains which might be an upper
C     part of a set of linked chains
C
      IF (NCNDU.EQ.0) GO TO 100
C
C---> Now find candidates for lower chains in up-down links
C     Criteria: 1) at most 2 hits in the outer 4 pad rows
C               2) no extrapolation outside the chamber in z
C               3) not being closer than 3 rows to the max radius
C
      NCNDL  = 0
c GAB:
c      DO 20 I=1,LROWS(KTCHA)
      DO 20 I=1,MY_NCHAIN
        IF (NCNDL.EQ.MCAND)THEN
         IF(TKSTDBG.GT.0) WRITE(*,*)'TFLNKS:MCAND too small'
         CALL TKVETO(1)
         GO TO 21
        ENDIF

C GAB:
C        ITPTES = ITABL(KTCHA,I,JTCHPI)
        ITPTES = MY_JTCHPI(I,NRTPC/32+1)

        NHITS = IAND(ITPTES,IOUT1)/IOUT1 +IAND(ITPTES,IOUT2)/IOUT2
     1         +IAND(ITPTES,IOUT3)/IOUT3 +IAND(ITPTES,IOUT4)/IOUT4
        IF (NHITS.GT.2)  THEN
         GO TO 20
        ENDIF
C GAB:
C        LSTPT = ITABL(KTCHA,I,JTCHOF) + ITABL(KTCHA,I,JTCHNC)
        LSTPT=MY_JTCHOF(I)+MY_JTCHNC(I)

        LB1PT = LSTPT - 1


        LSTCO = MY_TCTC(LSTPT)
        LB1CO = MY_TCTC(LB1PT)


        R1= TPCO_RV(LB1CO)
        R2= TPCO_RV(LSTCO)
        Z1= TPCO_ZV(LB1CO)
        Z2= TPCO_ZV(LSTCO)

C GAB:
        IF(R1.EQ.R2)  GO TO 20
        ZEXT   = ABS(Z2 + 2.0*DROW*(Z2-Z1)/(R2-R1))
        IF (ZEXT.GT.ZDRIFT)                           GO TO 20
        RMATK  = RTOP(I)

        DR1 = 0.25*(RMATK-R1)
        DR2 = 0.25*(RMATK-R2)
        IF (DR1.LT.DROW.AND.DR2.LT.DROW)          GO TO 20

C GAB:
        NCNDL = NCNDL+1
        IF(NCNDL.GT.MCAND)THEN
         IF(TKSTDBG.GT.0) WRITE(*,*)'TFLNKS:MCAND too small'
         CALL TKVETO(1)
         IER=1
         RETURN
        ENDIF

        IBUF2(NCNDL) = I
C GAB:
C        I2PAT(NCNDL) = ITPTES
        DO J=1,N32BITREG
         I2PAT(NCNDL,J) = MY_JTCHPI(I,J)
        ENDDO
        LB2(NCNDL)   = .FALSE.
   20 CONTINUE
   21 CONTINUE


      IF (NCNDL.LE.0)                                GO TO 100
C
C-------------------------------------------------------------------C
C              DEFINE UP/DOWN LINKS BETWEEN CANDIDATES              C
C              ***************************************              C
C-------------------------------------------------------------------C
C
      DO 40 ICNDU = 1,NCNDU

        JCH1=IBUF1(ICNDU)
C GAB:
C        ITSTU = I1PAT(ICNDU)
        DO J=1,N32BITREG
         ITSTU(J) = I1PAT(ICNDU,J)
        ENDDO
        FUDGE  = 1.E30
        LFUND = .FALSE.

        DO 30 ICNDL = 1,NCNDL
          IF (LB2(ICNDL))                      GO TO 30

          JCH2=IBUF2(ICNDL)
C GAB:
C          ITSTL = I2PAT(ICNDL)
          DO J=1,N32BITREG
           ITSTL(J) = I2PAT(ICNDL,J)
          ENDDO

C
C-->      Check consistency of hit pattern
C
C GAB:
C           IF (ITSTL.GT.ITSTU)                 GO TO 30


          IF(ITSTL(N32BITREG).GT.ITSTU(N32BITREG)) GOTO 30
                
          DO J=N32BITREG,2,-1
             IF( (ITSTU(J).EQ.0).AND.
     &            (ITSTL(J-1).GT.ITSTU(J-1)))THEN
                GO TO 30 
             ENDIF
          ENDDO

          
C          IF(ITSTL(7).GT.ITSTU(7))THEN
C           GO TO 30
C          ELSE IF( (ITSTU(7).EQ.0).AND.
C     &             (ITSTL(6).GT.ITSTU(6)))THEN
C           GO TO 30
C          ELSE IF( (ITSTU(6).EQ.0).AND.
C     &             (ITSTL(5).GT.ITSTU(5)))THEN
C           GO TO 30
C          ELSE IF( (ITSTU(5).EQ.0).AND.
C     &             (ITSTL(4).GT.ITSTU(4)))THEN
C           GO TO 30
C          ELSE IF( (ITSTU(4).EQ.0).AND.
C     &             (ITSTL(3).GT.ITSTU(3)))THEN
C           GO TO 30
C          ELSE IF( (ITSTU(3).EQ.0).AND.
C     &             (ITSTL(2).GT.ITSTU(2)))THEN
C           GO TO 30
C          ELSE IF( (ITSTU(2).EQ.0).AND.
C     &             (ITSTL(1).GT.ITSTU(1)))THEN
C           GO TO 30
C          ENDIF

C
          ITPTES = 0
C GAB:
C          ITPTES = IAND(ITSTL,ITSTU)

          DO J=1,N32BITREG
          ITEST(J) = IAND(ITSTL(J),ITSTU(J))
          ENDDO

C          ITEST1 = IAND(ITSTL(1),ITSTU(1))
C          ITEST2 = IAND(ITSTL(2),ITSTU(2))
C          ITEST3 = IAND(ITSTL(3),ITSTU(3))
C          ITEST4 = IAND(ITSTL(4),ITSTU(4))
C          ITEST5 = IAND(ITSTL(5),ITSTU(5))
C          ITEST6 = IAND(ITSTL(6),ITSTU(6))
C          ITEST7 = IAND(ITSTL(7),ITSTU(7))

C GAB:
C          IF (ITPTES.NE.0)                       GO TO 30
C          IF ((ITEST1.NE.0).OR.(ITEST2.NE.0).OR.(ITEST3.NE.0).OR.
C     &         (ITEST4.NE.0).OR.(ITEST5.NE.0).OR.(ITEST6.NE.0)
C     &         .OR.(ITEST7.NE.0))GO TO 30

          DO J=1,N32BITREG
             IF (ITEST(J).NE.0) GOTO 30
          ENDDO

C
C-->      Check consistency of track parameters
C
c GAB:
c          DPHI=RTABL(KTCHA,JCH1,JTCHP0)-RTABL(KTCHA,JCH2,JTCHP0)
          DPHI=MY_JTCHP0(JCH1)-MY_JTCHP0(JCH2)

          IF (ABS(DPHI) .GT. PI) DPHI=TWOPI-ABS(DPHI)

          CHI2 = (MY_JTCHIR(JCH1)
     1                        -MY_JTCHIR(JCH2))**2
     2          *EOMES(JCH1)*EOMES(JCH2)/(EOMES(JCH1)+EOMES(JCH2))
     3          + (MY_JTCHTL(JCH1)
     4                        -MY_JTCHTL(JCH2))**2
     5          *EALPS(JCH1)*EALPS(JCH2)/(EALPS(JCH1)+EALPS(JCH2))
     6          + DPHI**2
     7          *EPH0S(JCH1)*EPH0S(JCH2)/(EPH0S(JCH1)+EPH0S(JCH2))
     8          + (MY_JTCHD0(JCH1)
     9                        -MY_JTCHD0(JCH2))**2
     1          *EDD0S(JCH1)*EDD0S(JCH2)/(EDD0S(JCH1)+EDD0S(JCH2))
     2          + (MY_JTCHZ0(JCH1)
     3                        -MY_JTCHZ0(JCH2))**2
     4          *EZZ0S(JCH1)*EZZ0S(JCH2)/(EZZ0S(JCH1)+EZZ0S(JCH2))

C
          CHI2 = 0.2*CHI2

          IF (CHI2.LT.CHMAX) THEN

           DO I=1,N32BITREG
            FTMP(I)=DBLE(ITSTU(I))
C            IF(FTMP(I).LT.0) FTMP(I)=(FTMP(I)+FLOAT2)
            IF(FTMP(I).LT.0) FTMP(I)=(FTMP(I)+(2.D0**32))
           ENDDO

C this should become a loop

           FACTOR1 = 0
           DO I=1,N32BITREG
              FACTOR1 = FACTOR1 + FTMP(I)*(2.D0**(32*(I-1))) 
           ENDDO

C           FACTOR1=FTMP(1)  +  FTMP(2)*FLOAT2 +
C     &             FTMP(3)*FLOAT3 + FTMP(4)*FLOAT4 +
C     &             FTMP(5)*FLOAT5+FTMP(6)*FLOAT6+FTMP(7)*FLOAT7


           DO I=1,N32BITREG
            FTMP(I)=DBLE(I2PAT(ICNDL,I))
C            IF(FTMP(I).LT.0)FTMP(I)=(FTMP(I)+FLOAT2)
            IF(FTMP(I).LT.0)FTMP(I)=(FTMP(I)+(2.D0**32))
           ENDDO

C this should also become a loop

           FACTOR2 = 0
           DO I=1,N32BITREG
              FACTOR2 = FACTOR2 + FTMP(I)*(2.D0**(32*(I-1))) 
           ENDDO

C           FACTOR2=FTMP(1)  +  FTMP(2)*FLOAT2 +
C     &             FTMP(3)*FLOAT3 + FTMP(4)*FLOAT4 +
C     &             FTMP(5)*FLOAT5+FTMP(6)*FLOAT6+FTMP(7)*FLOAT7
     
           FUDG1 = CHI2*DLOG(0.5D0*FACTOR1/FACTOR2)

C
C---->      Minimise SQRT( DELTA(R/MAX) ) * CHI**2
C
            IF (FUDG1.LT.FUDGE) THEN
              LFUND = .TRUE.
              FUDGE  = FUDG1
              ITOP   = IBUF1(ICNDU)
              IBOT   = IBUF2(ICNDL)
              LOOPL  = ICNDL
            END IF

          END IF
   30   CONTINUE

C
C-----> Test via fit over all coordinates
C
        IF (LFUND) THEN
C GAB:
C          NQL = ITABL(KTCHA,IBOT,JTCHNC)
C          NQU = ITABL(KTCHA,ITOP,JTCHNC)

          NQL = MY_JTCHNC(IBOT)
          NQU = MY_JTCHNC(ITOP)

          DO 36 IPT = 1,NQL
C GAB:
C            LILLY(IPT) = IW(LIST(IBOT)+IPT)

C           WRITE(0,*) "LILLY(NQL+IPT)LILLY = ", (NQL+IPT)
C           WRITE(0,*) "IBOT = " , IBOT
C           WRITE(0,*) "MY_JTCHOF(IBOT)+IPT = " , MY_JTCHOF(IBOT)+IPT
            LILLY(IPT) = MY_TCTC(MY_JTCHOF(IBOT)+IPT)
   36     CONTINUE

          NLILL  = NQU+NQL
          DO 37 IPT = 1,NQU
C GAB:
C            LILLY(NQL+IPT) = IW(LIST(ITOP)+IPT)
            LILLY(NQL+IPT) = MY_TCTC(MY_JTCHOF(ITOP)+IPT)
   37     CONTINUE

c GAB:
c          CALL TFTHEL('TPCO',NLILL,LILLY,1,VAL,ERR,CHI2)

c(KH)      DO I=1,NLILL
           DO I=1,MIN(NLILL,MPT)
            XF(I)=DBLE(TPCO_XV(LILLY(I)))
            YF(I)=DBLE(TPCO_YV(LILLY(I)))
            RF(I)=TPCO_RV(LILLY(I))
            PF(I)=TPCO_PH(LILLY(I))
            WF(I)=1.D0/(DBLE(TPCO_SR(LILLY(I)))+1.D-9)
            ZF(I)=TPCO_ZV(LILLY(I))
            WZF(I)=1./(TPCO_SZ(LILLY(I))+1.E-6)
           ENDDO

           CALL TFITHL(NLILL,XF,YF,RF,PF,WF,ZF,WZF,1,
     1                 VAL,ERR,CHI2,CHISZ)
C
          CHI2 = CHI2/FLOAT(2*NLILL-5)
          LFUND =.FALSE.
          IF (CHI2.LT.CHFRE) LFUND=.TRUE.
        END IF

C
C-----> Store and disable
C
        IF (LFUND)THEN
          LB2(LOOPL)  = .TRUE.

          DO J=1,N32BITREG
           IDUM1 = ITSTU(J)
           IDUM2 = I2PAT(LOOPL,J)
           IDUM3 = IOR(IDUM1,IDUM2)
           I1PAT(ICNDU,J) = IDUM3
           I2PAT(LOOPL,J) = IDUM3
          ENDDO
C GAB:
          MY_JTLNUP(IBOT)=ITOP
          MY_JTLNDN(ITOP)=IBOT
        END IF
   40 CONTINUE
C
  100 CONTINUE

c-----------------------------------------------------------------

C
C-------------------------------------------------------------------C
C              SELECT CANDIDATES FOR TOP LINKS                      C
C              *******************************                      C
C-------------------------------------------------------------------C
C-----> Find seperately candidates for left and right parts
C       Criteria:  curvature positive
C                  max radius less than max chamber radiusc
C                  no upward link found
C
      NCNDL = 0
c GAB:
c      DO 110 I=1,LROWS(KTCHA)
      DO 110 I=1,MY_NCHAIN
        IF (NCNDL.EQ.MCAND)THEN
         IF(TKSTDBG.GT.0) WRITE(*,*)'TFLNKS:MCAND too small'
         CALL TKVETO(1)
         GO TO 110
        ENDIF
C GAB:
c        IF (RTABL(KTCHA,I,JTCHIR).LE.EPS)          GO TO 110
        IF (MY_JTCHIR(I).LE.EPS)          GO TO 110
C GAB:
C         RMTRK = 2.0/RTABL(KTCHA,I,JTCHIR) - RTABL(KTCHA,I,JTCHD0)
         RMTRK = 2.0/MY_JTCHIR(I) - MY_JTCHD0(I)

        IF (RMTRK.GT.TPCACRO)                       GO TO 110
C GAB:
c        IF (ITABL(KTLNK,I,JTLNUP).NE.0)            GO TO 110

        IF (MY_JTLNUP(I).NE.0)            GO TO 110

C GAB:
        NCNDL = NCNDL+1
        IF(NCNDL.GT.MCAND)THEN
         IF(TKSTDBG.GT.0) WRITE(*,*)'TFLNKS:MCAND too small'
         CALL TKVETO(1)
         IER=1
         RETURN
        ENDIF

        IBUF1(NCNDL) = I
  110 CONTINUE
      IF (NCNDL.LE.0)                              GO TO 200
C
C---> IBUF1 now contains the chains which might be the left
C     part of a spiral
c
c-----------------------------------------------------------------
c GAB: Join any such chains which have consistent parameters:
      DO I=1,NCNDL-1
111    JCH1=IBUF1(I)
       IF(JCH1.GT.0)THEN
        DO J=I+1,NCNDL
         JCH2=IBUF1(J)
         IF(JCH2.GT.0)THEN
          CHI2 = (MY_JTCHIR(JCH1)
     1                        -MY_JTCHIR(JCH2))**2
     2          *EOMES(JCH1)*EOMES(JCH2)/(EOMES(JCH1)+EOMES(JCH2))
     3          + (MY_JTCHTL(JCH1)
     4                        -MY_JTCHTL(JCH2))**2
     5          *EALPS(JCH1)*EALPS(JCH2)/(EALPS(JCH1)+EALPS(JCH2))
     6          + DPHI**2
     7          *EPH0S(JCH1)*EPH0S(JCH2)/(EPH0S(JCH1)+EPH0S(JCH2))
     8          + (MY_JTCHD0(JCH1)
     9                        -MY_JTCHD0(JCH2))**2
     1          *EDD0S(JCH1)*EDD0S(JCH2)/(EDD0S(JCH1)+EDD0S(JCH2))
     2          + (MY_JTCHZ0(JCH1)
     3                        -MY_JTCHZ0(JCH2))**2
     4          *EZZ0S(JCH1)*EZZ0S(JCH2)/(EZZ0S(JCH1)+EZZ0S(JCH2))

C
          CHI2 = 0.2*CHI2
          IF(CHI2.LT.CHMAX)THEN
C Get minimum distance between start of one chain and end of the other:
c order the combined chain accordingly
           IH1BEG=MY_TCTC(MY_JTCHOF(JCH1)+1)
           IH1END=MY_TCTC(MY_JTCHOF(JCH1)+MY_JTCHNC(JCH1))
           IH2BEG=MY_TCTC(MY_JTCHOF(JCH2)+1)
           IH2END=MY_TCTC(MY_JTCHOF(JCH2)+MY_JTCHNC(JCH2))
           D1=TPCO_RV(IH2BEG)-TPCO_RV(IH1END)
           D2=TPCO_RV(IH1BEG)-TPCO_RV(IH2END)
           IF((D1.LT.3.).OR.(D2.LT.3.))THEN
C Create a combined track:
            NLILL=0
            DO K=1,MY_JTCHNC(JCH1)
             NLILL=NLILL+1
             LILLY(NLILL)=MY_TCTC(MY_JTCHOF(JCH1)+K)
            ENDDO
            DO K=1,MY_JTCHNC(JCH2)
             NLILL=NLILL+1
             LILLY(NLILL)=MY_TCTC(MY_JTCHOF(JCH2)+K)
            ENDDO
C ripple sort, ordered by radius:
            KMIN=0
            DO WHILE(KMIN.LT.NLILL)
             RMIN=99999.
             KMIN=KMIN+1
             DO K=KMIN,NLILL
              IF(TPCO_RV(LILLY(K)).LT.RMIN)THEN
               ISWAP=K
               RMIN=TPCO_RV(LILLY(K))
              ENDIF
             ENDDO
             IHIT=LILLY(KMIN)
             LILLY(KMIN)=LILLY(ISWAP)
             LILLY(ISWAP)=IHIT
            ENDDO
c(KH)       DO K=1,NLILL
            DO K=1,MIN(NLILL,MPT)
             XF(K)=DBLE(TPCO_XV(LILLY(K)))
             YF(K)=DBLE(TPCO_YV(LILLY(K)))
             RF(K)=TPCO_RV(LILLY(K))
             PF(K)=TPCO_PH(LILLY(K))
             WF(K)=1.D0/(DBLE(TPCO_SR(LILLY(K)))+1.D-9)
             ZF(K)=TPCO_ZV(LILLY(K))
             WZF(K)=1./(TPCO_SZ(LILLY(K))+1.E-6)
            ENDDO
            CALL TFITHL(NLILL,XF,YF,RF,PF,WF,ZF,WZF,2,
     &                  VAL,ERR,CHI2,CHISZ)
C
            CHI2 = CHI2/FLOAT(2*NLILL-5)
            LFUND =.FALSE.
            IF (CHI2.LT.CHFRE)THEN
C----
c Re-arrange the TCTC bank to accomodate the merged chains:
             IBUF1(J)=0
             IF(JCH1.GT.JCH2)THEN
              ILAST=JCH1
              IFRST=JCH2
             ELSE
              ILAST=JCH2
              IFRST=JCH1
             ENDIF
             DO ICH=ILAST-1,IFRST+1,-1
              DO IHIT=MY_JTCHNC(ICH),1,-1
               MY_TCTC(MY_JTCHOF(ICH)+MY_JTCHNC(JCH2)+IHIT)=
     &          MY_TCTC(MY_JTCHOF(ICH)+IHIT)
             ENDDO
              MY_JTCHOF(ICH)=MY_JTCHOF(ICH)+MY_JTCHNC(JCH2)
             ENDDO
C Delete merged chain:
             MY_JTCHNC(JCH2)=0
             MY_JTCHNC(JCH2)=0
             DO K=1,15
              MY_JTCHER(JCH2,K)=0.
             ENDDO
             MY_JTCHIR(JCH2)=0.
             MY_JTCHTL(JCH2)=0.
             MY_JTCHP0(JCH2)=0.
             MY_JTCHD0(JCH2)=0.
             MY_JTCHZ0(JCH2)=0.
             MY_JTCHLC(JCH2)=0.
             MY_JTCHC1(JCH2)=0.
             MY_JTCHC2(JCH2)=0.

             MY_JTCHNC(JCH1)=NLILL
             DO IHIT=1,NLILL
              MY_TCTC(MY_JTCHOF(JCH1)+IHIT)=LILLY(IHIT)
             ENDDO
             DO K=1,15
              MY_JTCHER(JCH1,K)=ERR(K)
             ENDDO
             MY_JTCHIR(JCH1)=VAL(1)
             MY_JTCHTL(JCH1)=VAL(2)
             MY_JTCHP0(JCH1)=VAL(3)
             MY_JTCHD0(JCH1)=VAL(4)
             MY_JTCHZ0(JCH1)=VAL(5)
C             MY_JTCHLC(JCH1)=SCIRCL2 - SCIRCL1
             MY_JTCHC1(JCH1)=CHI2/FLOAT(NLILL-3)
             MY_JTCHC2(JCH1)=CHISZ/FLOAT(NLILL-2)
c Make sure we don't skip any potential links:
             GO TO 111
            ENDIF
           ENDIF
          ENDIF
         ENDIF
        ENDDO
       ENDIF
      ENDDO

C(KH) clean up buffer
      I=1
 12   IF (IBUF1(I).EQ.0) THEN
        DO J=I,NCNDL-1
          IBUF1(J)=IBUF1(J+1)
        END DO
        NCNDL=NCNDL-1
      ELSE
        I=I+1
      END IF
      IF (I.LE.NCNDL) GOTO 12

c-----------------------------------------------------------------
      NCNDR = 0
C GAB:
c      DO 120 I=1,LROWS(KTCHA)
      DO 120 I=1,MY_NCHAIN

        IF (NCNDR.EQ.MCAND)THEN
         IF(TKSTDBG.GT.0) WRITE(*,*)'TFLNKS:MCAND too small'
         CALL TKVETO(1)
         GO TO 120
        ENDIF

C GAB:
c        IF (RTABL(KTCHA,I,JTCHIR).GE.-EPS)         GO TO 120
        IF (MY_JTCHIR(I).GE.-EPS)         GO TO 120

C GAB:
c        RMTRK=-2.0/RTABL(KTCHA,I,JTCHIR) - RTABL(KTCHA,I,JTCHD0)
        RMTRK=-2.0/MY_JTCHIR(I) - MY_JTCHD0(I)

        IF (RMTRK.GT.TPCACRO)                       GO TO 120
C GAB:
c        IF (ITABL(KTLNK,I,JTLNUP).NE.0)            GO TO 120
        IF (MY_JTLNUP(I).NE.0)            GO TO 120
        NCNDR= NCNDR+1

C GAB:
        IF(NCNDR.GT.MCAND)THEN
         IF(TKSTDBG.GT.0) WRITE(*,*)'TFLNKS:MCAND too small'
         CALL TKVETO(1)
         IER=1
         RETURN
        ENDIF

        IBUF2(NCNDR)= I
        LB2(NCNDR)= .FALSE.
  120 CONTINUE
      IF (NCNDR.LE.0)                              GO TO 200


C---------------------------------------------------------------
c GAB: Join any such chains which have consistent parameters:
      DO I=1,NCNDR-1
222    JCH1=IBUF2(I)
       IF(JCH1.GT.0)THEN
        DO J=I+1,NCNDR
         JCH2=IBUF2(J)
         IF(JCH2.GT.0)THEN
          CHI2 = (MY_JTCHIR(JCH1)
     1                         -MY_JTCHIR(JCH2))**2
     2           *EOMES(JCH1)*EOMES(JCH2)/(EOMES(JCH1)+EOMES(JCH2))
     3           + (MY_JTCHTL(JCH1)
     4                         -MY_JTCHTL(JCH2))**2
     5           *EALPS(JCH1)*EALPS(JCH2)/(EALPS(JCH1)+EALPS(JCH2))
     6           + DPHI**2
     7           *EPH0S(JCH1)*EPH0S(JCH2)/(EPH0S(JCH1)+EPH0S(JCH2))
     8           + (MY_JTCHD0(JCH1)
     9                         -MY_JTCHD0(JCH2))**2
     1           *EDD0S(JCH1)*EDD0S(JCH2)/(EDD0S(JCH1)+EDD0S(JCH2))
     2           + (MY_JTCHZ0(JCH1)
     3                         -MY_JTCHZ0(JCH2))**2
     4           *EZZ0S(JCH1)*EZZ0S(JCH2)/(EZZ0S(JCH1)+EZZ0S(JCH2))
C
          CHI2 = 0.2*CHI2
          IF(CHI2.LT.CHMAX)THEN
C Get minimum distance between start of one chain and end of the other:
c order the combined chain accordingly
           IH1BEG=MY_TCTC(MY_JTCHOF(JCH1)+1)
           IH1END=MY_TCTC(MY_JTCHOF(JCH1)+MY_JTCHNC(JCH1))
           IH2BEG=MY_TCTC(MY_JTCHOF(JCH2)+1)
           IH2END=MY_TCTC(MY_JTCHOF(JCH2)+MY_JTCHNC(JCH2))
           D1=TPCO_RV(IH2BEG)-TPCO_RV(IH1END)
           D2=TPCO_RV(IH1BEG)-TPCO_RV(IH2END)
           IF((D1.LT.3.).OR.(D2.LT.3.))THEN
C Create a combined track:
            NLILL=0
            DO K=1,MY_JTCHNC(JCH1)
             NLILL=NLILL+1
             LILLY(NLILL)=MY_TCTC(MY_JTCHOF(JCH1)+K)
            ENDDO
            DO K=1,MY_JTCHNC(JCH2)
             NLILL=NLILL+1
             LILLY(NLILL)=MY_TCTC(MY_JTCHOF(JCH2)+K)
            ENDDO
C ripple sort, ordered by radius:
            KMIN=0
            DO WHILE(KMIN.LT.NLILL)
             RMIN=99999.
             KMIN=KMIN+1
             DO K=KMIN,NLILL
              IF(TPCO_RV(LILLY(K)).LT.RMIN)THEN
               ISWAP=K
               RMIN=TPCO_RV(LILLY(K))
              ENDIF
             ENDDO
             IHIT=LILLY(KMIN)
             LILLY(KMIN)=LILLY(ISWAP)
             LILLY(ISWAP)=IHIT
            ENDDO
c(KH)       DO K=1,NLILL
            DO K=1,MIN(NLILL,MPT)
             XF(K)=DBLE(TPCO_XV(LILLY(K)))      !(KH) K>200 -> crash
             YF(K)=DBLE(TPCO_YV(LILLY(K)))
             RF(K)=TPCO_RV(LILLY(K))
             PF(K)=TPCO_PH(LILLY(K))
             WF(K)=1.D0/(DBLE(TPCO_SR(LILLY(K)))+1.D-9)
             ZF(K)=TPCO_ZV(LILLY(K))
             WZF(K)=1./(TPCO_SZ(LILLY(K))+1.E-6)
            ENDDO
            CALL TFITHL(NLILL,XF,YF,RF,PF,WF,ZF,WZF,2,
     &                  VAL,ERR,CHI2,CHISZ)
C
            CHI2 = CHI2/FLOAT(2*NLILL-5)
            LFUND =.FALSE.
            IF (CHI2.LT.CHFRE)THEN
C----
c Re-arrange the TCTC bank to accomodate the merged chains:
             IBUF2(J)=0
             IF(JCH1.GT.JCH2)THEN
              ILAST=JCH1
              IFRST=JCH2
             ELSE
              ILAST=JCH2
              IFRST=JCH1
             ENDIF
             DO ICH=ILAST-1,IFRST+1,-1
              DO IHIT=MY_JTCHNC(ICH),1,-1
               MY_TCTC(MY_JTCHOF(ICH)+MY_JTCHNC(JCH2)+IHIT)=
     &          MY_TCTC(MY_JTCHOF(ICH)+IHIT)
             ENDDO
              MY_JTCHOF(ICH)=MY_JTCHOF(ICH)+MY_JTCHNC(JCH2)
             ENDDO
C Delete merged chain:
             MY_JTCHNC(JCH2)=0
             DO K=1,15
              MY_JTCHER(JCH2,K)=0.
             ENDDO
             MY_JTCHIR(JCH2)=0.
             MY_JTCHTL(JCH2)=0.
             MY_JTCHP0(JCH2)=0.
             MY_JTCHD0(JCH2)=0.
             MY_JTCHZ0(JCH2)=0.
             MY_JTCHLC(JCH2)=0.
             MY_JTCHC1(JCH2)=0.
             MY_JTCHC2(JCH2)=0.
             MY_JTCHNC(JCH1)=NLILL
             DO IHIT=1,NLILL
              MY_TCTC(MY_JTCHOF(JCH1)+IHIT)=LILLY(IHIT)
             ENDDO
             DO K=1,15
              MY_JTCHER(JCH1,K)=ERR(K)
             ENDDO
             MY_JTCHIR(JCH1)=VAL(1)
             MY_JTCHTL(JCH1)=VAL(2)
             MY_JTCHP0(JCH1)=VAL(3)
             MY_JTCHD0(JCH1)=VAL(4)
             MY_JTCHZ0(JCH1)=VAL(5)
C             MY_JTCHLC(JCH1)=SCIRCL2 - SCIRCL1
             MY_JTCHC1(JCH1)=CHI2/FLOAT(NLILL-3)
             MY_JTCHC2(JCH1)=CHISZ/FLOAT(NLILL-2)
c Make sure we don't skip any potential links:
             GO TO 222
            ENDIF
           ENDIF
          ENDIF
         ENDIF
        ENDDO
       ENDIF
      ENDDO

C(KH) clean up buffer
      I=1
 22   IF (IBUF2(I).EQ.0) THEN
        DO J=I,NCNDR-1
          IBUF2(J)=IBUF2(J+1)
          LB2(J)=LB2(J+1)
        END DO
        NCNDR=NCNDR-1
      ELSE
        I=I+1
      END IF
      IF (I.LE.NCNDR) GOTO 22

C
C-------------------------------------------------------------------C
C              DEFINE TOP LINKS BETWEEN CANDIDATES                  C
C              ***********************************                  C
C-------------------------------------------------------------------C
C
      DO 140 ICNDL = 1,NCNDL
        K1     = IBUF1(ICNDL)
        CHX    = CHTMX
        DO 130 ICNDR = 1,NCNDR
          IF (LB2(ICNDR))                               GO TO 130
          K2     = IBUF2(ICNDR)
C(KH)----bug workaround
          IF ((EALP(K1).EQ.0.0).OR.(EALP(K2).EQ.0.0)) THEN
             IF(TKSTDBG.GT.0) THEN
             WRITE(6,*) 'TFLNKS: internal error. k1,k2,ealp(k1),'//
     >                  'ealp(k2)=',k1,k2,ealp(k1),ealp(k2)
             ENDIF
             CALL TKVETO(1)
             GO TO 130
          END IF
C(KH)----end
          CHI2 = (ALP(K1)+ALP(K2))**2
     1         /( (1.0+CHF(K1))/EALP(K1) +(1.0+CHF(K2))/EALP(K2) )
C
          IF (CHI2.GT.(5.0*CHTMX))                       GO TO 130
C
          DPHI = AMOD(PH0(K1)-PH0(K2)-PI,TWOPI)
          IF (DPHI.LT.-5.0)DPHI=DPHI+TWOPI
          CHI2 = CHI2
     1      +   ( 2.0*(1.0/OME(K1)+1.0/OME(K2))-DD0(K1)+DD0(K2) )**2
     2        /( 4.0*(  1.0/(EOME(K1)*OME(K1)*OME(K1))
     3               +1.0/(EOME(K1)*OME(K2)*OME(K2)) )
     4        +1.0/EDD0(K1) + 1.0/EDD0(K2)               )
     5      +    DPHI**2
     6        /((2.0+CHF(K1))/EPH0(K1) + (2.0+CHF(K2))/EPH0(K2) )
     7    +   (ZZ0(K1)-ZZ0(K2)+PI*(ALP(K1)/OME(K1)+ALP(K2)/OME(K2)))**2
     8        /(   2.0/EZZ0(K1) + 2.0/EZZ0(K2)   )
          CHI2 = 0.25*CHI2

          IF (CHI2.LT.CHX) THEN
            LOOPR = ICNDR
            CHX = CHI2
          END IF
  130   CONTINUE

        IF (CHX.LT.CHTMX) THEN
          LB2(LOOPR) = .TRUE.
C GAB:
c          IW(KROW(KTLNK,IBUF1(ICNDL))+JTLNTP)=IBUF2(LOOPR)
c          IW(KROW(KTLNK,IBUF2(LOOPR))+JTLNTP)=IBUF1(ICNDL)
          MY_JTLNTP(IBUF1(ICNDL))=IBUF2(LOOPR)
          MY_JTLNTP(IBUF2(LOOPR))=IBUF1(ICNDL)
        END IF
  140 CONTINUE
C
  200 CONTINUE

C
C-------------------------------------------------------------------C
C              SELECT CANDIDATES FOR BOTTOM LINKS                   C
C                      WITHIN THE T.P.C.                            C
C              **********************************                   C
C-------------------------------------------------------------------C
C-----> Find seperately candidates for left and right parts
C       Criteria:  curvature positive
C                  minimum distance to beam line = R(ITC)
C                  no downward link found
C
      NCNDL = 0
C GAB:
C      DO 210 I=1,LROWS(KTCHA)
      DO 210 I=1,MY_NCHAIN
        IF (NCNDL.EQ.MCAND)THEN
         IF(TKSTDBG.GT.0) WRITE(*,*)'TFLNKS:MCAND too small'
         CALL TKVETO(1)
         GO TO 210
        ENDIF
C GAB:
C        IF (RTABL(KTCHA,I,JTCHIR).LE.EPS)          GO TO 210

        IF (MY_JTCHIR(I).LE.EPS)          GO TO 210
C GAB:
C        RMINI = ABS(RTABL(KTCHA,I,JTCHD0))
        RMINI = ABS(MY_JTCHD0(I))

        IF (RMINI.LT.TPPROW(1))                    GO TO 210

C GAB:
c        IF (ITABL(KTLNK,I,JTLNDN).NE.0)            GO TO 210
        IF (MY_JTLNDN(I).NE.0)            GO TO 210
        NCNDL = NCNDL+1

C GAB:
        IF(NCNDL.GT.MCAND)THEN
         IF(TKSTDBG.GT.0) WRITE(*,*)'TFLNKS:MCAND too small'
         CALL TKVETO(1)
         IER=1
         RETURN
        ENDIF

        IBUF1(NCNDL) = I
  210 CONTINUE

      IF (NCNDL.LE.0)                              GO TO 300
C
      NCNDR = 0
C GAB:
c      DO 220 I=1,LROWS(KTCHA)
      DO 220 I=1,MY_NCHAIN
        IF (NCNDR.EQ.MCAND)                         GO TO 220
C GAB:
c        IF (RTABL(KTCHA,I,JTCHIR).GE.-EPS)          GO TO 220
        IF (MY_JTCHIR(I).GE.-EPS)          GO TO 220
C GAB:
c        RMINI = ABS(RTABL(KTCHA,I,JTCHD0))
        RMINI = ABS(MY_JTCHD0(I))

        IF (RMINI.LT.TPPROW(1))                     GO TO 220


C GAB:
c        IF (ITABL(KTLNK,I,JTLNDN).NE.0)             GO TO 220
        IF (MY_JTLNDN(I).NE.0)             GO TO 220
        NCNDR = NCNDR+1

C GAB:
        IF(NCNDR.GT.MCAND)THEN
         IF(TKSTDBG.GT.0) WRITE(*,*)'TFLNKS:MCAND too small'
         CALL TKVETO(1)
         IER=1
         RETURN
        ENDIF

        IBUF2(NCNDR) = I
        LB2(NCNDR)   = .FALSE.
  220 CONTINUE
      IF (NCNDR.LE.0)                               GO TO 300

c-----------------------------------------------------------------



C
C-------------------------------------------------------------------C
C              DEFINE BOTTOM LINKS BETWEEN CANDIDATES               C
C              **************************************               C
C-------------------------------------------------------------------C
C
      DO 240 ICNDL = 1,NCNDL
        K1     = IBUF1(ICNDL)
        CHX    = CHTMX
        DO 230 ICNDR = 1,NCNDR
          IF (LB2(ICNDR))                               GO TO 230
          K2     = IBUF2(ICNDR)
C(KH) ---- bug workaround
          IF ((K1.EQ.0).OR.(K2.EQ.0)) THEN
             IF(TKSTDBG.GT.0) THEN
                WRITE(6,*) 'TFLNKS: internal error. k1,k2=',k1,k2
             ENDIF
             CALL TKVETO(1)
             GO TO 230
          END IF
          IF ((EDD0(K1).EQ.0.0).OR.(EDD0(K2).EQ.0.0)) THEN
             IF(TKSTDBG.GT.0) THEN
                WRITE(6,*) 'TFLNKS: internal error 230. '//
     >               'k1,k2,edd0(k1),edd0(k2)='
     >               ,k1,k2,edd0(k1),edd0(k2)
                ENDIF
             CALL TKVETO(1)
             GO TO 230
          END IF
C(KH) ---- end
          CHI2 = (DD0(K1)-DD0(K2))**2
     1        /((1.0+CHF(K1))/EDD0(K1) + (1.0+CHF(K2))/EDD0(K2))
C
          IF (CHI2.GT.(5.0*CHTMX))                       GO TO 230
C
          DPHI = AMOD(PH0(K1)-PH0(K2)-PI,TWOPI)
          IF (DPHI.LT.-5.0)DPHI=DPHI+TWOPI
          CHI2 = CHI2
     1       +   (ALP(K1)+ALP(K2))**2
     2       /( (1.0+CHF(K1))/EALP(K1) +(1.0+CHF(K2))/EALP(K2) )
     3       +     0.2*DPHI**2
     4       /((2.0+CHF(K1))/EPH0(K1) + (2.0+CHF(K2))/EPH0(K2) )
     5       +     (ZZ0(K1)-ZZ0(K2))**2
     6       /((1.0+CHF(K1))/EZZ0(K1) + (1.0+CHF(K2))/EZZ0(K2) )
          CHI2 = 0.25*CHI2
          IF (CHI2.LT.CHX) THEN
            LOOPR = ICNDR
            CHX = CHI2
          END IF
  230   CONTINUE
        IF (CHX.LT.CHTMX) THEN
          LB2(LOOPR) = .TRUE.
C GAB:
c          IW(KROW(KTLNK,IBUF1(ICNDL))+JTLNBM)=IBUF2(LOOPR)
c          IW(KROW(KTLNK,IBUF2(LOOPR))+JTLNBM)=IBUF1(ICNDL)
          MY_JTLNBM(IBUF1(ICNDL))=IBUF2(LOOPR)
          MY_JTLNBM(IBUF2(LOOPR))=IBUF1(ICNDL)

        END IF
  240 CONTINUE
C
  300 CONTINUE
C
C-------------------------------------------------------------------C
C              SELECT CANDIDATES FOR BOTTOM LINKS                   C
C                EXTRAPOLATING INTO THE I.T.C.                      C
C              **********************************                   C
C-------------------------------------------------------------------C
c-----> Find seperately candidates for left and right parts
C       Criteria:  curvature positive
C                  2* radius < TPC radius
C                  minimum distance to beam line < R(ITC)
C                  no downward link found
C
      NCNDL = 0
C GAB:
C      DO 310 I=1,LROWS(KTCHA)
      DO 310 I=1,MY_NCHAIN
        IF (NCNDL.EQ.MCAND)THEN
           IF(TKSTDBG.GT.0) WRITE(*,*)'TFLNKS:MCAND too small'
         CALL TKVETO(1)
         GO TO 310
        ENDIF
C GAB:
C        IF (RTABL(KTCHA,I,JTCHIR).LE.EPS)           GO TO 310
        IF (MY_JTCHIR(I).LE.EPS)           GO TO 310
        RLRGE = RTOP(I)
        IF (RLRGE.GT.TPCACRO)                        GO TO 310
C GAB:
C        RMINI = ABS(RTABL(KTCHA,I,JTCHD0))
        RMINI = ABS(MY_JTCHD0(I))
        IF (RMINI.GE.TPPROW(1))                     GO TO 310

C GAB:
C        IF (ITABL(KTLNK,I,JTLNDN).NE.0)             GO TO 310
        IF (MY_JTLNDN(I).NE.0)             GO TO 310
        NCNDL = NCNDL+1

C GAB:
        IF(NCNDL.GT.MCAND)THEN
         IF(TKSTDBG.GT.0) WRITE(*,*)'TFLNKS:MCAND too small'
         CALL TKVETO(1)
         IER=1
         RETURN
        ENDIF

        IBUF1(NCNDL) = I
  310 CONTINUE
      IF (NCNDL.LE.0)                               GO TO 400
C
      NCNDR = 0
C GAB:
c      DO 320 I=1,LROWS(KTCHA)
      DO 320 I=1,MY_NCHAIN
        IF (NCNDR.EQ.MCAND)                        GO TO 320
C GAB:
c        IF (RTABL(KTCHA,I,JTCHIR).GE.-EPS)         GO TO 320
        IF (MY_JTCHIR(I).GE.-EPS)         GO TO 320
        RLRGE = RTOP(I)
        IF (RLRGE.GT.TPCACRO)                       GO TO 320

C        RMINI = ABS(RTABL(KTCHA,I,JTCHD0))

        RMINI = ABS(MY_JTCHD0(I))
        IF (RMINI.GE.TPPROW(1))                    GO TO 320

C GAB:
C        IF (ITABL(KTLNK,I,JTLNDN).NE.0)            GO TO 320
        IF (MY_JTLNDN(I).NE.0)            GO TO 320
        NCNDR = NCNDR+1

C GAB:
        IF(NCNDR.GT.MCAND)THEN
         IF(TKSTDBG.GT.0) WRITE(*,*)'TFLNKS:MCAND too small'
         CALL TKVETO(1)
         IER=1
         RETURN
        ENDIF

        IBUF2(NCNDR)= I
        LB2(NCNDR)= .FALSE.
  320 CONTINUE
      IF (NCNDR.LE.0)                              GO TO 400


C-------------------------------------------------------------------C
C              DEFINE BOTTOM LINKS BETWEEN CANDIDATES               C
C              **************************************               C
C-------------------------------------------------------------------C
C
      DO 340 ICNDL = 1,NCNDL
        K1 = IBUF1(ICNDL)
        CHX = CHPMX
        DO 330 ICNDR = 1,NCNDR
          IF (LB2(ICNDR))                               GO TO 330
          K2 = IBUF2(ICNDR)
C(KH)----bug workaround
          IF ((EDD0(K1).EQ.0.0).OR.(EDD0(K2).EQ.0.0)) THEN
             IF(TKSTDBG.GT.0) THEN
                WRITE(6,*) 'TFLNKS: internal error 330. '//
     >               'k1,k2,edd0(k1),edd0(k2)='
     >               ,k1,k2,edd0(k1),edd0(k2)
             ENDIF
             CALL TKVETO(1)
             GO TO 330
          END IF
C(KH)----end
          CHI2 = 0.2 * (DD0(K1)-DD0(K2))**2
     1         /((1.0+CHF(K1))/EDD0(K1) + (1.0+CHF(K2))/EDD0(K2))
C
          IF (CHI2.GT.(5.0*CHTMX))                       GO TO 330
C
          DPHI = AMOD(PH0(K1)-PH0(K2)-PI,TWOPI)
          IF (DPHI.LT.-5.0)DPHI=DPHI+TWOPI
          CHI2 = CHI2
     1         +   (ALP(K1)+ALP(K2))**2
     2         /( (1.0+CHF(K1))/EALP(K1) +(1.0+CHF(K2))/EALP(K2) )
     3         +     0.2*DPHI**2
     4         /((2.0+CHF(K1))/EPH0(K1) + (2.0+CHF(K2))/EPH0(K2) )
     5         +    (ZZ0(K1)-ZZ0(K2))**2
     6         /((1.0+CHF(K1))/EZZ0(K1) + (1.0+CHF(K2))/EZZ0(K2) )
          CHI2 = 0.25*CHI2

          IF (CHI2.LT.CHX) THEN
            LOOPR = ICNDR
            CHX = CHI2
          END IF
  330   CONTINUE
        IF (CHX.LT.CHPMX) THEN
          LB2(LOOPR) = .TRUE.
C GAB:
c          IW(KROW(KTLNK,IBUF1(ICNDL))+JTLNBM)=IBUF2(LOOPR)
c          IW(KROW(KTLNK,IBUF2(LOOPR))+JTLNBM)=IBUF1(ICNDL)
          MY_JTLNBM(IBUF1(ICNDL))=IBUF2(LOOPR)
          MY_JTLNBM(IBUF2(LOOPR))=IBUF1(ICNDL)
        END IF
  340 CONTINUE
C
  400 CONTINUE
C
C-------------------------------------------------------------------C
C       Select Candidates for "long range + high pt" Links          C
C       ( for example muons passing sector boundaries )             C
C       ***************************************************         C
C-------------------------------------------------------------------C
C-----> Find separately candidates for top and bottom parts.
C       Criteria:  Curvature  > 5*R_tpc
C                  Only hits in the outer 6 pad rows
C                  no link found yet
C
      NCNDU = 0
C GAB:
C      DO 410 I=1,LROWS(KTCHA)
      DO 410 I=1,MY_NCHAIN
        IF (NCNDU.EQ.MCAND)THEN
         IF(TKSTDBG.GT.0) WRITE(*,*)'TFLNKS: MCAND too small'
         CALL TKVETO(1)
         GO TO 410
        ENDIF
C GAB:
c        IF (ABS(RTABL(KTCHA,I,JTCHIR)).GT.OCUT)      GO TO 410
        IF (ABS(MY_JTCHIR(I)).GT.OCUT)      GO TO 410
C GAB:
C        ITPTES = ITABL(KTCHA,I,JTCHPI)
        ITPTES = MY_JTCHPI(I,NRTPC/32+1)
        ISEL  = 0
        ISEL  = IAND(ITPTES,IF15)
        IF (ISEL.NE.0)                               GO TO 410
C GAB:
c        IF (ITABL(KTLNK,I,JTLNDN).NE.0)              GO TO 410
        IF (MY_JTLNDN(I).NE.0)              GO TO 410
        NCNDU = NCNDU + 1
        IBUF1(NCNDU) = I
C GAB:
C        I1PAT(NCNDU) = ITPTES
        DO J=1,N32BITREG
         I1PAT(NCNDU,J) = MY_JTCHPI(I,J)
        ENDDO
  410 CONTINUE

      IF (NCNDU.LE.0)                     GO TO 500
C
      NCNDL  = 0
C GAB:
c      DO 420 I=1,LROWS(KTCHA)
      DO 420 I=1,MY_NCHAIN
        IF (NCNDL.EQ.MCAND)THEN
         IF(TKSTDBG.GT.0) WRITE(*,*)'TFLNKS:MCAND too small'
         CALL TKVETO(1)
         GO TO 420
        ENDIF
C GAB:
c        IF (ABS(RTABL(KTCHA,I,JTCHIR)).GT.OCUT)       GO TO 420
        IF (ABS(MY_JTCHIR(I)).GT.OCUT)       GO TO 420
C GAB:
C        ITPTES = ITABL(KTCHA,I,JTCHPI)
        ITPTES = MY_JTCHPI(I,NRTPC/32+1)
        ISEL  = 0
        ISEL  = IAND(ITPTES,IL06)
        IF (ISEL.NE.0)                                GO TO 420
C GAB:
C        IF (ITABL(KTLNK,I,JTLNUP).NE.0)               GO TO 420
        IF (MY_JTLNUP(I).NE.0)               GO TO 420
C GAB:
c        LSTPT = ITABL(KTCHA,I,JTCHOF) + ITABL(KTCHA,I,JTCHNC)
        LSTPT = MY_JTCHOF(I) + MY_JTCHNC(I)
        LB1PT = LSTPT - 1
C GAB:
c        LSTCO = IW(KTCTC+2+LSTPT)
c        LB1CO = IW(KTCTC+2+LB1PT)
        LSTCO = MY_TCTC(LSTPT)
        LB1CO = MY_TCTC(LB1PT)

C GAB:
c        R1= RTABL(KTPCO,LB1CO,JTPCRV)
c        R2= RTABL(KTPCO,LSTCO,JTPCRV)
c        Z1= RTABL(KTPCO,LB1CO,JTPCZV)
c        Z2= RTABL(KTPCO,LSTCO,JTPCZV)

C        R1= TPCR(LB1CO)
C        R2= TPCR(LSTCO)
C        Z1= TPCHITS(3,LB1CO)
C        Z2= TPCHITS(3,LSTCO)
        R1= TPCO_RV(LB1CO)
        R2= TPCO_RV(LSTCO)
        Z1= TPCO_ZV(LB1CO)
        Z2= TPCO_ZV(LSTCO)
        ZEXT= ABS(Z2 + 2.0*DROW*(Z2-Z1)/(R2-R1))
        IF (ZEXT.GT.ZDRIFT)                             GO TO 420
        NCNDL = NCNDL + 1
        IF(NCNDL.GT.MCAND)THEN
         IF(TKSTDBG.GT.0) WRITE(*,*)'TFLNKS:MCAND too small'
         CALL TKVETO(1)
         IER=1
         RETURN
        ENDIF

        IBUF2(NCNDL) = I
        LB2(NCNDL)   = .FALSE.
  420 CONTINUE
      IF (NCNDL.LE.0)                      GO TO 500

c GAB: This is buggy - it does not forbid overlapping row patterns:
cC
cC-------------------------------------------------------------------C
cC              DEFINE UP/DOWN LINKS BETWEEN CANDIDATES              C
cC              ***************************************              C
cC-------------------------------------------------------------------C
cC
c      DO 440 ICNDU = 1,NCNDU
c        ITOP  = IBUF1(ICNDU)
cC GAB:
cc        NQU   = ITABL(KTCHA,ITOP,JTCHNC)
c        NQU   = MY_JTCHNC(ITOP)
c        DO 425 IPT = 1,NQU
cC GAB
cC          LILLY(IPT) = IW(LIST(ITOP)+IPT)
c          LILLY(IPT) = MY_TCTC(MY_JTCHOF(ITOP)+IPT)
c  425   CONTINUE
c        chi2_best = 1.e30
c        link_best = 0
c        DO 430 ICNDL = 1,NCNDL
cc
c
c          IBOT   = IBUF2(ICNDL)
c          IF (LB2(ICNDL))           GO TO 430
cC GAB:
cC          NQL   = ITABL(KTCHA,IBOT,JTCHNC)
c          NQL   = MY_JTCHNC(IBOT)
c          DO 426 IPT = 1,NQL
cC GAB:
cc            LILLY(NQU+IPT) = IW(LIST(IBOT)+IPT)
c          LILLY(NQU+IPT) = MY_TCTC(MY_JTCHOF(IBOT)+IPT)
c  426     CONTINUE
c          NLILL  = NQU+NQL
cC GAB:
cc          CALL TFTHEL('TPCO',NLILL,LILLY,1,VAL,ERR,CHI2)
c           DO I=1,NLILL
c            XF(I)=DBLE(TPCO_XV(LILLY(I)))
c            YF(I)=DBLE(TPCO_YV(LILLY(I)))
c            RF(I)=TPCO_RV(LILLY(I))
c            PF(I)=TPCO_PH(LILLY(I))
c            WF(I)=1.D0/(DBLE(TPCO_SR(LILLY(I)))+1.D-9)
c            ZF(I)=TPCO_ZV(LILLY(I))
c            WZF(I)=1./(TPCO_SZ(LILLY(I))+1.E-6)
c           ENDDO
c

c           CALL TFITHL(NLILL,XF,YF,RF,PF,WF,ZF,WZF,1,
c     1                 VAL,ERR,CHI2,CHISZ)
c
c          CHI2 = CHI2/FLOAT(2*NLILL-5)
c          IF (CHI2.LT.CHFRE .and. chi2.lt.chi2_best) THEN
c            chi2_best = chi2
c            link_best = icndl
c          endif
c  430   CONTINUE
cC+ Select best link rather than the first
c        if (link_best.gt.0) then
c            LB2(link_best)   = .TRUE.
c            ibot = ibuf2(link_best)
cC GAB:
cc            IW(KROW(KTLNK,IBOT)+JTLNUP) =ITOP
cc            IW(KROW(KTLNK,ITOP)+JTLNDN) =IBOT
c            MY_JTLNUP(IBOT) =ITOP
c            MY_JTLNDN(ITOP) =IBOT
c        endif
c  440 CONTINUE
C

  500 CONTINUE

      RETURN
C  998 CONTINUE
C      RETURN
C  999 CONTINUE
Cc      CALL RERROR('TFLNKS',1,'No room in BOS for named banks.')
C      IER=1
C      RETURN
      END
C----------------------------------------------------------------------------
