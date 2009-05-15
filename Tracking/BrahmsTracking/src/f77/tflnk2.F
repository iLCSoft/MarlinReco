*  
* $Id: tflnk2.F,v 1.6 2009-05-15 15:12:30 aplin Exp $
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
* Revision 1.4  2007/09/05 09:47:29  rasp
* Updated version
*
* Revision 1.2  2005/07/29 15:06:17  aplin
* Initialised number of tracks found to 0 in tpcrun.
*
* Revision 1.1.1.1  2005/07/11 15:15:48  aplin
* Initial version
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
      SUBROUTINE TFLNK2(IER)
C****************************************************************************
c* Version number: 0.0      Date of last modification  26.08.99             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
C* Main modifiers from ALEPH to BRAHMS code:                                *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
C-----------------------------------------------------
C! Finds links between pairs of tracks in the TPC.
C  If the chi**2 of the combined fit is too bad the point
C  with the largest contribution to the chi**2 is rejected
C  and the fit is repeated. A maximum of 2 points can be removed.
C
C  Input banks:  T1FT, T1TL, T1CL
C  Output banks: TGFT, TGTL, TGCL
C
C  Author : LL. Garrido 30/5/1989
C  Modified:  R. Johnson   29-05-91
C  Modified:  P. Colrain   01-08-95
C             Treatment of track candidates (2 linked tracks)
C             is now consistent with that in TFITTK, the main
C             TPC track fitting routine
C             D. Casper    09-07-97
c             Check for bogus fits
C-----------------------------------------------------
      IMPLICIT NONE

      INTEGER IER
#include "include/ufitcb.inc"
c#include "fconjj.h"
C
c----------------------------------------------------------------------
c my stuff:
#include "include/tpgeom.inc"
#include "include/fndparms.inc"
#include "include/tpco.inc"
#include "include/my_tcal.inc"
#include "include/my_t1ft.inc"
#include "include/my_tgft.inc"
#include "include/my_t1tl.inc"
#include "include/my_tgtl.inc"
#include "include/tpcdims.inc"
#include "include/my_t1cl.inc"
#include "include/my_tgcl.inc"
c#include "gcflag.inc"
#include "include/tflnwb.inc"
#include "include/alcons.inc"
#include "include/tpcmatpar.inc"
#include "include/my_tcha.inc"
C
C---> Steering constants follow:
C
#include "include/strcons.inc"

C-------------------------------------
C GAB:tmp:
c#include "include/geantrk.inc"
C(KH) INTEGER IGNT(100)
c      INTEGER IGNT(MAXCND+1)
C--------------------------------------

      INTEGER INDEX
      INTEGER MY_NBMWRK(1000)
      INTEGER N1A,N1B,N2A,N2B
CC GAB NEW:
C      INTEGER L1,L2,L3,N3A,N3B,NTHRSM,ITHRSM(100,3)
c----------------------------------------------------------------------
C
c GAB:
      INTEGER MXGD
      REAL VV1,VV2,VV0,EE0,COV
      PARAMETER (MXGD=400)
c(kh)      PARAMETER (MXGD=200)
      DIMENSION VV1(10),VV2(10),VV0(6),EE0(21),COV(21)
C GAB
C      DIMENSION LIST1(21),LIST2(2),LIST3(2)
      INTEGER LIST1,LIST2,LIST3
      DIMENSION LIST1(LTPDRO),LIST2(2),LIST3(2)
C GAB:
C      DIMENSION LGOOD(MXGD,27),QGOOD(MXGD,21)
      INTEGER LGOOD
      REAL QGOOD
c(kh) DIMENSION LGOOD(MXGD,127),QGOOD(MXGD,22)
      DIMENSION LGOOD(MXGD,LTPDRO),QGOOD(MXGD,22)
C GAB:
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      REAL RF
      DATA RF/75./
C


c GAB:
      INTEGER IX
      INTEGER MAXGDHT,I,J,K,NTRI,KP1,NPT,NPT1,NPT2,ICO,NPAD,NG,I1,I2
      INTEGER IFF,IT1F1,IT1F2,IK,J1,J2,IJ,LOST,IOPT
      INTEGER NITC,NVDET,NPTPC1,NPTPC2,NTPC,NMULT,IERR,IMAX,IKK
      INTEGER IT1,IT2,ICH1,ICH2,NPREG,NTASS,IT,NT,I0,NFL
      INTEGER ITGFT,ITGTL,NRT1,NRT2,NRT,KP2,IT1FT,IT1TL,KP,NTOT,IC
      INTEGER NCL,ITK,IOFF,N1,II
      REAL CHIP,DPHI,CHI2,PPP,CMAX

      INTEGER JBIT
      EXTERNAL JBIT
      REAL PROB
      EXTERNAL PROB

c      INTEGER NMAXMLT,NMXINMT
c      PARAMETER(NMAXMLT=10)

c      PARAMETER(NMXINMT=10)
c      INTEGER NMLNK(NMAXMLT),MULTLNK(NMAXMLT,NMXINMT)
c      INTEGER IMULTI,ILINK
c      LOGICAL LNEW,LTHCND(3)

      SAVE



C-------------------------------------------------------------------------

      IF(FIRST) THEN
        FIRST=.FALSE.
        MAXGDHT=6+NRTPC
C        NBPWRK=0
c        NTPCO=NAMIND('TPCO')
c        NFCON=NAMIND('FCON')
c        NATGFT=NAMIND('TGFT')
c        NATGTL=NAMIND('TGTL')
c        NATGCL=NAMIND('TGCL')
c        NT1FT=NAMIND('T1FT')
c        NT1TL=NAMIND('T1TL')
c        NT1CL=NAMIND('T1CL')
      ENDIF
      DO I=1,1000
       DO J=1,N32BITREG
        NBPWRK(I,J)=0
       ENDDO
      ENDDO

cC
cC++   Check if the input banks exist
cC
c      KT1FT=IW(NT1FT)
c      IF(KT1FT.EQ.0) THEN
c        IER=0
c        RETURN
c      ENDIF
c GAB:
c      NTRI=LROWS(KT1FT)
      NTRI=NT1FT

      IF(NTRI.EQ.0)THEN
       WRITE(*,*)'TRLNK2: no track candidates'
       RETURN
      ENDIF
c      KT1TL=IW(NT1TL)
c      KT1CL=IW(NT1CL)
c      IF (KT1TL.EQ.0 .OR. KT1CL.EQ.0) THEN
c        IER=5
c        CALL RERROR('TFLNK2',IER,'T1TL and T1CL banks are missing.')
c        RETURN
c      ENDIF
c      KTPCO=IW(NTPCO)
c      IF (KTPCO.EQ.0) THEN
c        IER=6
c        CALL RERROR('TFLNK2',IER,'TPCO coordinate banks is missing.')
c        RETURN
c      ENDIF
cC
cC++   Drop the output banks if they already exist
cC
c      KTGFT=NDROP('TGFT',1)
c      KTGTL=NDROP('TGTL',1)
c      KTGCL=NDROP('TGCL',1)
cC
cC++   Create the output banks
cC
c      LEN=NTRI*LTGFTA+LMHLEN
c      CALL AUBOS('TGFT',1,LEN,KTGFT,IGARB)
c      IF (IGARB.EQ.2) THEN
c        IER=8
c        CALL RERROR('TFLNK2',IER,'No room in BOS for bank TGFT')
c        RETURN
c      ENDIF
c      IW(KTGFT+LMHCOL)=LTGFTA
c      LEN=NTRI*LTGTLA+LMHLEN
c      CALL AUBOS('TGTL',1,LEN,KTGTL,IGARB)
c      IF (IGARB.EQ.2) THEN
c        IER=9
c        KTGFT=NDROP('TGFT',1)
c        CALL RERROR('TFLNK2',IER,'No room in BOS for bank TGTL')
c        RETURN
c      ENDIF
c      IF (IGARB.NE.0) THEN
c        KTGFT=NLINK('TGFT',1)
c      ENDIF
c      IW(KTGTL+LMHCOL)=LTGTLA
c      KT1CL=IW(NT1CL)
c      LEN=IW(KT1CL)
c      CALL AUBOS('TGCL',1,LEN,KTGCL,IGARB)
c      IF (IGARB.EQ.2) THEN
c        IER=10
c        KTGFT=NDROP('TGFT',1)
c        KTGTL=NDROP('TGTL',1)
c        CALL RERROR('TFLNK2',IER,'No room in BOS for bank TGCL')
c        RETURN
c      ENDIF
c      IF (IGARB.NE.0) THEN
c        KTGFT=NLINK('TGFT',1)
c        KTGTL=NLINK('TGTL',1)
c      ENDIF
c      IW(KTGCL+LMHCOL)=LTGCLA
c      CALL BLIST(IW,'T+','TGFTTGTLTGCL')
cC
c      KT1FT=IW(NT1FT)
c      KT1TL=IW(NT1TL)
c      KT1CL=IW(NT1CL)

C
C     Calculating bit pattern for each track in a workbank
C
c      IW(1)=1
c      CALL WBANK(IW,NBPWRK,NTRI,*848)

c      write(*,*) "NTRI = " , NTRI

      DO 1 I=1,NTRI
C GAB:
c        IW(NBPWRK+I)=0
          MY_NBMWRK(I)=0
C GAB:
c        KP1=ITABL(KT1TL,I,JT1TIO)
c        NPT=ITABL(KT1TL,I,JT1TN1)
        KP1=MY_JT1TIO(I)
        NPT=MY_JT1TN1(I)

C        write(*,*) "NPT=", NPT
        DO 11 J=1,NPT

C GAB:
C          ICO=IW(KT1CL+LMHLEN+KP1+J)
C          NPAD=ITABL(KTPCO,ICO,JTPCIN)/100000+10
          ICO=MY_T1CLT(KP1+J)
          NPAD=1+INT((TPCO_RV(ICO)-TPCACRI)/TPCPADR)

C GAB:
C           IF (NPAD.GT.10.AND.NPAD.LT.33) CALL SBIT1(IW(NBPWRK+I),NPAD)

          DO INDEX=1,N32BITREG
             IF(NPAD.LE.INDEX*32) THEN
                CALL SBIT1(NBPWRK(I,INDEX),NPAD-((INDEX-1)*32))
             ENDIF
          END DO

C           IF(NPAD.LE.32)THEN
C            CALL SBIT1(NBPWRK(I,1),NPAD)
C           ELSE IF (NPAD.LE.64)THEN
C            CALL SBIT1(NBPWRK(I,2),NPAD-32)
C           ELSE IF (NPAD.LE.96)THEN
C            CALL SBIT1(NBPWRK(I,3),NPAD-64)
C           ELSE IF (NPAD.LE.128)THEN
C            CALL SBIT1(NBPWRK(I,4),NPAD-96)
C           ELSE IF (NPAD.LE.160)THEN
C            CALL SBIT1(NBPWRK(I,5),NPAD-128)
C           ELSE IF (NPAD.LE.192)THEN
C            CALL SBIT1(NBPWRK(I,6),NPAD-160)
C           ELSE
C            CALL SBIT1(NBPWRK(I,7),NPAD-192)
C           ENDIF

   11   CONTINUE

c gab: tmp
c SJA  can't find this used anywhere
ccc        IGNT(I)=IGNTRK(ICO)

    1 CONTINUE
C
C      write(*,*) "finished loops"
C
C     Loop over pairs of tracks to find candidates to link
C
      NG=0
      DO 91 I1=1,NTRI-1
c GAB:
c        NPT1=ITABL(KT1TL,I1,JT1TN1)
        NPT1=MY_JT1TN1(I1)

        IF (NPT1.EQ.0) GOTO 91
        IFF=0
        DO 89 I2=I1+1,NTRI
C GAB:
c          NPT2=ITABL(KT1TL,I2,JT1TN1)
          NPT2=MY_JT1TN1(I2)
          IF (NPT2.EQ.0) THEN
           GOTO 89
          ENDIF

C GAB:
C          IX=IAND(IW(NBPWRK+I1),IW(NBPWRK+I2))
          IX=0
          DO I=1,7
           IX=IX+IAND(NBPWRK(I1,I),NBPWRK(I2,I))
          ENDDO
          IF (IX.NE.0) THEN
           GOTO 89
          ENDIF

          IF (IFF.EQ.0) THEN
            IFF=1
C GAB:
c            IT1F1=KROW(KT1FT,I1)
            IT1F1=I1
C GAB:
C            DO 22 I=JT1FIR,JT1FZ0
c              VV1(I)=RW(IT1F1+I)
c   22       CONTINUE
              VV1(1)=MY_JT1FIR(I1)
              VV1(2)=MY_JT1FTL(I1)
              VV1(3)=MY_JT1FP0(I1)
              VV1(4)=MY_JT1FD0(I1)
              VV1(5)=MY_JT1FZ0(I1)
C GAB:
c            VV1( 6)=SQRT(ABS(RW(IT1F1+JT1FEM)))
c            VV1( 7)=SQRT(ABS(RW(IT1F1+JT1FEM+ 2)))
c            VV1( 8)=SQRT(ABS(RW(IT1F1+JT1FEM+ 5)))
c            VV1( 9)=SQRT(ABS(RW(IT1F1+JT1FEM+ 9)))
c            VV1(10)=SQRT(ABS(RW(IT1F1+JT1FEM+14)))
            VV1( 6)=SQRT(ABS(MY_JT1FEM(I1,1)))
            VV1( 7)=SQRT(ABS(MY_JT1FEM(I1,3)))
            VV1( 8)=SQRT(ABS(MY_JT1FEM(I1,6)))
            VV1( 9)=SQRT(ABS(MY_JT1FEM(I1,10)))
            VV1(10)=SQRT(ABS(MY_JT1FEM(I1,15)))
          ENDIF

C GAB:
c          IT1F2=KROW(KT1FT,I2)
          IT1F2=I2
c          DO 23 I=JT1FIR,JT1FZ0
c            VV2(I)=RW(IT1F2+I)
c  23      CONTINUE
c          VV2( 6)=SQRT(ABS(RW(IT1F2+JT1FEM)))
c          VV2( 7)=SQRT(ABS(RW(IT1F2+JT1FEM+ 2)))
c          VV2( 8)=SQRT(ABS(RW(IT1F2+JT1FEM+ 5)))
c          VV2( 9)=SQRT(ABS(RW(IT1F2+JT1FEM+ 9)))
c          VV2(10)=SQRT(ABS(RW(IT1F2+JT1FEM+14)))
            VV2(1)=MY_JT1FIR(I2)
            VV2(2)=MY_JT1FTL(I2)
            VV2(3)=MY_JT1FP0(I2)
            VV2(4)=MY_JT1FD0(I2)
            VV2(5)=MY_JT1FZ0(I2)
            VV2( 6)=SQRT(ABS(MY_JT1FEM(I2,1)))
            VV2( 7)=SQRT(ABS(MY_JT1FEM(I2,3)))
            VV2( 8)=SQRT(ABS(MY_JT1FEM(I2,6)))
            VV2( 9)=SQRT(ABS(MY_JT1FEM(I2,10)))
            VV2(10)=SQRT(ABS(MY_JT1FEM(I2,15)))
C
C         Check the consistency of the helix parameters of the 2 tracks
C

          CHIP=0.
          DPHI=VV1(3)-VV2(3)
C GAB:
C          IF (ABS(DPHI) .GT. ACOS(-1.)) DPHI=2*ACOS(-1.)-ABS(DPHI)
          IF (ABS(DPHI) .GT. PI) DPHI=TWOPI-ABS(DPHI)

          DO 67 IK=1,5
            IF (IK.EQ.3) THEN    ! PHI
              CHIP=CHIP+DPHI**2/(VV1(IK+5)**2+VV2(IK+5)**2)
            ELSE
              CHIP=CHIP+(VV1(IK)-VV2(IK))**2/(VV1(IK+5)**2+VV2(IK+5)**2)
            ENDIF
   67     CONTINUE
          IF (NPT1.EQ.3.OR.NPT2.EQ.3) CHIP=CHIPM-1.

          IF (CHIP.GT.CHIPM) THEN
           GOTO 89
          ENDIF
C
C         Assemble the TPC coordinate list in increasing pad row number
C
C GAB:
c          DO J1=1,21
          DO J1=1,NRTPC
            LIST1(J1)=0
          ENDDO
C GAB:
c          NPTPC1 =ITABL(KT1TL,I1,JT1TIO)
C          DO J1=1,NPT1
C            ICO=IW(KT1CL+LMHLEN+NPTPC1+J1)
C            NPAD=ITABL(KTPCO,ICO,JTPCIN)/100000
C            IF(NPAD.GT.0.AND.NPAD.LT.22) LIST1(NPAD)=ICO
C          ENDDO
          NPTPC1 =MY_JT1TIO(I1)
          DO J1=1,NPT1
            ICO=MY_T1CLT(NPTPC1+J1)
            NPAD=1+INT((TPCO_RV(ICO)-TPCACRI)/TPCPADR)
            IF(NPAD.GT.0.AND.NPAD.LE.NRTPC) LIST1(NPAD)=ICO
          ENDDO

C GAB:
c          NPTPC2 =ITABL(KT1TL,I2,JT1TIO)
c          DO J2=1,NPT2
c            ICO=IW(KT1CL+LMHLEN+NPTPC2+J2)
c            NPAD=ITABL(KTPCO,ICO,JTPCIN)/100000
c            IF(NPAD.GT.0.AND.NPAD.LT.22) LIST1(NPAD)=ICO
c          ENDDO
          NPTPC2 =MY_JT1TIO(I2)
          DO J2=1,NPT2
            ICO=MY_T1CLT(NPTPC2+J2)
            NPAD=1+INT((TPCO_RV(ICO)-TPCACRI)/TPCPADR)
            IF(NPAD.GT.0.AND.NPAD.LE.NRTPC) LIST1(NPAD)=ICO
          ENDDO


          NTPC=0
C GAB:
c          DO IJ=1,21
          DO IJ=1,NRTPC
            IF(LIST1(IJ).NE.0) THEN
              NTPC=NTPC+1
              LIST1(NTPC)=LIST1(IJ)
            ENDIF
          ENDDO
C GAB:
C          DO IJ=NTPC+1,21
          DO IJ=NTPC+1,NRTPC
            LIST1(IJ)=0
          ENDDO



C
C         Fit the track candidate. Remove upto LOMAX points to get a
C         good fit, as in TFITTK.
C
          LOST=0
          IOPT=2
c gab: multiple scattering terms (1=on, 0=off)
          NMULT=1
c          NMULT=0

c GAB: I don't know why this was put in... so I am removing it...
c          NTPC=NPT1+NPT2
          NITC=0
          NVDET=0
  333     CONTINUE

          CALL UFITMS(NTPC,NITC,NVDET,LIST1,LIST2,LIST3,IOPT,NMULT
     +                              ,RMS,RADLEN,VV0,EE0,COV,CHI2,IERR)


C SJA:Commented out the hard coded CHI2 condition inserted below, and reinstated the original value of PPP
          IF ( CHI2.GT.0 .AND. (2*NPT-5).GT.0 .AND. NPT.GT.0 ) THEN
             PPP = PROB(CHI2,(2*NPT-5))
          ELSE
             WRITE(*,*) 
     +            "tflnk2: CHI2 or NDF sent to CERNLIB PROB() are zero"
             PPP = 1e-25
          ENDIF
          
c MVL          PPP=PROB(CHI2,2*NTPC-5)
C      IF(CHI2.GT.0.0.and.CHI2.lT.300.0.and.2*NTPC-5.GT.0)THEN
C        PPP=max(PROB(CHI2,(2*NTPC-5)),1.e-25)
C      ELSE
C        PPP=1.e-25
C      ENDIF
c MVL
          IF (PPP.LT.PMIN) THEN
            IF (NTPC.GT.MINPT.AND.LOST.LT.LOMAX) THEN
              LOST=LOST+1
              CMAX=0.
              IMAX=0
              DO 234 IKK=1,NTPC
                IF(CHIACU(IKK).GT.CMAX) THEN
                  CMAX=CHIACU(IKK)
                  IMAX=IKK
                ENDIF
  234         CONTINUE
              J2=0
              DO J1=1,NTPC
                IF (J1.NE.IMAX) THEN
                  J2=J2+1
                  LIST1(J2)=LIST1(J1)
                ENDIF
              ENDDO
              LIST1(NTPC)=0
              NTPC=NTPC-1
              GOTO 333
            ELSE
              GOTO 89
            ENDIF
          ENDIF
C  Check for garbage fit
          CALL TFCHEK(VV0,.FALSE.,IERR)

          IF (IERR.EQ.0) THEN
C
C  Pair found, store relevant parameters
C
           IF(NG.LT.MXGD) THEN
             NG=NG+1
             LGOOD(NG,1)=I1
             LGOOD(NG,2)=I2
             LGOOD(NG,3)=LOST
             LGOOD(NG,4)=NINT(1000000.0*CHI2)
             LGOOD(NG,5)=0
             LGOOD(NG,6)=NTPC

C GAB
C             DO I=1,21
             DO I=1,NRTPC
              LGOOD(NG,6+I)=LIST1(I)
             ENDDO
             DO I=1,5
              QGOOD(NG,I)=VV0(I)
             ENDDO
             DO I=1,15
              QGOOD(NG,5+I)=COV(I)
             ENDDO
             QGOOD(NG,21)=CHI2
             QGOOD(NG,22)=PPP
            ELSE
             WRITE(*,*)'TFLNK2: Overflow of buffer for track pairs'
             WRITE(*,*)'Some TPC track links lost'
             CALL TKVETO(1)
            ENDIF
           ENDIF
 89      CONTINUE
 91    CONTINUE

C------------------------------------------------------------
C Reset the TPC track counters:
      NT=0
      NCL=0

C
C     If 2 links have a track in common, flag the link with the
C     largest chi**2.
C
      IF (NG.GT.1) THEN
        DO 107 IT1=1,NG-1
C GAB:
C          IF (LGOOD(IT1,5).EQ.1) GOTO 107
          IF ((LGOOD(IT1,5).EQ.1).OR.(LGOOD(IT1,5).LT.0))GOTO 107
          N1A=LGOOD(IT1,1)
          N1B=LGOOD(IT1,2)
          ICH1=LGOOD(IT1,4)
          DO 7 IT2=IT1+1,NG
            IF (LGOOD(IT2,5).EQ.1) GOTO 7
            N2A=LGOOD(IT2,1)
            N2B=LGOOD(IT2,2)
            IF (N1A.NE.N2A .AND. N1A.NE.N2B .AND. N1B.NE.N2A
     +                                   .AND. N1B.NE.N2B) GOTO 7
            ICH2=LGOOD(IT2,4)
            IF (ICH1.LT.ICH2) THEN
              LGOOD(IT2,5)=1
            ELSE
              LGOOD(IT1,5)=1
              ICH1=ICH2
            ENDIF
    7     CONTINUE
  107   CONTINUE
      ENDIF
      NPREG=0
      NTASS=0
      DO 123 IT=1,NG
        IF (LGOOD(IT,5).EQ.0) THEN
          NTASS=NTASS+1
          NPREG=NPREG+LGOOD(IT,3)
        ENDIF
  123 CONTINUE

C
       NTGFT = NTRI - NTASS
c GAB:
c      NTGCL = LROWS(KT1CL) - NPREG
      NTGCL = NT1CL - NPREG

C(KH)----begin
      IF (NTGCL.GT.MAXCND) THEN
        WRITE(6,*) 'TFLNK2: MAXCND TOO SMALL. NEED ',NTGCL
        CALL TKVETO(1)
        NTGCL=MAXCND
      END IF
C(KH)----end

c      IW(KTGFT+LMHROW)=NTGFT
c      IW(KTGTL+LMHROW)=NTGFT
c      IW(KTGCL+LMHROW)=NTGCL

C
C++   Loop over all tracks
C
      



      DO 567 I0=1,NTRI
C
C++     Look for a link
C
        DO 577 IT=1,NG
          I1=LGOOD(IT,1)
          I2=LGOOD(IT,2)
          IF (I1.NE.I0.AND.I2.NE.I0) GOTO 577
C
C++       Here is a link with another track
C
          NFL=LGOOD(IT,5)
C GAB:
          IF(NFL.LT.0)GOTO 577

C
C++       Skip this link if a better (smaller chi**2) link exists
C
          IF (NFL.EQ.1) GOTO 577
C
C++       Skip this track if it already was linked to a previous one
C
          IF (NFL.EQ.2) GOTO 567
C
C++       Copy fit results to the new bank
C
          LGOOD(IT,5)=2
C
          NTPC=LGOOD(IT,6)
C GAB
C          DO I=1,21
          DO I=1,NRTPC
            LIST1(I)=LGOOD(IT,6+I)
          ENDDO
          DO I=1,5
            VV0(I)=QGOOD(IT,I)
          ENDDO
          DO I=1,15
            COV(I)=QGOOD(IT,5+I)
          ENDDO
          CHI2=QGOOD(IT,21)
C
          NT=NT+1

C GAB:
c          IW(KROW(KT1FT,I1)+JT1FTG) = NT
c          IW(KROW(KT1FT,I2)+JT1FTG) = NT
c          ITGFT=KROW(KTGFT,NT)
          MY_JT1FTG(I1) = NT
          MY_JT1FTG(I2) = NT
          ITGFT=NT
C GAB:
c          CALL UCOPY(VV0,RW(ITGFT+JTGFIR),5)
c          CALL UCOPY(COV,RW(ITGFT+JTGFEM),15)
c          RW(ITGFT+JTGFCD)=CHI2
c          IW(ITGFT+JTGFDF)=2*NTPC-5
c          IW(ITGFT+JTGFTC)=0  ! Pointer to fitted FRFT track

         MY_JTGFIR(NT)=VV0(1)
         MY_JTGFTL(NT)=VV0(2)
         MY_JTGFP0(NT)=VV0(3)
         MY_JTGFD0(NT)=VV0(4)
         MY_JTGFZ0(NT)=VV0(5)
         DO K=1,15
           MY_JTGFEM(NT,K)=COV(K)
         ENDDO
         MY_JTGFCD(NT)=CHI2
         MY_JTGFDF(NT)=2*NTPC-5
         MY_JTGFTC(NT)=0  ! Pointer to fitted FRFT track

C GAB:
c          ITGTL=KROW(KTGTL,NT)
c          IW(ITGTL+JTGTIO)=NCL
c          IW(ITGTL+JTGTN1)=NTPC
          ITGTL=NT
          MY_JTGTIO(ITGTL)=NCL
          MY_JTGTN1(ITGTL)=NTPC

C GAB:
c          NRT1=ITABL(KT1TL,I1,JTGTNR)
c          NRT2=ITABL(KT1TL,I2,JTGTNR)
c          NRT=NRT1+NRT2
c          IW(ITGTL+JTGTNR)=NRT
          NRT1=MY_JT1TNR(I1)
          NRT2=MY_JT1TNR(I2)
          NRT=NRT1+NRT2
          MY_JTGTNR(ITGTL)=NRT
C GAB:
c          CALL UCOPY(LIST1,IW(KTGCL+LMHLEN+NCL+1),NTPC)
          CALL UCOPY(LIST1,MY_TGCLT(NCL+1),NTPC)

          IF (NRT1.NE.0) THEN
C GAB:
c            KP1=ITABL(KT1TL,I1,JTGTIO)
c            NPT1=ITABL(KT1TL,I1,JTGTN1)
c            CALL UCOPY(IW(KT1CL+LMHLEN+KP1+NPT1+1),
c     +                       IW(KTGCL+LMHLEN+NCL+NTPC+1),NRT1)
            KP1=MY_JT1TIO(I1)
            NPT1=MY_JT1TN1(I1)
            CALL UCOPY(MY_T1CLT(KP1+NPT1+1),
     +                       MY_TGCLT(NCL+NTPC+1),NRT1)
          ENDIF
          IF (NRT2.NE.0) THEN
C GAB:
c            KP2=ITABL(KT1TL,I2,JTGTIO)
c            NPT2=ITABL(KT1TL,I2,JTGTN1)
c            CALL UCOPY(IW(KT1CL+LMHLEN+KP2+NPT2+1),
c     +                       IW(KTGCL+LMHLEN+NCL+NTPC+NRT1+1),NRT2)
            KP2=MY_JT1TIO(I2)
            NPT2=MY_JT1TN1(I2)
            CALL UCOPY(MY_T1CLT(KP2+NPT2+1),
     +                       MY_TGCLT(NCL+NTPC+NRT1+1),NRT2)
          ENDIF
          NCL=NCL+NTPC+NRT
          GOTO 567
  577   CONTINUE
C
C++     Copy the unchanged tracks
C

        NT=NT+1
        
C     GAB:
c     IT1FT=KROW(KT1FT,I0)
c     IW(IT1FT+JT1FTG) = NT
c     ITGFT=KROW(KTGFT,NT)
c     CALL UCOPY(RW(IT1FT+JTGFIR),RW(ITGFT+JTGFIR),5)
c     CALL UCOPY(RW(IT1FT+JTGFEM),RW(ITGFT+JTGFEM),15)
c     RW(ITGFT+JTGFCD)=RW(IT1FT+JTGFCD)
c     IW(ITGFT+JTGFDF)=IW(IT1FT+JTGFDF)
c     IW(ITGFT+JTGFTC)=IW(IT1FT+JT1FTC)
        IT1FT=I0
        MY_JT1FTG(IT1FT) = NT
        ITGFT=NT
        MY_JTGFIR(ITGFT)=MY_JT1FIR(IT1FT)
        MY_JTGFTL(ITGFT)=MY_JT1FTL(IT1FT)
        MY_JTGFP0(ITGFT)=MY_JT1FP0(IT1FT)
        MY_JTGFD0(ITGFT)=MY_JT1FD0(IT1FT)
        MY_JTGFZ0(ITGFT)=MY_JT1FZ0(IT1FT)
        DO  K=1,15
           MY_JTGFEM(ITGFT,K)=MY_JT1FEM(IT1FT,K)
        ENDDO
        MY_JTGFCD(ITGFT)=MY_JT1FCD(IT1FT)
        MY_JTGFDF(ITGFT)=MY_JT1FDF(IT1FT)
        MY_JTGFTC(ITGFT)=MY_JT1FTC(IT1FT)
        
C GAB:
c     IT1TL=KROW(KT1TL,I0)
c     ITGTL=KROW(KTGTL,NT)
c     IW(ITGTL+JTGTIO)=NCL
c     IW(ITGTL+JTGTN1)=IW(IT1TL+JTGTN1)
c     IW(ITGTL+JTGTNR)=IW(IT1TL+JTGTNR)
c     KP=ITABL(KT1TL,I0,JTGTIO)
c     NTOT=IW(ITGTL+JTGTN1)+IW(ITGTL+JTGTNR)
c     CALL UCOPY(IW(KT1CL+LMHLEN+KP+1),IW(KTGCL+LMHLEN+NCL+1),NTOT)
c     NCL=NCL+NTOT
        IT1TL=I0
        ITGTL=NT
        MY_JTGTIO(ITGTL)=NCL
        MY_JTGTN1(ITGTL)=MY_JT1TN1(IT1TL)
        MY_JTGTNR(ITGTL)=MY_JT1TNR(IT1TL)
        KP=MY_JT1TIO(I0)
        NTOT=MY_JTGTN1(ITGTL)+MY_JTGTNR(ITGTL)
        CALL UCOPY(MY_T1CLT(KP+1),MY_TGCLT(NCL+1),NTOT)
        NCL=NCL+NTOT
C     

 567  CONTINUE
c      write(*,*) "got this far"
C
C Update the pointers in TPCO to the TGFT bank
C
C GAB:
c      DO 705 IC=1,LROWS(KTPCO)
c        IW(KROW(KTPCO,IC)+JTPCTN)=0
c  705 CONTINUE
c      DO 756 ITK=1,LROWS(KTGTL)
c        IOFF=ITABL(KTGTL,ITK,JTGTIO)+LMHLEN
c        N1=ITABL(KTGTL,ITK,JTGTN1)
c        NTOT=N1+ITABL(KTGTL,ITK,JTGTNR)
c        DO 732 II=1,N1
c          IC=IW(KTGCL+IOFF+II)
c          IW(KROW(KTPCO,IC)+JTPCTN)=ITK
c  732   CONTINUE
c        DO 742 II=N1+1,NTOT
c          IC=IW(KTGCL+IOFF+II)
c          IW(KROW(KTPCO,IC)+JTPCTN)=-ITK
c  742   CONTINUE
c  756 CONTINUE



      DO 705 IC=1,NTPCO
        TPCO_TN(IC)=0
  705 CONTINUE
      DO 756 ITK=1,NTGCL
        IOFF=MY_JTGTIO(ITK)
        N1=MY_JTGTN1(ITK)
        NTOT=N1+MY_JTGTNR(ITK)
        DO 732 II=1,N1
          IC= MY_TGCLT(IOFF+II)
          TPCO_TN(IC)=ITK
  732   CONTINUE

        DO 742 II=N1+1,NTOT
          IC=MY_TGCLT(IOFF+II)
c(kh)     TPCO_TN(IC)=-ITK
c(kh)--bug workaround
          IF (IC.GT.0) THEN
             TPCO_TN(IC)=-ITK
          ELSE
c             WRITE(6,*) 'TFLNK2: internal error. IC=',IC
             CALL TKVETO(1)
          END IF
c(kh)--end
  742   CONTINUE
  756 CONTINUE

C GAB:
      NTGFT=ITGFT


cC
cC Compress the new banks
cC
c      CALL AUBPRS('TGFTTGTLTGCL')
cC
c      IER=0
c      CALL WDROP(IW,NBPWRK)
c      RETURN
c  848 CONTINUE
c      IER=11
c      CALL RERROR('TFLNK2',IER,'No room to create a workbank.')
c      KTGFT=NDROP('TGFT',1)
c      KTGTL=NDROP('TGTL',1)
c      KTGCL=NDROP('TGCL',1)

      END
C-----------------------------------------------------------------
