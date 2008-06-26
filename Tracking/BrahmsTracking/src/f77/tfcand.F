*  
* $Id: tfcand.F,v 1.4 2008-06-26 09:35:56 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.3  2008/01/23 17:22:22  aplin
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
      SUBROUTINE TFCAND(IER)
C****************************************************************************
c* Version number: 0.0      Date of last modification  26.08.99             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
C* Main modifiers from ALEPH to BRAHMS code:                                *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
C-----*****************---------------------------------------------
      IMPLICIT NONE
C! Assemble tpc chains to track candidates
C!
C!    AUTHOR        :- MARTIN POPPE 87/5/25
C!    Modified      :- Robert Johnson 88/9/14
C!
C!    OUTPUT: IER  = ERROR CODE: =0 IF OK, =1 IF INCOMPLETE
C!
C!    INPUT  BANKS  : 'TCHA' , 'TCTC' , 'TLNK'
C!    OUTPUT BANKS  : 'TCAL' , 'TTCC' , 'TARC', 'TATC'
C!
C!    DESCRIPTION
C!    ===========
C!
C!    According to the table of links between chains provided by
C!    the bank 'tlnk', sets of chains are assembled to form the
C!    so called track candidates. Track candidates are sets of
C!    tpc pad coordinates which are assumed to originate in the
C!    same particle. the coordinates are ordered according to
C!    their path of flight. For non spiraling tracks, this path
C!    is assumed to point away from the beam line. For spiraling
C!    tracks, the order is deduced from the assumption that
C!    the particles loose energy as they traverse the chamber.
C!    All pad coordinates that belong to a chain will also
C!    belong to a track candidate. Every chain will be part
C!    of exactly one track candidate.
C!    Also, the TARC and TATC banks are filled with the chains
C!    which are not in the first arc of a track candidate.
C!    This is done in order to make available the helix parameters
C!    of the extra loops of track candidates.
C!
C-------------------------------------------------------------------
C
#include "include/fndparms.inc"
#include "include/tarcjj.inc"
#include "include/tlnkjj.inc"
#include "include/tpgeom.inc"
c#include "marlin_tpcgeom.h"
c---------------------------------------------------------------------
c my stuff:
#include "include/tpco.inc"
#include "include/my_tcha.inc"
#include "include/my_tlnk.inc"
#include "include/my_tctc.inc"
#include "include/my_tatc.inc"
#include "include/my_tcal.inc"
#include "include/my_tarc.inc"
#include "include/my_ttcc.inc"
c#include "gcflag.inc"
#include "include/tcanwb.inc"

      INTEGER MY_ITCALW(NCHMAX),MY_JPLNKW(NCHMAX),MY_JTYPEW(NCHMAX)
      INTEGER IBMAP(N32BITREG),INEW
c------------------------------------------------------------------------
      INTEGER MXFIT
      PARAMETER (MXFIT=LTPDRO+2)
      INTEGER IBUFF
      REAL XFT,YFT,EFT
      DIMENSION IBUFF(LTPDRO),XFT(MXFIT),YFT(MXFIT),EFT(MXFIT)
cSJA      DIMENSION IBUFF(500),XFT(MXFIT),YFT(MXFIT),EFT(MXFIT)
cSJA IBUFF should be the same dimention as the max number of tpcrows      
      LOGICAL FIRST,TAKEN,FSWAP
      INTEGER CDUMY,CTYPE

      INTEGER I,J,II,IER,NCHAN,LTATC,NPUT,NSTOR,ICHAN,NLNKS
      INTEGER KC1,KC2,KC3,KC4,JCH,NJCH,NTB,KCH,NFIT,LAST,IDUMY
      INTEGER KKTCTC,IPLNK,IORDR,ISTRT,IPICE,JCHA1,NARC1,NPT,KSCHC,IPT
      INTEGER LPT,KTAKE,JCHAN,NARCS,LMOFF,LM,N1ARC

      REAL PT,SPT,RNPT,SX,SX2,SY,SXY,SS,B,SB,W,DET
      REAL D01,Z01,D02,Z02

      REAL SCUT,D0CUT,Z0CUT
      SAVE SCUT,D0CUT,Z0CUT

      DATA FIRST/.TRUE./
C
C---> The following parameters should eventually go on data base:
C
      DATA SCUT/60.0/,D0CUT/100./,Z0CUT/2./
C
C     SCUT= significance cut on the rate of change of curvature of
C           a track.  If the measured rate of change divided by
C           the estimated uncertainty is less in absolute value than
C           SCUT, then this is assumed not to be a significant
C           measure of which direction the particle is going.
C           If the change of curvature is not significant, then
C           take the end of the helix closest to the origin in
C           z if the two ends are further apart in z than Z0CUT
C           and at least one end comes within D0CUT of the beamline.
C           If those conditions aren't satified, use the curvature to
C           decide, even though it is not significant.
C
C
C     Description of workbank ITCALW:
C
C     Words 1..nchains        Flag for used chains
C     nchains+1..2*nchains    list of chains belonging to candidate
C     2*nchains+1..3*nchains  link type to next chain
C
C---------------------- Statement Functions -------------------------
C
C#include "bmacro.h"
C
C GAB:
cC---> Offset in TCTC for first coordinate of track II
C GAB:
c      KKTCTC(II) = KTCTC+LMHLEN+ITABL(KTCHA,II,JTCHOF)
       KKTCTC(II) = MY_JTCHOF(II)
cC---> Functions to reference the workbank ITCALW
c      TAKEN(II) = IW(ITCALW+II).EQ.1
C GAB:
       TAKEN(II) = MY_ITCALW(II).EQ.1

cC---> IPLNK = Array of chains belonging to the track candidate
C GAB:
c      IPLNK(II) = IW(JPLNKW+II)
      IPLNK(II) = MY_JPLNKW(II)
cC---> CTYPE = Link types to next chain for each chain of candidate
C GAB:
c      CTYPE(II) = IW(JTYPEW+II)
      CTYPE(II) = MY_JTYPEW(II)
C
C--------------------------------------------------------------------
C
      IF (FIRST) THEN
        ITCALW=0
c        NTCHA = NAMIND('TCHA')
c        NTCTC = NAMIND('TCTC')
c        NTLNK = NAMIND('TLNK')
c        NTCAL = NAMIND('TCAL')
c        NTTCC = NAMIND('TTCC')
c        NTARC = NAMIND('TARC')
c        NTATC = NAMIND('TATC')
        FIRST = .FALSE.
      END IF
      IER=0

c GAB:
c      KTCHA = IW(NTCHA)
c      IF (KTCHA.LE.0)                               GO TO 998
      IF (MY_NCHAIN.LE.0)                               GO TO 998
C      KTCTC = IW(NTCTC)
C      IF (KTCTC.LE.0)                               GO TO 998
C      KTLNK = IW(NTLNK)
C      IF (KTLNK.LE.0)                               GO TO 998
C GAB:
C      NCHAN = LROWS(KTCHA)
      NCHAN = MY_NCHAIN
C
C-----> create banks of max possible size / drop old ones if necessary
C
c GAB:
c      IF (IW(NTTCC).NE.0) CALL BDROP(IW,'TTCC')
c      IF (IW(NTCAL).NE.0) CALL BDROP(IW,'TCAL')
c      IF (IW(NTARC).NE.0) CALL BDROP(IW,'TARC')
c      IF (IW(NTATC).NE.0) CALL BDROP(IW,'TATC')
C
c      CALL BLIST(IW,'T+','TARCTTCCTCALTATC')
c      NCOOD = LROWS(KTCTC)
c      LTTCC = LMHLEN+NCOOD
c      IW(1)=1
c      CALL AUBOS('TTCC',1,LTTCC,KTTCC,IGARB)
c      IF (IGARB.EQ.2) GO TO 999
c      IW(KTTCC+LMHCOL) = 1
c      LTCAL = LMHLEN + LTCALA*NCHAN
c      IW(1)=1
c      CALL AUBOS('TCAL',1,LTCAL,KTCAL,IGARB)
c      IF (IGARB.EQ.2) GO TO 999
c      IW(KTCAL+LMHCOL) = LTCALA
C
c      LTARC= NCHAN*LTARCA + LMHLEN
c      IW(1)=1
c      CALL AUBOS('TARC',0,LTARC,KTARC,IGARB)
c      IF (IGARB.EQ.2) GO TO 999
c      IW(KTARC+LMHCOL)=LTARCA
      NARCS=0                      ! KH uncommented this line
c      LTATC= IW(KTCTC)
c      IW(1)=1
c      CALL AUBOS('TATC',0,LTATC,KTATC,IGARB)
c      IF (IGARB.EQ.2) GO TO 999
c      IW(KTATC+LMHCOL)=1

      LTATC=0
C
c      KTCHA = IW(NTCHA)
c      KTCTC = IW(NTCTC)
c      KTLNK = IW(NTLNK)
c      KTCAL = IW(NTCAL)
c      KTARC = IW(NTARC)
c      KTTCC = IW(NTTCC)
cC
cC---> Create workbank for temporary storage of candidates
cC
c      IW(1)=1
c      CALL WBANK(IW,ITCALW,3*NCHAN,*995)
c      IW(ITCALW-3)=INTCHA('TCAW')
c      CALL VZERO(IW(ITCALW+1),NCHAN)
c      JPLNKW = ITCALW+NCHAN
c      JTYPEW = ITCALW+2*NCHAN
C GAB:
       CALL VZERO(MY_ITCALW(1),NCHMAX)
C
c GAB:
c      NPUT  = KTTCC+LMHLEN
      NPUT=0
C          (= absolute offset for next coordinates to be stored)
      NSTOR = 0
C          (= number of coordinates stored already )
      NCAND = 0
C          (= number of track candidates)

C-------------------------------------------------------------------C
C              TREATMENT OF CHAINS WITH A SINGLE LINK               C
C              **************************************               C
C-------------------------------------------------------------------C

      DO 400 ICHAN = 1,NCHAN
       IF (TAKEN(ICHAN))                       GO TO 400

c GAB:
C        KSLNK = KROW(KTLNK,ICHAN)
        NLNKS  = 0
C GAB:
C        KC1    = IW(KSLNK+JTLNUP)
C        KC2    = IW(KSLNK+JTLNDN)
C        KC3    = IW(KSLNK+JTLNTP)
C        KC4    = IW(KSLNK+JTLNBM)

        KC1    = MY_JTLNUP(ICHAN)
        KC2    = MY_JTLNDN(ICHAN)
        KC3    = MY_JTLNTP(ICHAN)
        KC4    = MY_JTLNBM(ICHAN)

        IF (KC1.GT.0 .AND. (.NOT.TAKEN(KC1))) NLNKS = NLNKS +1
        IF (KC2.GT.0 .AND. (.NOT.TAKEN(KC2))) NLNKS = NLNKS +1
        IF (KC3.GT.0 .AND. (.NOT.TAKEN(KC3))) NLNKS = NLNKS +1
        IF (KC4.GT.0 .AND. (.NOT.TAKEN(KC4))) NLNKS = NLNKS +1

        IF (NLNKS.NE.1)                          GO TO 400
        NCAND = NCAND + 1

C GAB:
C        KSCAL = KROW(KTCAL,NCAND)
C        IW(KSCAL+JTCAIO) = NSTOR
        MY_JTCAIO(NCAND)=NSTOR

        JCH = ICHAN
        NJCH= 1
        NTB = 0
C GAB:
C        IW(JPLNKW+1) = ICHAN
C        IW(JPLNKW+NJCH) = JCH
        MY_JPLNKW(1) = ICHAN
        MY_JPLNKW(NJCH) = JCH
C
C-----> Mark the chain used.  Then loop over linked chains.
C
C GAB:
C        IW(ITCALW+JCH) = 1
        MY_ITCALW(JCH) = 1

  110   CONTINUE
C GAB:
c          KSLNK = KROW(KTLNK,JCH)
          INEW=JCH
          KCH = 0
C GAB:
c          JCH = IW(KSLNK+JTLNUP)
          JCH = MY_JTLNUP(INEW)
C


          IF (JCH.NE.0 .AND. (.NOT.TAKEN(JCH))) THEN
C GAB:
c            IW(JTYPEW+NJCH)=JTLNUP
             MY_JTYPEW(NJCH)=JTLNUP
            KCH = JCH
          END IF
C GAB:
c          JCH = IW(KSLNK+JTLNDN)
          JCH = MY_JTLNDN(INEW)
          IF (JCH.NE.0 .AND. (.NOT.TAKEN(JCH))) THEN
C GAB:
c            IW(JTYPEW+NJCH)=JTLNDN
            MY_JTYPEW(NJCH)=JTLNDN
            KCH = JCH
          END IF
C GAB
c          JCH = IW(KSLNK+JTLNTP)
          JCH = MY_JTLNTP(INEW)
          IF (JCH.NE.0 .AND. (.NOT.TAKEN(JCH))) THEN
C GAB:
c            IW(JTYPEW+NJCH)=JTLNTP
            MY_JTYPEW(NJCH)=JTLNTP
            KCH = JCH
            NTB = NTB+1
          END IF
C GAB
C          JCH = IW(KSLNK+JTLNBM)
          JCH = MY_JTLNBM(INEW)
          IF (JCH.NE.0 .AND. (.NOT.TAKEN(JCH))) THEN
C GAB:
c            IW(JTYPEW+NJCH)=JTLNBM
            MY_JTYPEW(NJCH)=JTLNBM
            KCH = JCH
            NTB = NTB+1
          END IF
          IF (KCH.EQ.0)           GO TO 111

          JCH  = KCH
          NJCH = NJCH + 1
C GAB
c          IW(JPLNKW+NJCH)=JCH
          MY_JPLNKW(NJCH)=JCH
C
C-------> Mark the chain used, then go to the next link
C
C GAB:
c          IW(ITCALW+JCH) = 1
          MY_ITCALW(JCH) = 1
          GO TO 110
  111   CONTINUE

c GAB:
c        IW(JTYPEW+NJCH)=-1
        MY_JTYPEW(NJCH)=-1
C
C--->   NCAND = current track candidate number
C       NTB   = number of tob / bottom links
C       NJCH  = number of chains making the spiral

C
C--->   Now decide in which order chains are stored and reverse the
C       order of the list if necessary.   The algorithm is:
C       Fit a straight line to the points x=link number vs. y=1/R.
C       If the slope is positive, then the list is in the right order,
C       otherwise it must be reversed.  If the slope is consistant with
C       zero and the links pass close to the beam line, then the
C       first arc is that with the smallest z0.
C
        FSWAP=.FALSE.
        IF (NJCH.GT.1 .AND. NTB.GT.0) THEN

          NFIT=0
          KCH=1
  113     IF (KCH.GT.NJCH) GO TO 118
            PT=0.
            SPT=0.
            RNPT=0.
  114       CONTINUE

C GAB:
C              PT=PT+ABS(RTABL(KTCHA,IPLNK(KCH),JTCHIR))
c              SPT=SPT + 1./RTABL(KTCHA,IPLNK(KCH),JTCHER)
              PT=PT+ABS(MY_JTCHIR(IPLNK(KCH)))
C(KH)         SPT=SPT + 1./MY_JTCHER(IPLNK(KCH),1)
C(KH)---dirty workaround to prevent division by zero errors
              IF (MY_JTCHER(IPLNK(KCH),1).NE.0.0) THEN
                SPT=SPT + 1./MY_JTCHER(IPLNK(KCH),1)
              ELSE
                SPT=SPT+1.E10
              END IF
C(KH)---end
              RNPT=RNPT+1.
              IF (CTYPE(KCH).NE.JTLNUP .AND. CTYPE(KCH).NE.JTLNDN) THEN
                GO TO 115
              ENDIF
              KCH=KCH+1
              GO TO 114
  115       CONTINUE

            PT=PT/RNPT
            SPT=SPT/RNPT
            IF (NFIT.LT.MXFIT) THEN
              NFIT=NFIT+1
              XFT(NFIT)=FLOAT(NFIT)
              YFT(NFIT)=PT
              EFT(NFIT)=SPT
            ENDIF
            KCH=KCH+1
            GO TO 113
  118     CONTINUE
C
C--->     Weighted fit to a straight line
C
          SX=0.
          SX2=0.
          SY=0.
          SXY=0.
          SS=0.
          DO 119 I=1,NFIT
            W=1./EFT(NFIT)
            SX=SX+XFT(I)*W
            SY=SY+YFT(I)*W
            SXY=SXY+XFT(I)*YFT(I)*W
            SX2=SX2+(XFT(I)**2)*W
            SS=SS+W
  119     CONTINUE
          DET= SS*SX2 - SX**2
          IF (DET.LE.0.) THEN
            B=1.0
            SB=100.
          ELSE
            B= (SS*SXY - SX*SY)/DET
            SB= SQRT(SS/DET)
          ENDIF
C
C--->     Decide which end of the helix is the start
C
          IF (B/SB.LT.-SCUT) THEN
            FSWAP=.TRUE.
          ELSEIF (B/SB.LT.SCUT) THEN
C GAB:
c            D01=ABS(RTABL(KTCHA,IPLNK(1),JTCHD0))
c            Z01= RTABL(KTCHA,IPLNK(1),JTCHZ0)
c            D02=ABS(RTABL(KTCHA,IPLNK(NJCH),JTCHD0))
c            Z02= RTABL(KTCHA,IPLNK(NJCH),JTCHZ0)

            D01=ABS(MY_JTCHD0(IPLNK(1)))
            Z01= MY_JTCHZ0(IPLNK(1))
            D02=ABS(MY_JTCHD0(IPLNK(NJCH)))
            Z02= MY_JTCHZ0(IPLNK(NJCH))

            IF (ABS(Z02).LT.ABS(Z01) .AND. D02.LT.D0CUT) THEN
              FSWAP=.TRUE.
            ELSEIF (D01.GE.D0CUT) THEN
              IF (B.LT.0.) FSWAP=.TRUE.
            ENDIF
          ENDIF
        ENDIF
C
C--->   Reverse the order of chains if necessary
C

        IF (FSWAP) THEN
          LAST = NJCH/2
          DO 120 I=1,LAST
            IDUMY = IPLNK(I)
            CDUMY = CTYPE(I)
C GAB:
c            IW(JPLNKW+I) = IPLNK(1+NJCH-I)
c            IW(JTYPEW+I) = CTYPE(NJCH-I)
c            IW(JPLNKW+1+NJCH-I) = IDUMY
c            IW(JTYPEW+NJCH-I) = CDUMY

            MY_JPLNKW(I) = IPLNK(1+NJCH-I)
            MY_JTYPEW(I) = CTYPE(NJCH-I)
            MY_JPLNKW(1+NJCH-I) = IDUMY
            MY_JTYPEW(NJCH-I) = CDUMY
  120     CONTINUE
          DO 121 I=1,NJCH
            CDUMY = CTYPE(I)
C GAB:
c            IF (CDUMY.EQ.JTLNUP) IW(JTYPEW+I) = JTLNDN
c            IF (CDUMY.EQ.JTLNDN) IW(JTYPEW+I) = JTLNUP
            IF (CDUMY.EQ.JTLNUP) MY_JTYPEW(I) = JTLNDN
            IF (CDUMY.EQ.JTLNDN) MY_JTYPEW(I) = JTLNUP
  121     CONTINUE
        END IF
C
C--->   Decide order of points in the first chain
C
        IORDR = 1
        IF (CTYPE(1).EQ.JTLNDN .OR. CTYPE(1).EQ.JTLNBM) IORDR= -1
        IF (NTB.EQ.0) IORDR= 1
C
C--->   Reset intermediates: ISTRT = start addresse in buffer
C                            IPICE = 1 -> NJCH
C                            JCHA1 = 1st chain in IBUFF
        ISTRT = 0
        IPICE = 0
        JCHA1 = IPLNK(1)

C GAB
C        IBMAP = 0
        DO I=1,N32BITREG
         IBMAP(I) = 0
        ENDDO
        NARC1 = 1
C
C-----> Loop over chains of the candidate and fill the TTCC list
C       of coordinates in order of direction of flight of the particle
C
  150   CONTINUE
          IPICE = IPICE + 1
          JCH = IPLNK(IPICE)
C GAB:
C          IBMAP = IOR(IBMAP,ITABL(KTCHA,JCH,JTCHPI))
c          NPT = ITABL(KTCHA,JCH,JTCHNC)
          DO I=1,N32BITREG
           IBMAP(I) = IOR(IBMAP(I),MY_JTCHPI(JCH,I))
          ENDDO
          NPT = MY_JTCHNC(JCH)

ccc gab for info: C---> Offset in TCTC for first coordinate of track II
cccc gab for info:      KKTCTC(II) = MY_JTCHOF(II)
          KSCHC = KKTCTC(JCH)

          DO 200 IPT =1,NPT
C GAB:
C            IBUFF(ISTRT+IPT) = IW(KSCHC+IPT)
            IBUFF(ISTRT+IPT) = MY_TCTC(KSCHC+IPT)
  200     CONTINUE

          IF (ISTRT.NE.0) THEN
C
C------->   Order buffer by ripple sort in case of up/down link
C
            LPT = NPT+ISTRT
            DO 240 I=1,LPT-1
              FSWAP=.FALSE.
              DO 230 J=1,LPT-1
                IF (IBUFF(J+1).LT.IBUFF(J)) THEN
                  KTAKE = IBUFF(J+1)
                  IBUFF(J+1) = IBUFF(J)
                  IBUFF(J)   = KTAKE
                  FSWAP=.TRUE.
                END IF
  230         CONTINUE
              IF (.NOT.FSWAP) GO TO 241
  240       CONTINUE
  241       CONTINUE
          END IF
C
C-------> End of buffer reordering
C
          IF (CTYPE(IPICE).EQ.JTLNUP.OR.CTYPE(IPICE).EQ.JTLNDN) THEN
            ISTRT = ISTRT+NPT
            GO TO 150
C                  |
C      <-----------'
          ENDIF
          LPT    = ISTRT + NPT

C
C--->     Fill the arc banks.  Use helix parameters from the
C--->     chain closest to the origin.  1st arc only goes in TCAL.
C
          IF (NARC1.GT.1) THEN
            IF (CTYPE(IPICE).EQ.JTLNTP) THEN
              JCHAN=JCHA1
            ELSE
              JCHAN=JCH
            ENDIF
            NARCS=NARCS+1
C GAB:
C            KSARC=KROW(KTARC,NARCS)
c            IW(KSARC+JTARTN)=NCAND
c            IW(KSARC+JTAROF)=LTATC
c            IW(KSARC+JTARNC)=LPT
c            IW(KSARC+JTARAN)=NARC1*IORDR
c            DO 317 IPT =1,LPT
c              IW((KTATC+LMHLEN+LTATC)+IPT) = IBUFF(IPT)
c  317       CONTINUE

            MY_JTARTN(NARCS)=NCAND
            MY_JTAROF(NARCS)=LTATC
            MY_JTARNC(NARCS)=LPT
            MY_JTARAN(NARCS)=NARC1*IORDR
            DO 317 IPT =1,LPT
              MY_TATC(LTATC+IPT) = IBUFF(IPT)
  317       CONTINUE

            LTATC=LTATC+LPT
C GAB:
C            IW(KSARC+JTARPI)=IBMAP
            DO I=1,N32BITREG
             MY_JTARPI(NARCS,I)=IBMAP(I)
            ENDDO
C
C------->   Convert from TASSO sign convention (used in TCHA) to
C           the more standard convention (as used in TGFT, FRFT)
C
C GAB:
c            RW(KSARC+JTARIR)=-RTABL(KTCHA,JCHAN,JTCHIR)
c            RW(KSARC+JTARTL)= RTABL(KTCHA,JCHAN,JTCHTL)
c            RW(KSARC+JTARP0)= RTABL(KTCHA,JCHAN,JTCHP0)
c            RW(KSARC+JTARD0)= SIGN(1.0,RW(KSARC+JTARIR))
c     &                           * RTABL(KTCHA,JCHAN,JTCHD0)
c            RW(KSARC+JTARZ0)= RTABL(KTCHA,JCHAN,JTCHZ0)
c            RW(KSARC+JTARC1)= RTABL(KTCHA,JCHAN,JTCHC1)
c            RW(KSARC+JTARC2)= RTABL(KTCHA,JCHAN,JTCHC2)

            MY_JTARIR(NARCS)=-MY_JTCHIR(JCHAN)
            MY_JTARTL(NARCS)= MY_JTCHTL(JCHAN)
            MY_JTARP0(NARCS)= MY_JTCHP0(JCHAN)
            MY_JTARD0(NARCS)= SIGN(1.0,MY_JTARIR(NARCS))
     &                           * MY_JTCHD0(JCHAN)
            MY_JTARZ0(NARCS)= MY_JTCHZ0(JCHAN)
            MY_JTARC1(NARCS)= MY_JTCHC1(JCHAN)
            MY_JTARC2(NARCS)= MY_JTCHC2(JCHAN)

            LMOFF = 0
            DO 305 LM=1,5
              LMOFF = LMOFF + LM
C GAB:
c              RW(KSARC+JTARER-1+LM)=RTABL(KTCHA,JCHAN,JTCHER-1+LMOFF)
              MY_JTARER(LM)=MY_JTCHER(JCHAN,LMOFF)
  305       CONTINUE
          ENDIF
          NARC1=NARC1+1
C
          ISTRT = 0
          IF (IORDR.EQ.1) THEN
            DO 310 IPT =1,LPT
C GAB:
c              IW(NPUT+IPT) = IBUFF(IPT)
               IF(NPUT+IPT.LE.NTPCOMX)THEN
                MY_JTCCIT(NPUT+IPT)=IBUFF(IPT)
               ENDIF
  310       CONTINUE
          ELSE
            DO 320 IPT =1,LPT
C GAB:
C              IW(NPUT+IPT) = IBUFF(1+LPT-IPT)
              IF(NPUT+IPT.LE.NTPCOMX)THEN
               MY_JTCCIT(NPUT+IPT)= IBUFF(1+LPT-IPT)
              ENDIF
  320       CONTINUE

          END IF

          NPUT= NPUT + LPT
          NSTOR= NSTOR+ LPT
          IORDR=-IORDR
          JCHA1= IPLNK(IPICE+1)
          DO I=1,N32BITREG
           IBMAP(I)= 0
          ENDDO

        IF (IPICE.LT.NJCH)         GO TO 150
C                                        |
C      <---------------------------------'
C-----> Now fill 'TCAL' columns 2 and 3

C
C GAB:
C        N1ARC = ITABL(KTCHA,IPLNK(1),JTCHNC)
        N1ARC = MY_JTCHNC(IPLNK(1))
        DO 330 I=1,NJCH
          IF (CTYPE(I).NE.JTLNUP .AND. CTYPE(I).NE.JTLNDN) GO TO 331
C GAB:
c          N1ARC = N1ARC+ITABL(KTCHA,IPLNK(I+1),JTCHNC)
          N1ARC = N1ARC+MY_JTCHNC(IPLNK(I+1))
  330   CONTINUE
  331   CONTINUE

C GAB:
c        IW(KSCAL+JTCANF)=N1ARC
c        IW(KSCAL+JTCANR)=NSTOR-IW(KSCAL+JTCAIO)-N1ARC
        MY_JTCANF(NCAND)=N1ARC
        MY_JTCANR(NCAND)=NSTOR-MY_JTCAIO(NCAND)-N1ARC

  400 CONTINUE

C-------------------------------------------------------------------C
C              TREATMENT OF REMAINING CHAINS                        C
C              *****************************                        C
C-------------------------------------------------------------------C
      DO 500 ICHAN = 1,NCHAN


C GAB:
c        KSLNK = KROW(KTLNK,ICHAN)
        IF (TAKEN(ICHAN))              GO TO 500

        NCAND = NCAND + 1
C GAB:
c        NPT    = ITABL(KTCHA,ICHAN,JTCHNC)
        NPT    = MY_JTCHNC(ICHAN)

        KSCHC = KKTCTC(ICHAN)
C GAB:
c        KSCAL = KROW(KTCAL,NCAND)

C GAB:
c        IW(KSCAL+JTCAIO) = NSTOR
c        IW(KSCAL+JTCANF) = NPT   1

c        IW(KSCAL+JTCANR) = 0
         MY_JTCAIO(NCAND)=NSTOR
         MY_JTCANF(NCAND)=NPT
         MY_JTCANR(NCAND)=0

        DO 490 IPT =1,NPT
C GAB:
C          IW(NPUT+IPT) = IW(KSCHC+IPT)
           IF(NPUT+IPT.LE.NTPCOMX)THEN
            MY_JTCCIT(NPUT+IPT)=MY_TCTC(KSCHC+IPT)
           ENDIF

  490   CONTINUE
        NPUT   = NPUT + NPT
        NSTOR  = NSTOR+ NPT
  500 CONTINUE

C GAB:
cC
cC---> Define correct lengths of banks
cC
c      IW(KTTCC+LMHROW) = NSTOR
c      IW(KTCAL+LMHROW) = NCAND
c      IW(KTARC+LMHROW) = NARCS
c      IW(KTATC+LMHROW) = LTATC
c      CALL AUBPRS('TCALTARCTATC')
      GO TO 994
C
C  995 CONTINUE
C GAB:
C      CALL RERROR('TFCAND',4,'No room to book workbank')
C      IER=4
C      CALL BDROP(IW,'TCALTTCCTARCTATC')
C      GO TO 994
  998 CONTINUE
C     input banks not found
      IER=2
      GO TO 994
  999 CONTINUE
C GAB:
C      CALL RERROR('TFCAND',3,'Not enough space to create banks')
C      IER=3
C      CALL BDROP(IW,'TCALTTCCTARCTATC')
  994 CONTINUE
C
C---> Release the BOS space used by the workbank
C
C GAB:
C      CALL WDROP(IW,ITCALW)
      RETURN
      END

c-------------------------------------------------------------------
