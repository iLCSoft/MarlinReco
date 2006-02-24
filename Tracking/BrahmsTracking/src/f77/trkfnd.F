
*  
* $Id: trkfnd.F,v 1.6 2006-02-24 16:36:55 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.5  2005/11/03 15:16:14  aplin
* Added the Trackstring creation and the biulding of full Track candiates (TK's) which have passed the Delphi Ambiguity resolver fxambi. The material description of the vtx detector, as for the TPC, is hard coded in setmat. Presently the VTX and SIT resolutions are hard coded in LEPTrackingProcessor. The debug output has been reduced and can be controlled via TKSTDBG etc. in tkinit.F. delsolve contains the delphi ambuguity resolver written in C and is contained in the directory named C. The Tk's are written back into the C++ side in tktrev. The corresponding Tk bank structure analogous to the TE bank structure has been added in tktkbank whilst the access wrapper functions are contained in LEPTracking.
*
* Revision 1.4  2005/08/08 07:09:13  aplin
* Made f77 tracking code use GEAR to define the geomtery of the TPC. LTPDRO now defines the maximum number of rows is used to define the size of arrays, this is limited to 224 due the use of 7 '32 bit' bit registers in trkfnd.F increased, though at present it is not likely that anybody would want more. The number of TPC padrows is defined at run time by NRTPC which should of course not exceed LTPDRO, although this is checked and the programe exits with a verbose error message. A wrapper function gettpcgeom is used to pass the GEAR TPC parameters from C++ to f77. MarlinUtil/include/marlin_tpcgeom.h have MarlinUtil/src/marlin_tpcgeom.cc consequently been removed as they are no longer needed.
*
* Revision 1.3  2005/08/04 13:17:48  aplin
* added more info to error message on too many links
*
* Revision 1.2  2005/08/04 12:54:51  aplin
* *** empty log message ***
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
      SUBROUTINE TRKFND(IER)
      IMPLICIT NONE
C****************************************************************************
c* Version number: 0.0      Date of last modification  26.08.99             *
c*------------------------------------------------------------------------- *
c* Main track finding routine                                               *
c* List of modifications since last version:                                *
C* Main modifiers from ALEPH to BRAHMS code:                                *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
C--------------------------------------------------------------------
C! Track-finding in the TPC
C! using 3-D chain method (ALEPH-TPCDET Note 84-69)
C!
C!   Author:     M.Mermikides    06-10-86
C!   Modified:   R.Johnson        8-04-88
C!   Modified:   M.Mermikides    19-08-89
C!
C!   Output:  IER    /I    Error code, nonzero if error has occurred
C!   Input banks:
C!            'TPCO'       TPC coordinates
C!            'TCRL'       TPC coordinate row list
C!   Output banks:
C!            'TCHA'       Chain parameters
C!            'TCTC'       Chain -> coordinate reference
C!
C!
C!   Called by TPCRUN, UNMERGE
C!
C!   Strategy:
C!   ========
C!
C!     Three-point links are formed from successive padrows,
C!  starting from the outermost row.(There can be a maximum gap
C!  of MAXGAP between successive rows). Links are sorted in order
C!  of decreasing radius of curvature and propagated inwards
C!  using a local point search to form longest possible chain.
C!  Chains are verified by a helix fit and fitted parameters
C!  are used to project incomplete chains both inwards and
C!  outwards to pick extra coordinates.
C!
C!  Workbanks used internally:
C!                INDCXY    X-Y coordinates corresp. to TPCO
C!                INDLNK    Padrow and hit references and
C!                          helix parameters of links from
C!                          current staring row
C!                JSORTW    Workspace for sorting links
C!                INDBIN    Date structure for binned coordinates
C!
C!  NOTE:  This routine and the subroutines which it calls represent
C!         a MODULE within JULIA.
C-------------------------------------------------------------------
C
#include "include/fndparms.inc"
#include "include/tpco.inc"
c#include "include/gcflag.inc"
#include "include/alcons.inc"
#include "include/tchajj.inc"
#include "include/tpwork.inc"
#include "include/trfdct.inc"
#include "include/tpgeom.inc"
c#include "marlin_tpcgeom.h"
#include "include/tchain.inc"
#include "include/tcbins.inc"
#include "include/my_jndb.inc"
#include "include/my_tcrl.inc"
#include "include/tpcdims.inc"
#include "include/my_link.inc"
#include "include/my_tctc.inc"
#include "include/my_tcha.inc"
C
C----------------------------------------------------------------------
C
C GAB:
      INTEGER NDCHA,NXCHA,NWCHC,NWCXY,MNLNK,NWLNK,NDLNK
      INTEGER I,J,K,IROW,NCO1,NCO2,NCO3,IBIN,IPTR,IRWM,IRWI,NL,NLNKS
      INTEGER KK,ILINK,NHITS,IC,IPTBW,IPTFW

      DATA NDCHA/800/, NXCHA/200/
      DATA NWCHC/1/,NWCXY/2/,MNLNK/7/
      DATA NWLNK/16/, NDLNK/300/
      LOGICAL LFRST
      DATA LFRST/.TRUE./

#include "include/unmrg.inc"

C
C----------------------------------------------------------------------
      INTEGER LNSRT,IER
      PARAMETER (LNSRT=100)

      REAL MY_WORK(NMAXLNK)
      INTEGER MY_ISORT(NMAXLNK)

      SAVE
c----------------------------------------------------------------------

      NLAST=0
      MY_NCHAIN=0
      DO I=1,NTPCOMX
       MY_TCTC(I)=0
      ENDDO

      DO I=1,NCHMAX
       DO J=1,7
        MY_JTCHPI(I,J)=0
       ENDDO
      ENDDO

      NLNK=0
C
      IF (LFRST) THEN
        LFRST=.FALSE.
        INDLNK=0
        INDCXY = 0
        JSORTW = 0
      ENDIF

C
      IER=0

C     check that input bank is present
      IF (NTPCO.EQ.0) GO TO 800

C
C  Bin the coordinates
C

      CALL TCOBIN(IER)
      IF (IER.NE.0) THEN
         WRITE(*,*)'TRKFND:Error in binning of coordinates'
         IER=5
         GO TO 999
      ENDIF

C  Start looking for 3-point links from outermost padrow.
C  Stop when we reach padrow 4.
C
C GAB: stop at padrow MINMRG when looking for merged hits:
C GAB:
      IF(LUNMRG)THEN
       MINROW=MINMRG
      ELSE
       MINROW=4
      ENDIF

      write(*,*) "MINROW = ", MINROW

      DO 30 IROW=NRTPC, MINROW, -1

         NCO1=MY_JTCRNC(IROW)

         IF(NCO1.EQ.0) GO TO 30
C
C  Loop over every other phi bin in this row (this will get all
C  coordinates on this row, since they are double binned).
C
         DO 20 IBIN=1,LMXBIN,2
C
C  Loop through the linked list of coordinates for this bin
C
            IPTR= IBNPTR(IBIN,IROW)

c(kh)  513       IF (IPTR.EQ.0) GO TO 514
  513       IF (IPTR.LE.0) GO TO 514
C
C  Test all combinations of 3 padrows with max gap MAXGAP between
C  successive ones.
C
               DO 25 J=1,MAXGAP+1
                  IRWM = IROW - J
                  IF (IRWM.LT.3) GO TO 27

                  NCO2 = MY_JTCRNC(IRWM)

                  IF (NCO2.EQ.0) GO TO 25
C
                  DO 26 K=1,MAXGAP+1
                     IRWI = IRWM - K
                     IF (IRWI.LT.2) GO TO 26

                     NCO3 = MY_JTCRNC(IRWI)

                     IF (NCO3.EQ.0) GO TO 26
C
C  Search for 3-point link(s) from this hit
C
                    CALL TRSLNK(IROW,IRWM,IRWI,IPTR,NL,IER)
                     IF (IER.NE.0) GO TO 999
C
C  Exit the double loop as soon as a valid link is found
C
C     SJAWRT
c                     write(*,*) "the number of links found = ", NL

                     IF (NL.GT.0) GO TO 27
   26             CONTINUE
   25          CONTINUE
   27          CONTINUE
C
C  Find the next coordinate in the linked list for this phi bin
C
               IPTR=MY_JNDBFW(IPTR)
               GO TO 513
  514       CONTINUE
   20    CONTINUE

C
C   End of links originating from current pad row.
C   Sort links in decreasing radius of curvature so we get
C   fast tracks first.  Extend workbank if necessary.

c GAB:
          NLNKS=NLNK

         IF(NLNKS.GT.NMAXLNK)THEN
          WRITE(*,*)'TRKFND: NLNKS too big;tracks may be lost'
          CALL TKVETO(1)
         ENDIF
         IF (NLNKS.GE.1) THEN
             DO 123 I = 1,NLNKS
              MY_WORK(I) = ABS(RLNK(I,7))
  123       CONTINUE

           CALL SORTZV(MY_WORK,MY_ISORT,NLNKS,1,1,0)

         ENDIF

C
C   Try to start a chain from each valid link
C
         DO 10 KK= 1,NLNKS
C GAB:
            ILINK=MY_ISORT(KK)
C
C   Copy the chain's parameters into the current chain (common/TCHAIN/)
C   If any hits in this link have already been used by a validated
C   chain, then skip it.
C
            DO I=1,3
             IHCHAI(I)=INT(RLNK(ILINK,3+I))
            ENDDO

            IF (MY_JNDBFW(IHCHAI(1)).LT.0) GO TO 10
            IF (MY_JNDBFW(IHCHAI(2)).LT.0) GO TO 10
            IF (MY_JNDBFW(IHCHAI(3)).LT.0) GO TO 10

            NHCHAI = 3

             DO I=1,3
              IRCHAI(I)=INT(RLNK(ILINK,I))
             ENDDO

             DO I=1,6
              PCHAIN(I)=RLNK(ILINK,6+I)
             ENDDO

            PCHAIN(JTCHC1) = 999.
            PCHAIN(JTCHC2) = 999.
            CALL VZERO(PCHAIN(JTCHER),15)

            D2CHAI(1) = RLNK(ILINK, 14)
            D2CHAI(2) = RLNK(ILINK, 15)
            D2CHAI(3) = RLNK(ILINK, 16)
            ISTCHA = 0
C

C   Move inward row-by-row to add new points to the chain
            CALL TRSEEK

C  If we have enough points fit the chain to validate.
            IF (NHCHAI.GE.MNLNK)THEN
             CALL TRFTCH
            ENDIF
C
C  If fit fails and there are enough points, try to rescue chain
C
            IF (ISTCHA.EQ.0 .AND. NHCHAI.GE.8) THEN
              CALL TRRESC
            ENDIF

C
C  Drop chain if unsuccessful; otherwise try to project across all padro
C  to pick up additional hits.
C
            IF (ISTCHA.EQ.0)         GO TO 10

C     SJAWRT
c            write(*,*) "number of hits in the chain = " , NHCHAI


C
            IF (NHCHAI.GE.MNLNK) CALL TRPROJ
            
               MY_NCHAIN=MY_NCHAIN+1
               MY_JTCHIR(MY_NCHAIN)=1./PCHAIN(1)
               MY_JTCHTL(MY_NCHAIN)=PCHAIN(JTCHTL)
               MY_JTCHP0(MY_NCHAIN)=PCHAIN(JTCHP0)
               MY_JTCHD0(MY_NCHAIN)=PCHAIN(JTCHD0)
               MY_JTCHZ0(MY_NCHAIN)=PCHAIN(JTCHZ0)
               MY_JTCHLC(MY_NCHAIN)=PCHAIN(JTCHLC)
               MY_JTCHC1(MY_NCHAIN)=PCHAIN(JTCHC1)
               MY_JTCHC2(MY_NCHAIN)=PCHAIN(JTCHC2)
               DO J=1,15
                  MY_JTCHER(MY_NCHAIN,J)=PCHAIN(8+J)
               ENDDO
               
               write(10,100) ABS(PCHAIN(1)), NHCHAI
               
 100           FORMAT(F10.3,I6)
               
               DO J=1,7
                  MY_JTCHPI(MY_NCHAIN,J)=0
               ENDDO
               

C
C   Remove associated coordinates from the phi bins, enter
C   hits in TCTC bank and define offset in TCHA.

            NHITS = 0

            DO 11 J=NHCHAI,1,-1
               IF (MCHAIN(J).EQ.0) THEN
                  IPTR= IHCHAI(J)

                   IF(IRCHAI(J).LE.32)THEN
                    MY_JTCHPI(MY_NCHAIN,1)=
     &               IBSET(MY_JTCHPI(MY_NCHAIN,1),IRCHAI(J)-1)
                   ELSE IF(IRCHAI(J).LE.64)THEN
                    MY_JTCHPI(MY_NCHAIN,2)=
     &               IBSET(MY_JTCHPI(MY_NCHAIN,2),IRCHAI(J)-33)
                   ELSE IF(IRCHAI(J).LE.96)THEN
                    MY_JTCHPI(MY_NCHAIN,3)=
     &               IBSET(MY_JTCHPI(MY_NCHAIN,3),IRCHAI(J)-65)
                   ELSE IF (IRCHAI(J).LE.128)THEN
                    MY_JTCHPI(MY_NCHAIN,4)=
     &               IBSET(MY_JTCHPI(MY_NCHAIN,4),IRCHAI(J)-97)
                   ELSE IF (IRCHAI(J).LE.160)THEN
                    MY_JTCHPI(MY_NCHAIN,5)=
     &               IBSET(MY_JTCHPI(MY_NCHAIN,5),IRCHAI(J)-129)
                   ELSE IF (IRCHAI(J).LE.192)THEN
                    MY_JTCHPI(MY_NCHAIN,6)=
     &               IBSET(MY_JTCHPI(MY_NCHAIN,6),IRCHAI(J)-161)
                   ELSE IF (IRCHAI(J).LE.224)THEN
                    MY_JTCHPI(MY_NCHAIN,7)=
     &               IBSET(MY_JTCHPI(MY_NCHAIN,7),IRCHAI(J)-193)
                   ELSE
                    WRITE(6,*) 'FATAL INTERNAL ERROR IN TRKFND:'
                    WRITE(6,*) ' TOO MANY PAD ROWS!'
                    STOP
                   ENDIF

                  NHITS = NHITS + 1
c GAB:
                   IC=MY_JNDBCO(IPTR)
                   MY_TCTC(NLAST+NHITS)=IC
C
C  Each coordinate must be removed from two phi bins
C
                  IPTBW= MY_JNDBBW(IPTR)
                  IPTFW= MY_JNDBFW(IPTR)

                  IF (IPTBW.GT.0) THEN
C GAB:
                     MY_JNDBFW(IPTBW)= IPTFW
                  ELSE
                     IBNPTR(-IPTBW,IRCHAI(J))= IPTFW
                  ENDIF

C GAB:
c(kh)             IF (IPTFW.NE.0)MY_JNDBBW(IPTFW)= IPTBW
                  IF (IPTFW.GT.0)MY_JNDBBW(IPTFW)= IPTBW

C GAB:
                  MY_JNDBFW(IPTR)= -1


C
C  The other entry for the same coordinate is adjacent to the first
C  in workbank INDBIN, but we do not know to which side.
C
                  IPTR= IPTR-1
                  IF (IPTR.NE.0) THEN
C GAB:
c(kh)               IF (MY_JNDBCO(IPTR).NE.IC) THEN
                    IF (MY_JNDBCO(IPTR).NE.IC.and.iptr.lt.ntpcomx-1)THEN
                      IPTR= IPTR+2
                    ENDIF
c(kh)             ELSE
                  ELSE if (iptr.lt.ntpcomx-1) then
                    IPTR= IPTR+2
                  ENDIF

C GAB:
                  IPTBW= MY_JNDBBW(IPTR)
                  IPTFW= MY_JNDBFW(IPTR)

                  IF (IPTBW.GT.0) THEN
C GAB:
                     MY_JNDBFW(IPTBW)= IPTFW
                  ELSE
                     IBNPTR(-IPTBW,IRCHAI(J))= IPTFW
                  ENDIF
C GAB
c(kh)             IF (IPTFW.NE.0) MY_JNDBBW(IPTFW)= IPTBW
                  IF (IPTFW.GT.0) MY_JNDBBW(IPTFW)= IPTBW
                  MY_JNDBFW(IPTR)= -1
               ENDIF
   11       CONTINUE

C GAB:
             MY_JTCHNC(MY_NCHAIN)=NHITS
             MY_JTCHOF(MY_NCHAIN)=NLAST
             NLAST=NLAST+NHITS

   10    CONTINUE
C
   30 CONTINUE


      GO TO 999
C
C Set error code and print error messages
C
  800 CONTINUE
      IER=1
      WRITE(*,*)'TPCFND: No TPC hits'
      GO TO 999


  999 CONTINUE


C      SJAWRT
c      write(*,*) "finished trkfnd with No. chains = ", MY_NCHAIN 

      RETURN
      END

C-----------------------------------------------------------------------

