*  
* $Id: tpcdigi.F,v 1.2 2005-08-08 07:09:13 aplin Exp $
*  
* $Log: not supported by cvs2svn $
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
* Revision 1.2  2004/07/22 16:34:01  aplin
* changed random number generator from RANMAR to GRNDM
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      SUBROUTINE TPCDIGI(ITPTRA,ITPDGT,IER)
C****************************************************************************
c* Version number: 0.1      Date of last modification   6.11.99             *
c*------------------------------------------------------------------------- *
C* Routine to superpose random hits with occupancy TPC_OCC                  *
c* List of modifications since last version:                                *
C* Main Authors:                                                            *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
      IMPLICIT NONE
#include "include/fndparms.inc"
#include "include/tpcdims.inc"
#include "include/alcons.inc"
#include "include/tpgeom.inc"
c#include "marlin_tpcgeom.h"
#include "include/detnums.inc"
c#include "include/tpco.inc"
#include "include/tpcres.inc"
#include "include/tphtcns.inc"
c#include "gcflag.inc"
c#include "bkgmode.inc"
c#include "gcbank.inc"

#include "include/tpcfunc2.inc"


      INTEGER IER,ITPTRA(*),ITPDGT(*)

      INTEGER IROW,NPHI_PIX(LTPDRO),NZ_PIX,IHIT
      INTEGER IPHI,IZ,IPIX

      REAL RROW(LTPDRO),PIX_R,PIX_RP2,R,PHI,MAXAMP
      LOGICAL FIRST,LNEW

      INTEGER TPCIHIT(NTPMAX),N_ROW,I_ROW,IH_1,IH_2

      DATA FIRST/.TRUE./

      INTEGER  IP,NPIX,PHPX(NTPCOMX),ROWPX(NTPCOMX),ZPX(NTPCOMX),I
      INTEGER HITPX(NTPCOMX),IH1,IHIT1,IHIT2,IP1,IP2,IMGDWTH(NTPCOMX)
      REAL CHTOT,RVEC(3),RHIT,ZHIT,PHIHIT,HPOS(3)

      LOGICAL LMERGED,LGO

      REAL TKTPCR !(KH)
      
      external writetpccpp
      integer writetpccpp
      real tvalue
      integer tattribute, thit
      integer try

      SAVE

#include "include/tpcfunc.inc"

C------------------------------------------------------
      IER=0

      IF(FIRST) THEN

CCC FIXME: this is normally done at the end of tpcini.F
C TPC readout pixel size (cm) (from CDR)
       PIX_RP=0.33
       PIX_Z=1.4


       FIRST=.FALSE.
CCC set TPC_OCC to 0.01 normally set in cardsfile
       TPC_OCC = 0.01
       IF((TPC_OCC.LT.0).OR.(TPC_OCC.GE.1.))THEN
        WRITE(*,*)' TPCDIGI: Invalid occupancy'
        IER=1
        RETURN
       ENDIF
       IF(NRTPC.GE.LTPDRO )THEN
          WRITE(*,*)'TPCDIGI: Inconsistent number of TPC padrows'
          WRITE(*,*)'TPCDIGI: NRTPC greater than LTPDRO'
          WRITE(*,*)'TPCDIGI: Increase LTPDRO and run gmake clean bin'
          WRITE(*,*)'TPCDIGI: LTPDRO = ' , LTPDRO
          WRITE(*,*)'TPCDIGI: Presently LTPDRO limited to a maximum 
     +         value of 224'
          IER=2
          STOP
       ENDIF
       PIX_R=TPCPADR
       PIX_RP2=PIX_RP**2
C Define an index to divide up the TPC volume into 'pixels'
       DO IROW=1,NRTPC
        RROW(IROW)=TPCACRI +(REAL(IROW)-0.5)*TPCPADR
        NPHI_PIX(IROW)=INT(TWOPI*RROW(IROW)/PIX_RP)
        NZ_PIX=INT(2.*ZDRIFT/PIX_Z)
       ENDDO
      ENDIF

CSJAWRT
c      write(*,*) "the number of hits in tpcdigi = " , NTPHITS

C Assume initially that there is no merging:
      DO IHIT=1,NTPHITS
         ITPDGT(IHIT)=ITPTRA(IHIT)
      ENDDO
      
C     Translate each physics hit into a 'pixel' hit:
      NPIX=0
      
      DO IHIT=1,NTPHITS

         R=SQRT(TPCHITS(1,IHIT)**2+TPCHITS(2,IHIT)**2)
         PHI=ATAN2(TPCHITS(2,IHIT),TPCHITS(1,IHIT))
         IF(PHI.LT.0.)PHI=PHI+TWOPI
         IROW=1+INT((R-TPCACRI)/TPCPADR)
         IF (IROW.LT.1) IROW=1
         IF (IROW.GT.NRTPC) IROW=NRTPC
         IPHI=1+INT(REAL(NPHI_PIX(IROW))*PHI/TWOPI)
         IZ=1+INT(REAL(NZ_PIX)*(ZDRIFT+TPCHITS(3,IHIT))/(2.*ZDRIFT))
         LNEW=.TRUE.
         DO IP=1,NPIX
            IF((PHPX(IP).EQ.IPHI).AND.
     &           (ROWPX(IP).EQ.IROW).AND.(ZPX(IP).EQ.IZ))THEN
               LNEW=.FALSE.
               IPIX=IP
            ENDIF
         ENDDO
         IF(LNEW)THEN
            NPIX=NPIX+1
            IPIX=NPIX
            PHPX(NPIX)=IPHI
            ROWPX(NPIX)=IROW
            ZPX(NPIX)=IZ
         ENDIF
         HITPX(IHIT)=IPIX
      ENDDO
      
C Generate background hits:
      DO IPIX=1,NPIX
         
***   SJA*** Use GRNDM instead
*     CALL RANMAR(RVEC,1)
c         CALL GRNDM(RVEC,1)     
c     use ranecu as grndm needs geant
            CALL RANECU(RVEC,1,0)

         IF(RVEC(1).LT.TPC_OCC)THEN
            IPHI=PHPX(IPIX)
            IROW=ROWPX(IPIX)
            IZ=ZPX(IPIX)
C We have a hit pixel, so generate a random point within the pixel:

***SJA*** Use GRNDM instead  
*        CALL RANMAR(RVEC,3)
C            CALL GRNDM(RVEC,3)
c     use ranecu as grndm needs geant
            CALL RANECU(RVEC,3,0)
            LGO=.TRUE.
            IHIT=0
            DO WHILE(LGO.AND.(IHIT.LT.NTPHITS))
               IHIT=IHIT+1
               IF (HITPX(IHIT).EQ.IPIX)LGO=.FALSE.
            ENDDO
            RHIT=RROW(IROW)+(RVEC(1)-0.5)*TPCPADR
            ZHIT=ZDRIFT*
     &           ((REAL(IZ)+RVEC(2)-1.)/REAL(NZ_PIX)-0.5)
            PHIHIT=(REAL(IPHI)+RVEC(3)-1.)
     &           *TWOPI/REAL(NPHI_PIX(IROW))


 
CCCc  Replace the original hit with a 'saturation' pulse
CCC        Q(LBRTP1 + 7*(IHIT-1)+1+1) = RHIT*COS(PHIHIT)      ! x
CCC        Q(LBRTP1 + 7*(IHIT-1)+1+2) = RHIT*SIN(PHIHIT)      ! y
CCC        Q(LBRTP1 + 7*(IHIT-1)+1+3) = ZHIT                  ! z
CCC        Q(LBRTP1 + 7*(IHIT-1)+1+4) = TPCHITS(4,IHIT)+TPSAT ! energy
CCC        Q(LBRTP1 + 7*(IHIT-1)+1+5) = REAL(ID_TPC)          ! detector ID
CCCc(kh)   Q(LBRTP1 + 7*(IHIT-1)+1+6) = TPRP                  ! rphi resolution
CCC        Q(LBRTP1 + 7*(IHIT-1)+1+6) = TKTPCR(ZHIT)          ! rphi resolution
CCC        Q(LBRTP1 + 7*(IHIT-1)+1+7) = TPZ                   ! z resolution
CCCc Update the corresponding hit in TK banks
CCC        CALL TKUPDA('HIT',ID_TPC,IHIT,4,Q(LBRTP1 + 7*(IHIT-1)+1+1))

c            try = writetpccpp((RHIT*COS(PHIHIT)),1,IHIT1)
c            try = writetpccpp((RHIT*SIN(PHIHIT)),2,IHIT1)
c            try = writetpccpp(ZHIT,3,IHIT1)
c            try = writetpccpp(-CHTOT,4,IHIT1) 
c            try = writetpccpp(REAL(ID_TPC),5,IHIT1)
c            try = writetpccpp(TKTPCR(ZHIT),6,IHIT1)
c            try = writetpccpp(TPZ,7,IHIT1)


         ENDIF
      ENDDO

c loop over all rows here
      DO I_ROW = 1,NRTPC
c get all hits for this row
	N_ROW = 0
	CALL VZERO (TPCIHIT,NTPMAX)
        DO IHIT = 1,NTPHITS
           IF (ROWPX(HITPX(IHIT)).EQ.I_ROW) THEN
              N_ROW = N_ROW + 1
              TPCIHIT(N_ROW) = IHIT
           ENDIF
	ENDDO

        
c Look to merge neighbouring ('overlapping') hits, as defined
c by the effective pixel dimensions.
	CALL VZERO(IMGDWTH,NTPHITS)
	DO IH_1=1,N_ROW-1
           IHIT1 = TPCIHIT(IH_1)

           IF(ABS(TPCHITS(4,IHIT1)).GT.TPTHRS)THEN
              DO IH_2=IH_1+1,N_ROW
                 IHIT2 = TPCIHIT(IH_2)

                 LMERGED=.FALSE.
                 IF(ABS(TPCHITS(4,IHIT2)).GT.TPTHRS)THEN
                    IP1=HITPX(IHIT1)
                    IP2=HITPX(IHIT2)
                    IF(IP1.EQ.IP2)THEN
c     Merge any hits within the same pixel:
                       LMERGED=.TRUE.
                    ELSE IF (ROWPX(IP2).EQ.ROWPX(IP1)) THEN
                       IROW=ROWPX(IP1)
                       IF (ABS(ZPX(IP2)-ZPX(IP1)).LE.1)THEN
                          IF ((ABS(PHPX(IP2)-PHPX(IP1)).LE.1).OR.
     &                         (ABS(PHPX(IP2)-PHPX(IP1)).EQ.
     &                         (NPHI_PIX(IROW)-1)))THEN
C     Check whether hits within the same padrow with adjacent phi and z
c     should be merged:
                             IF(ABS(TPCHITS(3,IHIT1)-TPCHITS(3,IHIT2))
     &                            .LT.PIX_Z)THEN
                            IF(((TPCHITS(1,IHIT1)-TPCHITS(1,IHIT2))**2+
     &                          (TPCHITS(2,IHIT1)-TPCHITS(2,IHIT2))**2)
     &                          .LT.PIX_RP2) then 
                                   LMERGED=.TRUE.
                                   write(*,*) "found a merged hit"
                                ENDIF
                             ENDIF
                          ENDIF
                       ENDIF
                    ENDIF
                 ENDIF
                 IF(LMERGED)THEN
                    IF(IMGDWTH(IHIT1).NE.0)THEN
                       IH1=IMGDWTH(IHIT1)
                    ELSE
                       IH1=IHIT1
                    ENDIF
                    IMGDWTH(IHIT2)=IH1
                 ENDIF
              ENDDO
           ENDIF
        ENDDO
        
c Now actually merge the hits:
        DO IH_1=1,N_ROW-1
           IHIT1 = TPCIHIT(IH_1)
           
           LMERGED=.FALSE.
           CHTOT=0.
           CALL VZERO(HPOS,3)
           ITPDGT(IHIT1)=ITPTRA(IHIT1)
           MAXAMP=TPCHITS(4,IHIT1)
           
           DO IH_2=IH_1+1,N_ROW
              IHIT2 = TPCIHIT(IH_2)
              
c Assign the geant track to that hit with the highest amplitude:
              IF(IMGDWTH(IHIT2).EQ.IHIT1)THEN
                 LMERGED=.TRUE.
                 IF(TPCHITS(4,IHIT2).GT.MAXAMP)THEN
                    MAXAMP=TPCHITS(4,IHIT2)
                    ITPDGT(IHIT1)=ITPTRA(IHIT2)
                 ENDIF
                 CHTOT=CHTOT+TPCHITS(4,IHIT2)
c weight the position according to the relative charges of the adjacent hits:
                 DO I=1,3
                    HPOS(I)=HPOS(I)+ TPCHITS(4,IHIT2)*TPCHITS(I,IHIT2)
                 ENDDO
              ENDIF
           ENDDO
c Flag the resulting merged hit with a negative amplitude:
           IF(LMERGED)THEN
              CHTOT=CHTOT+TPCHITS(4,IHIT1)              
              DO I=1,3
                 HPOS(I)=HPOS(I)+ TPCHITS(4,IHIT1)*TPCHITS(I,IHIT1)
              ENDDO
              DO I=1,3
CCC   Q(LBRTP1+7*(IHIT1-1)+1+I) = HPOS(I)/CHTOT

C here is the place to try to write back to the bank
c     this will be removed for now to keep brahms and marlin consistent
c                 try = writetpccpp((HPOS(I)/CHTOT),I,IHIT1)
                 
              ENDDO

CCC   Q(LBRTP1+7*(IHIT1-1)+1+4) = -CHTOT
c     this will be removed for now to keep brahms and marlin consistent
c              try = writetpccpp(-CHTOT,4,IHIT1)

           ENDIF
        ENDDO
      ENDDO
      
c remove the absorbed merged hits by setting their amplitude to zero:
      DO IHIT2=1,NTPHITS
         IF(IMGDWTH(IHIT2).NE.0) THEN
CCC   Q(LBRTP1+7*(IHIT2-1)+1+4)=0.
c     this will be removed for now to keep brahms and marlin consistent
c            try = writetpccpp(0.,4,IHIT2)
C FIXME:SJA: TKUPDA is updating the tkhit bank for the merged hit. 
C     At present only hits from the tpc are used so this will be done later, 
C     IHPOINT and IHNUMB will be needed to as the hold index info on the hits in the bank  
CCC   CALL TKUPDA('HIT',ID_TPC,IHIT2,4,Q(LBRTP1 + 7*(IHIT2-1)+1+1))
         END IF
      ENDDO



      RETURN
      END
      
