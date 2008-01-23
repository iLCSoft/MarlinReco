*  
* $Id: tcobin.F,v 1.4 2008-01-23 17:22:22 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.3  2005/08/09 10:27:50  aplin
* increased the maximum number of hits in one row to 1000 ( NMAXROW = 1000 in src/f77/include/my_tcrl.inc ) , and ensured the tracking is skipped in the event when this number of hits is exceeded (tcobin.F).
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
      SUBROUTINE TCOBIN(IER)
      IMPLICIT NONE
C****************************************************************************
c* Version number: 0.0      Date of last modification  26.08.99             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
C* Main modifiers from ALEPH to BRAHMS code:                                *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
C
C--------------------------------------------------------------------
C! Bin coordinates for pattern recognition
C!
C!    Author:  R. Johnson     4-12-87
C!
C!    Output:   IER     /I     Error return code=0 for good return
C!
C!    Called by TRKFND
C!
C!-------------------------------------------------------------------
C
c
#include "include/temp.inc"
#include "include/fndparms.inc"
#include "include/tpco.inc"
c#include "gcflag.inc"
#include "include/alcons.inc"
#include "include/tpgeom.inc"
c#include "marlin_tpcgeom.h"

#include "include/tcbins.inc"
#include "include/tpcdims.inc"
#include "include/my_tcrl.inc"
#include "include/my_jndb.inc"

#include "include/unmrg.inc"
      INTEGER IER,I,K,I2NMX,IROW,IHIT,NC,IPTR,II,IBN1,IBN2
      REAL RHIT,PHI
      DIMENSION IPTR(LMXBIN)
      LOGICAL FIRST
      SAVE

      DATA FIRST/.TRUE./
C
C---------------------- Statement Function Definitions --------------
C
C++   Initialize the common TCBINS on first call
C

      IF (FIRST) THEN
        FIRST=.FALSE.
        INDBIN=0
C
C++     Fill array used to wrap around the branch 0,2pi when calculating
C++     the bin number for a given value of phi
C
        DO 20 I=1,LMXBIN
          IXWRAP(I+2*NMXWRP+1)=I
   20   CONTINUE
        I2NMX=2*NMXWRP
        DO 30 I=1,I2NMX
          IXWRAP(I+LMXBIN+2*NMXWRP+1)=I
   30   CONTINUE
        IXWRAP(LNWRAP)=IXWRAP(LNWRAP-1)+1
        IF (IXWRAP(LNWRAP).GT.LMXBIN) THEN
          IXWRAP(LNWRAP)=IXWRAP(LNWRAP)-LMXBIN
        ENDIF
        K=LMXBIN+1
        DO 40 I=I2NMX,1,-1
          K=K-1
          IXWRAP(I+1)=K
   40   CONTINUE
        IXWRAP(1)=IXWRAP(2)-1
        IF (IXWRAP(1).LT.1) IXWRAP(1)=IXWRAP(1)+LMXBIN
C
        TCBINW= 2.*TWOPI/FLOAT(LMXBIN)
        TCDELT= TCBINW/4.
        TC2DEL= 2.*TCDELT
        TCOFFS= 1.5 + FLOAT(2*NMXWRP)
      ENDIF

C GAB:
      IF(LUNMRG)THEN
       MINROW=MINMRG
      ELSE
       MINROW=1
      ENDIF

      IF (NTPCO.EQ.0) THEN
        WRITE(*,*)'TCOBIN: no tpc hits'
        IER=2
        GO TO 999
      ENDIF

      IBLAST=1
      CALL VZERO(IBNPTR,LMXBIN*NRTPC)

      DO IROW=1,NRTPC
       MY_JTCRNC(IROW)=0
       DO IHIT=1,NMAXROW
        IHITROW(IROW,IHIT)=0
       ENDDO
      ENDDO

      NC=0
      DO IHIT=1,NTPCO
c     GAB: allow for hits to be removed by a merge-finder:
         IF(.NOT.LBLOT(IHIT))THEN
            RHIT=TPCO_RV(IHIT)
            IROW=1+INT((RHIT-TPCACRI)/TPCPADR)
            IF (IROW.LT.1) IROW=1
            IF (IROW.GT.NRTPC) IROW=NRTPC
C     GAB:
            IF(IROW.GE.MINROW)THEN
               MY_JTCRNC(IROW)=MY_JTCRNC(IROW)+1
               IF(MY_JTCRNC(IROW).LE.NMAXROW)THEN
                  IHITROW(IROW,MY_JTCRNC(IROW))=IHIT
               ELSE
                  WRITE(*,*)'TPCOBIN: NMAXROW too small ',
     &        'This is FATAL: TPC Tracking in event Skipped  !',                   
     &        'Increase NMAXROW in my_tcrl.inc and recomplile'
                  CALL TKVETO(1)
                  IER = 1
                  GO TO 999
               ENDIF
            ENDIF
         ENDIF
      ENDDO

C GAB:
      DO 900 IROW=MINROW,NRTPC
        CALL VZERO(IPTR,LMXBIN)
C
C++     Loop over all coordinates in the padrow
C
        NC=MY_JTCRNC(IROW)


        DO 800 II=1,NC

           IHIT=IHITROW(IROW,II)
           PHI=TPCO_PH(IHIT)
C
C++       For each coordinate, fill pointers into two overlapping bins
C
          IBN1= IXWRAP(2*(INT(PHI/TCBINW)+NMXWRP)+2)

          IF (IPTR(IBN1).EQ.0) THEN

            IBNPTR(IBN1,IROW)=IBLAST
C GAB:
               MY_JNDBBW(IBLAST)=-IBN1
          ELSE
C GAB:
               MY_JNDBFW(IPTR(IBN1))=IBLAST
               MY_JNDBBW(IBLAST)=IPTR(IBN1)
          ENDIF
C GAB:
          MY_JNDBCO(IBLAST)=IHIT
          MY_JNDBFW(IBLAST)=0
          IPTR(IBN1)=IBLAST
C
C++       Note that the two entries for a single coordinate always lie
C++       adjacent to each other in workbank INDBIN.
C
          IBLAST=IBLAST+1

C GAB:
          IF(IBLAST.GT.NTPCOMX)THEN
           WRITE(*,*)'TCOBIN: NTPCOMX too small. ', iblast,
     &        'This is FATAL: Tracking in event Skipped !'
           CALL TKVETO(1)
           IBLAST=IBLAST-1
           IER = 1
           GO TO 999
          ENDIF

          IBN2= IXWRAP(2*(INT((PHI/TCBINW)+0.5)+NMXWRP)+1)

          IF (IPTR(IBN2).EQ.0) THEN
            IBNPTR(IBN2,IROW)=IBLAST
C GAB:
            MY_JNDBBW(IBLAST)= -IBN2
          ELSE
C GAB:
           MY_JNDBFW(IPTR(IBN2))=IBLAST
C GAB:
           MY_JNDBBW(IBLAST)=IPTR(IBN2)
          ENDIF

C GAB:
          MY_JNDBCO(IBLAST)=IHIT
          MY_JNDBFW(IBLAST)=0
          IPTR(IBN2)=IBLAST
          IBLAST=IBLAST+1
C GAB:
          IF(IBLAST.GT.NTPCOMX)THEN
           WRITE(*,*)'TCOBIN: NTPCOMX = ',iblast,' too small.',
     &        'This is FATAL: Tracking in event Skipped !'
           CALL TKVETO(1)
           IBLAST=IBLAST-1
           IER = 1
           GO TO 999
        ENDIF

  800   CONTINUE
  900 CONTINUE
      IER=0
      GO TO 999

  999 CONTINUE
      RETURN
      END

C-------------------------------------------------------------------------

