*  
* $Id: tpcset.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
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
      SUBROUTINE TPCSET
C****************************************************************************
c* Version number: 0.0      Date of last modification  26.08.99             *
c*------------------------------------------------------------------------- *
c* Set up the various constants required for the pattern recognition        *
c* List of modifications since last version:                                *
C* Main Authors:                                                            *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
      IMPLICIT NONE
#include "include/fndparms.inc"
#include "include/alcons.inc"
#include "include/coildims.inc"
#include "include/rconds.inc"
c#include "include/tpgeom.inc"
#include "marlin_tpcgeom.h"
#include "include/tpccon.inc"
#include "include/tpcdims.inc"
#include "include/tpgeop.inc"
#include "include/strcons.inc"
#include "include/rflags.inc"
#include "include/tphtcns.inc"
#include "include/tpcmatpar.inc"

      INTEGER I

C Set the field to be the same as in Brahms:
      FIELRC=BFIELD
C Convert the B-field to natural units
      BCFGEV = BFIELD*CLGHT*1.E-5

c convert radial pad spacing:
      DO I=1,LTSTYP
       TPDRST(I)=TPCPADR
      ENDDO

c Enter the tpc row radii
      write(*,*) "TPCACRI = " , TPCACRI
      write(*,*) "TPCPADR = " , TPCPADR

      TPPROW(1)=TPCACRI+0.5*TPCPADR
      DO I=2,LTPDRO
       TPPROW(I)=TPPROW(1)+(I-1)*TPCPADR
      ENDDO

C TPC material parameters for track fit (KH)
      RADLEN=(TPCACRO-TPCACRI)/10971.
      RMS=29.6

      RETURN

      END

C----------------------------------------------------------------------

