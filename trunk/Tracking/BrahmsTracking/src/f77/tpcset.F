*  
* $Id: tpcset.F,v 1.3 2008-04-10 21:53:32 aplin Exp $
*  
* $Log: not supported by cvs2svn $
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
#include "include/tpgeom.inc"
c#include "marlin_tpcgeom.h"
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
c      write(*,*) "TPCACRI = " , TPCACRI
c      write(*,*) "TPCPADR = " , TPCPADR

      TPPROW(1)=TPCACRI+0.5*TPCPADR
      DO I=2,NRTPC
       TPPROW(I)=TPPROW(1)+(I-1)*TPCPADR
      ENDDO

C TPC material parameters for track fit (KH)
C     FIXME:SJA The material used for the RADLENGTH is Hard coded here !
      RADLEN=(TPCACRO-TPCACRI)/10971.
      RMS=29.6

      RETURN

      END

C----------------------------------------------------------------------

