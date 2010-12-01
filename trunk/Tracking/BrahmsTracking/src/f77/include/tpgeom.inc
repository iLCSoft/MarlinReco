*  
* $Id: tpgeom.inc,v 1.3 2008-01-23 17:22:22 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.2  2005/08/08 07:09:13  aplin
* Made f77 tracking code use GEAR to define the geomtery of the TPC. LTPDRO now defines the maximum number of rows is used to define the size of arrays, this is limited to 224 due the use of 7 '32 bit' bit registers in trkfnd.F increased, though at present it is not likely that anybody would want more. The number of TPC padrows is defined at run time by NRTPC which should of course not exceed LTPDRO, although this is checked and the programe exits with a verbose error message. A wrapper function gettpcgeom is used to pass the GEAR TPC parameters from C++ to f77. MarlinUtil/include/marlin_tpcgeom.h have MarlinUtil/src/marlin_tpcgeom.cc consequently been removed as they are no longer needed.
*
* Revision 1.1.1.1  2005/07/11 15:15:48  aplin
* Initial version
*
* Revision 1.1.1.1  2005/06/18 10:05:36  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:29  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:48  hvogt
* Brahms V308
*
*  

#include "padrow.inc"

C      INTEGER LTPDRO, LTTROW, LTSROW, LTWIRE, LTSTYP, LTSLOT, LTCORN,
C     +        LTSECT, LTTPAD, LMXPDR, LTTSRW

      INTEGER LTTROW, LTSROW, LTWIRE, LTSTYP, LTSLOT, LTCORN,
     +        LTSECT, LTTPAD, LMXPDR, LTTSRW


C      PARAMETER (LTPDRO=256,LTTROW=19,LTSROW=12,LTWIRE=200,LTSTYP=3,
C     +           LTSLOT=12,LTCORN=6,LTSECT=LTSLOT*LTSTYP,LTTPAD=4,
C     +           LMXPDR=150,LTTSRW=11)

      PARAMETER (LTTROW=19,LTSROW=12,LTWIRE=200,LTSTYP=3,
     +           LTSLOT=12,LTCORN=6,LTSECT=LTSLOT*LTSTYP,LTTPAD=4,
     +           LMXPDR=150,LTTSRW=11)

      INTEGER NTSECT, NTPCRN, ITPTYP, ITPSEC, IENDTP
      REAL RTPCMN, DRTPMN, DRTPMX, DZTPMX, TPFRDZ,
     &     TPFRDW, TPAVDZ, TPFOF1, TPFOF2, TPFOF3, TPPROW, TPTROW,
     &     TPCORN, TPPHI0, TPCPH0, TPSPH0
      COMMON /TPGEOM/RTPCMN,DRTPMN,DRTPMX,DZTPMX,
     &               TPFRDZ,TPFRDW,TPAVDZ,TPFOF1,TPFOF2,TPFOF3,
     &               TPPROW(LTPDRO),TPTROW(LTTROW),NTSECT,
     &               NTPCRN(LTSTYP),TPCORN(2,LTCORN,LTSTYP),
     &               TPPHI0(LTSECT),TPCPH0(LTSECT),TPSPH0(LTSECT),
     &               ITPTYP(LTSECT),ITPSEC(LTSECT),IENDTP(LTSECT)

