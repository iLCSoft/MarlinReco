*  
* $Id: tpcfunc2.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:36  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:29  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:49  hvogt
* Brahms V308
*
*  
#include "brlink.inc"
      REAL TPCHITS

* SJA ahit and property used for the C++ interface
      INTEGER ATTRIBUTE, AHIT
* SJA stfreadtpchits will be used as the statement function agrument
      REAL STFREADTPCHITS

* FIXME: Replaced zebra bank with array for now with maximum size NTPCOMX=50000.
*        Should probably be done as stl vector in MARLIN
c      INTEGER HITLNK

      INTEGER HITLNK(NTPCOMX)


      INTEGER NTPHITS
      INTEGER ITPCOFF, ITPCSET
      INTEGER NTPC
      PARAMETER (NTPC = 7 )

      COMMON /HITLNKCOM/ HITLNK
      SAVE   /HITLNKCOM/

      COMMON / CNTPC / NTPHITS
      SAVE   / CNTPC /

