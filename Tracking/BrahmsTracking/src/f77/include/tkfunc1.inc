*
* $Id: tkfunc1.inc,v 1.3 2005-07-21 16:37:00 aplin Exp $
*
* $Log: not supported by cvs2svn $
* Revision 1.2  2005/07/20 17:46:38  aplin
* The itedat bank has been added with read and write interface. Though
* only the write function is used (tkstor). The read will have to be
* done through a statement function.
*
* Revision 1.1.1.1  2005/07/11 15:15:48  aplin
* Initial version
*
* Revision 1.1.1.1  2005/06/18 10:05:35  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:29  aplin
* *** empty log message ***
*
* Revision 1.7  2004/07/22 13:52:00  aplin
* specification statements added separately for statement functions in tkfunc.inc
*
*



#include "brlink.inc"
      REAL RTKHIT, TKMCTR
      INTEGER ITKHIT, ITKMCT
      INTEGER IOFFSET, IHITOFFSET

*     TRACK TE variables
      REAL RTE
      INTEGER ITE, ITEDAT

* SJA ahit used for the C++ interface
      INTEGER ATTRIBUTE, AHIT
* SJA stfreadtkhits will be used as the statement function agrument
      REAL RSTFREADTKHITS
      INTEGER ISTFREADTKHITS

* SJA ate used for the C++ interface
      INTEGER ATE
* SJA stfreadtkhits will be used as the statement function agrument
      REAL RSTFREADTKTE
      INTEGER ISTFREADTKTE
      INTEGER STFREADTKITEDAT
	