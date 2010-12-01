*
* $Id: tkfunc1.inc,v 1.5 2008-01-24 17:21:50 aplin Exp $
*
* $Log: not supported by cvs2svn $
* Revision 1.4  2008/01/23 17:22:22  aplin
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
* Revision 1.3  2005/07/21 16:37:00  aplin
* The itedat bank has been added with read and write interface.
*
* TKADDH has now been added so that the silicon hits are included as
* TE's.
*
* Added FSFSTR from delsearch to steer the creation of
* TS's
*
* The following include files need to be added for FSFSTR
* include/fsfdebug.inc
* include/fsfreject.inc
* include/fsfflags.inc
*
* FSFSTR is now called by TKTREV although all action is commented out.
*
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
	
