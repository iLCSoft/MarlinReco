*  
* $Id: tflnwb.inc,v 1.2 2008-01-23 17:22:22 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/07/11 15:15:48  aplin
* Initial version
*
* Revision 1.1.1.1  2005/06/18 10:05:35  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:29  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:49  hvogt
* Brahms V308
*
*  

#include "padrow.inc"

C      INTEGER N32BITREG3
C      PARAMETER(N32BITREG3=8)

      INTEGER*4 NBPWRK
      COMMON/TFLNWB/NBPWRK(1000,N32BITREG)

