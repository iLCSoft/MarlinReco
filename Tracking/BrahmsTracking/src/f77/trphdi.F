*  
* $Id: trphdi.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
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
      FUNCTION TRPHDI(PHI1,PHI2)
C
C! Difference between PHI1 ind PHI2 in range -pi to + pi
C  (PHI1, PHI2 are in range 0 - 2pi)
C
#include "include/alcons.inc"
C
      DPHI = PHI1 - PHI2
      IF(DPHI.GT. PI) DPHI =-TWOPI + DPHI
      IF(DPHI.LT.-PI) DPHI = TWOPI + DPHI
      TRPHDI = DPHI
      RETURN
      END

c---------------------------------------------------------------------------
