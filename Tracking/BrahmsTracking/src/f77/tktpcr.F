*  
* $Id: tktpcr.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      FUNCTION TKTPCR(Z)
C****************************************************************************
C* PARAMETRIC RPHI RESOLUTION FOR TPC HITS                                  *
C*                                                                          *
C* Kristian Harder, February 2001                                           *
C****************************************************************************
      IMPLICIT NONE

      REAL Z,TKTPCR
#include "include/tpcdims.inc"

      TKTPCR=0.019-ABS(Z)/ZDRIFT*0.010

      RETURN
      END


