*  
* $Id: fsfand.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      SUBROUTINE FSFAND(P1,P2,P3,I)
*******************************************************************************
*
*
*
*     calculate smallest difference between two angels
*
*
*
*******************************************************************************
*
cd MVL can leads to infinite loop
c#include "include/gconst.inc"
#include "include/alcons.inc"

      I  = 0
      P3 = P1-P2
   10 IF(ABS(P3).LE.PI) GOTO 999
        I  = I+NINT(SIGN(1.,P3))
        P3 = P3-SIGN(TWOPI,P3)
      GOTO 10
  999 CONTINUE
      RETURN
      END

*   all hit reconstruction routines of the tracking detectors are
*   collected in this patch


