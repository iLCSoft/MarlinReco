*  
* $Id: fsfftdm.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      LOGICAL FUNCTION FSFFTDM(IDAT)
********************************************************************
* selection function to find FTD multipoints.
*
* INPUT:  IDAT   - TE structure
*                                         written by   Markus Elsing
*                                         changes by Kristian Harder
********************************************************************
      IMPLICIT NONE
#include "include/fsfparams.inc"
      INTEGER IDAT(*)
*
      IF (IDAT(1).EQ.I_FTD) THEN
        FSFFTDM=.TRUE.
      ELSE
        FSFFTDM=.FALSE.
      END IF
      END

