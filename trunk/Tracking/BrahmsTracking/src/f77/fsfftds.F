*  
* $Id: fsfftds.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      LOGICAL FUNCTION FSFFTDS(IDAT)
********************************************************************
* selection function to find FTD singlepoints.
*
* INPUT:  IDAT   - TE structure
*                                         written by   Markus Elsing
*                                         changes by Kristian Harder
********************************************************************
      IMPLICIT NONE
#include "include/fsfparams.inc"
      INTEGER IDAT(*)
*
      IF (IDAT(1).GE.I_FTD1.AND.IDAT(1).LE.I_FTD7) THEN
        FSFFTDS=.TRUE.
      ELSE
        FSFFTDS=.FALSE.
      END IF
      END

