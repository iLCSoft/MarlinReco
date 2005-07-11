*  
* $Id: deldet.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:18  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      INTEGER FUNCTION DELDET(BRADET)
C****************************************************************************
C* CONVERT BRAHMS DETECTOR ID TO DELPHI CONVENTION FOR TK2FIT CALL          *
C*                                                                          *
C* Kristian Harder, October 1999                                            *
C****************************************************************************
      IMPLICIT NONE
#include "include/detnums.inc"
#include "include/fsfparams.inc"

      INTEGER BRADET,DELID,ID

      ID=DELID(BRADET)
      IF (ID.GT.0.AND.ID.LE.DETMAX) THEN
        DELDET=DETID(DELID(BRADET))
      ELSE
        DELDET=0
      END IF

      RETURN
      END

