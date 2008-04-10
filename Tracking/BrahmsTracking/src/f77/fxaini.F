*  
* $Id: fxaini.F,v 1.2 2008-04-10 21:53:32 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/07/11 15:15:48  aplin
* Initial version
*
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
       SUBROUTINE FXAINI
******************************************************************
#include "include/bkgmode.inc"
#include "include/tkconst.inc"
#include "include/fxapar.inc"
       CALL VZERO(ISTAT,1+NWSTAT)
       CALL VZERO(KSUMER,NECODE)
       CALL VZERO(KSUMWRN,NWCODE)
C       WRITE(6,1000)
 1000  FORMAT(' Ambiguity Processor initialised. Version ',
     +40HDELAMBI    0.5  991123  11.30                   ,
     +/)
ccc       CALL fxsolvversion
       RETURN
       END

