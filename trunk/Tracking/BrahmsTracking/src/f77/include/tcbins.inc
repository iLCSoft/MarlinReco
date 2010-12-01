*  
* $Id: tcbins.inc,v 1.2 2008-06-26 09:35:56 aplin Exp $
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

      INTEGER IBNPTR,INDBIN,IXWRAP,LNBBIG,LNBSML,LMXBIN
      INTEGER NMXWRP,LNWRAP,JNDBCO,JNDBFW,JNDBBW,LNDBIN
      REAL TCBINW,TCDELT,TC2DEL,TCOFFS
      INTEGER IBLAST
      PARAMETER (LNBBIG=10,LNBSML=1,LMXBIN=LNBBIG*LNBSML)
      PARAMETER (NMXWRP=LNBSML/2,LNWRAP=LMXBIN+4*NMXWRP+2)
      PARAMETER (JNDBCO=1,JNDBFW=2,JNDBBW=3,LNDBIN=3)
      COMMON/TCBINS/ IBNPTR(LMXBIN,LTPDRO),INDBIN,IXWRAP(LNWRAP),
     &               TCBINW,TCDELT,IBLAST,TC2DEL,TCOFFS

