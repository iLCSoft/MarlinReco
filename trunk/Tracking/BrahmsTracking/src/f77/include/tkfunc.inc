*  
* $Id: tkfunc.inc,v 1.3 2005-07-21 16:37:00 aplin Exp $
*  
* $Log: not supported by cvs2svn $
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
* Revision 1.2  2004/07/22 14:41:07  aplin
* moved specification statements to tkfunc1.inc
*
* Revision 1.1.1.1  2003/05/23 13:17:48  hvogt
* Brahms V308
*
*  
*#include "brlink.inc"
*      REAL RTKHIT, TKMCTR
*      INTEGER ITKHIT, ITKMCT
*      INTEGER IOFFSET, IHITOFFSET
*     TRACK TE variables
*      REAL RTE
*      INTEGER ITE, ITEDAT

*     define the statement functions which replace the
*     static arrays
*     hits array

cSJA	use cpp fortan interface STFREADTKHITS

      RTKHIT(ATTRIBUTE,AHIT) = RSTFREADTKHITS(ATTRIBUTE,AHIT)	
      ITKHIT(ATTRIBUTE,AHIT) = ISTFREADTKHITS(ATTRIBUTE,AHIT)	

	
c      RTKHIT(IOFFSET,IHITOFFSET) =
c     &  Q(LBRHITS1 + MXHT*(IHITOFFSET-1)+IOFFSET)
c      ITKHIT(IOFFSET,IHITOFFSET) =
c     &  IQ(LBRHITS1+ MXHT*(IHITOFFSET-1)+IOFFSET)

*     MC tree array

c      TKMCTR(IOFFSET,IHITOFFSET) =
c     &  Q(LBRMC1 + 11*(IHITOFFSET-1)+IOFFSET)
c      ITKMCT(IOFFSET,IHITOFFSET) =
c     &  IQ(LBRMC1 + 11*(IHITOFFSET-1)+IOFFSET)

cSJA	use cpp fortan interface STFREADTKTE

      RTE(ATTRIBUTE,ATE) = RSTFREADTKTE(ATTRIBUTE,ATE)
      ITE(ATTRIBUTE,ATE) = ISTFREADTKTE(ATTRIBUTE,ATE)

      ITEDAT(ATTRIBUTE,ATE) = STFREADTKITEDAT(ATTRIBUTE,ATE)	 

*     TE bank
c      RTE(IOFFSET,IHITOFFSET) =
c     &  Q(LBRTE+MXTE*(IHITOFFSET-1)+IOFFSET)
c      ITE(IOFFSET,IHITOFFSET) =
c     &  IQ(LBRTE+MXTE*(IHITOFFSET-1)+IOFFSET)
c      ITEDAT(IOFFSET,IHITOFFSET) =
c     &  IQ(LBRT2+5*(IHITOFFSET-1)+IOFFSET)


