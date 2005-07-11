*  
* $Id: tkfunc.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
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

      RTKHIT(ATTRIBUTE,AHIT) = STFREADTKHITS(ATTRIBUTE,AHIT)	

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

      RTE(IOFFSET,IHITOFFSET) = STFREADTKTE(ATTRIBUTE,ATE)

*     TE bank
c      RTE(IOFFSET,IHITOFFSET) =
c     &  Q(LBRTE+MXTE*(IHITOFFSET-1)+IOFFSET)
c      ITE(IOFFSET,IHITOFFSET) =
c     &  IQ(LBRTE+MXTE*(IHITOFFSET-1)+IOFFSET)
c      ITEDAT(IOFFSET,IHITOFFSET) =
c     &  IQ(LBRT2+5*(IHITOFFSET-1)+IOFFSET)


