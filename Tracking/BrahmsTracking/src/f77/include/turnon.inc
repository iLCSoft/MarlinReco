*  
* $Id: turnon.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:36  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:29  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:48  hvogt
* Brahms V308
*
*  
      INTEGER ITURNON(20), IKINOPT
      INTEGER ILAT, JLAT
      INTEGER IBCAL, IECAL
      LOGICAL LVXDP,LVXDC,LITC,LFTD,LSIT,LTPC,LBCAL,LECAP,LCOIL
      LOGICAL LYOKE,LLAT,LMASK,LFCH, LBMP, LBMPE
      LOGICAL LQUAD, LHALL, LBMPL
      INTEGER CCDTYPE
      COMMON/TURNON/ITURNON, IKINOPT,
     &              LVXDP,LVXDC,LITC,LFTD,LSIT,LTPC,LBCAL,LECAP,LCOIL,
     &              LYOKE,LLAT,LBMP,LBMPE,LMASK,LFCH,LQUAD,LHALL,LBMPL,
     &              CCDTYPE, ILAT, JLAT, IBCAL, IECAL
      SAVE /TURNON/

