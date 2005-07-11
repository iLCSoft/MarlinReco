*  
* $Id: fkexts.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
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
*  description of inter- and extrapolation surfaces
*
       integer   nexs,nexhpc
       real*4    rzsurf,zrmin,zrmax
       COMMON /FKEXTS/ NEXS,RZSURF(NEXSMX),ZRMIN(NEXSMX),ZRMAX(NEXSMX),
     &                 NEXHPC
*
       CHARACTER*4 TYPX
       COMMON /FKEXTY/ TYPX(NEXSMX)
