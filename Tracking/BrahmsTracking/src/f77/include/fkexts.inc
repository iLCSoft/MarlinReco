*  
* $Id: fkexts.inc,v 1.3 2008-11-27 17:08:22 engels Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.2  2006/11/02 11:49:01  rasp
*
* Modified files adopted to FullLDCTracking processor
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
* Revision 1.1.1.1  2003/05/23 13:17:49  hvogt
* Brahms V308
*
*  
*  description of inter- and extrapolation surfaces
*
       integer   nexs,nexhpc,itexts
       real*4    rzsurf,zrmin,zrmax
       COMMON /FKEXTS/ NEXS,RZSURF(NEXSMX),ZRMIN(NEXSMX),ZRMAX(NEXSMX),
     &                 ITEXTS(NEXSMX),NEXHPC
*
       CHARACTER*4 TYPX
       COMMON /FKEXTY/ TYPX(NEXSMX)

