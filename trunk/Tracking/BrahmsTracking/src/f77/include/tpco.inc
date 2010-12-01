*  
* $Id: tpco.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      INTEGER NTPCO,TPCO_TN
      REAL TPCO_XV,TPCO_YV,TPCO_ZV,TPCO_RV,TPCO_PH
      REAL TPCO_SR,TPCO_SZ
      COMMON/TPCO/NTPCO,TPCO_TN(NTPCOMX),TPCO_XV(NTPCOMX),
     & TPCO_YV(NTPCOMX),TPCO_ZV(NTPCOMX),TPCO_RV(NTPCOMX),
     & TPCO_PH(NTPCOMX),TPCO_SR(NTPCOMX),TPCO_SZ(NTPCOMX)

