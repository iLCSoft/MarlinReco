*  
* $Id: fkcons.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
       real*4     prtot,prsmo
       real*4     raddeg,chi2mx
       integer    lhist,ibafo
       PARAMETER (RADDEG=57.2958,CHI2MX=10.,LHIST=37000)
       COMMON /FKCONS/PRTOT(0:1),PRSMO(0:1),IBAFO
