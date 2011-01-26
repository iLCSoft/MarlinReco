*  
* $Id: tkcnst1.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:35  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:29  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:48  hvogt
* Brahms V308
*
*  
      INTEGER NHMAX,MXEVTE
* --- maximum number of hits per event
      PARAMETER (NHMAX=(1-IFULLBKG)*25000+IFULLBKG*250000)
c      PARAMETER (NHMAX=20000)
* --- maximum number of TEs per event
      PARAMETER (MXEVTE=(1-IFULLBKG)*3000+IFULLBKG*150000)
c      PARAMETER (MXEVTE=3000)


