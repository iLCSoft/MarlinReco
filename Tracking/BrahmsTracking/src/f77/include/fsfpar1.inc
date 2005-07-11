*  
* $Id: fsfpar1.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      integer    spacesize,temax,tsmax
      parameter (spacesize=(1-IFULLBKG)*250000+IFULLBKG*4000000)
      parameter (temax=(1-IFULLBKG)*1000+IFULLBKG*50000)   !max number of TE per detector we can handle.
      parameter (tsmax=(1-IFULLBKG)*2000+IFULLBKG*100000)  !max number of TS we can handle.
c      parameter (spacesize=250000)
c      parameter (temax=1000)   !max number of TE per detector we can handle.
c      parameter (tsmax=2000)   !max number of TS we can handle.


