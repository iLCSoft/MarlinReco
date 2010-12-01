*  
* $Id: fkddes2.inc,v 1.1 2007-11-23 08:50:56 rasp Exp $
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
*   description of material surfaces (cylinders and planes)
*
       INTEGER NPLMAT

       REAL XPLMAT(NPLMAMX),YPLMAT(NPLMAMX)
       REAL ZPLMAT(NPLMAMX)
       REAL WIDPLMAT(NPLMAMX),LENPLMAT(NPLMAMX)
       REAL PHIPLMAT(NPLMAMX)
       REAL XRLPL(NPLMAMX),XELOSPL(NPLMAMX)

       COMMON /FKDDES2/ NPLMAT,XPLMAT,YPLMAT,ZPLMAT,
     & WIDPLMAT,LENPLMAT,PHIPLMAT,
     & XRLPL,XELOSPL

