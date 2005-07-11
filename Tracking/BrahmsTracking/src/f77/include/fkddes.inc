*  
* $Id: fkddes.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
*   description of material surfaces (cylinders and planes)
*
       INTEGER NCMAT,NPMAT

       REAL RCMAT(NCMAMX),ZCMAX(NCMAMX),XRLC(NCMAMX)
       REAL ZPMAT(NPMAMX),RPMIN(NPMAMX),RPMAX(NPMAMX),XRLP(NPMAMX)
       REAL XELOSC(NCMAMX),XELOSP(NPMAMX),ZCMIN(NCMAMX)

       COMMON /FKDDES/ NCMAT,RCMAT,ZCMAX,XRLC,
     & NPMAT,ZPMAT,RPMIN,RPMAX,XRLP,
     & XELOSC,XELOSP,ZCMIN

