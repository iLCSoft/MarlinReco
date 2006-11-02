*  
* $Id: fkddes1.inc,v 1.1 2006-11-02 12:06:20 rasp Exp $
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
       INTEGER NCONMAT

       REAL R1CONMAT(NCMAMX),R2CONMAT(NCMAMX)
       REAL Z1CONMAT(NCMAMX),Z2CONMAT(NCMAMX)	
       REAL XRL1CON(NCMAMX),XRL2CON(NCMAMX)
       REAL XEL1CON(NCMAMX),XEL2CON(NCMAMX)
	

       COMMON /FKDDES1/ NCONMAT,
     & Z1CONMAT,Z2CONMAT,R1CONMAT,R2CONMAT,
     & XRL1CON,XRL2CON,XEL1CON,XEL2CON

